using JuMP
using SCIP
using LinearAlgebra
using CSV         # Aggiunto per scrivere CSV
using DataFrames  # Aggiunto per gestire i dati tabellari

"""
    solve_lqn_optimization(
        M::Int,
        N::Int,
        J::Matrix{Int},
        P_matrix::Matrix{Float64},
        PossibleServiceRates::Dict{Int, Vector{Float64}},
        ServiceRateCosts::Dict{Int, Vector{Float64}},
        MinCores::Vector{Float64},
        CoreCosts::Vector{Float64},
        MaxTotalCores::Union{Float64, Nothing}=nothing,
        TotalQueueLength::Union{Float64, Nothing}=nothing,
        SteadyStateTolerance::Float64=1e-6,
        CostWeight::Float64=0.7,
        PerformanceWeight::Float64=0.3
    )

Solve the LQN optimization problem using JuMP and SCIP as a non-linear mixed integer problem (MINLP)
based on a Markov population process approach.

Parameters:
- M: Number of stations/tasks/entries
- N: Number of transitions in the network
- J: N×M jump matrix describing the evolution of clients within the LQN
       Negative values (-1) indicate source stations, positive values (1) indicate destination stations
- P_matrix: M×M routing matrix (p_ij)
- PossibleServiceRates: Dictionary mapping each station to possible service rates
- ServiceRateCosts: Dictionary mapping each station to costs for each possible service rate
- MinCores: M×1 vector of minimum cores per station
- CoreCosts: M×1 vector of costs per core for each station
- MaxTotalCores: Optional maximum total cores across all stations (excluding station 1)
- TotalQueueLength: Optional total steady-state queue length across all stations
- SteadyStateTolerance: Tolerance for steady-state condition (default: 1e-6)
- CostWeight: Weight for the cost component in the objective function (default: 0.7)
- PerformanceWeight: Weight for the performance component in the objective function (default: 0.3)
"""
function solve_lqn_optimization(
    M::Int,
    N::Int,
    J::Matrix{Int},
    P_matrix::Matrix{Float64},
    PossibleServiceRates::Dict{Int, Vector{Float64}},
    ServiceRateCosts::Dict{Int, Vector{Float64}},
    MinCores::Vector{Float64},
    CoreCosts::Vector{Float64},
    MaxTotalCores::Union{Float64, Nothing}=nothing,
    TotalQueueLength::Union{Float64, Nothing}=nothing,
    SteadyStateTolerance::Float64=1e-6,
    CostWeight::Float64=0.7,
    PerformanceWeight::Float64=0.3
)
    # Create the optimization model with SCIP
    model = Model(SCIP.Optimizer)
    
    # Set SCIP parameters for MINLP
    set_attribute(model, "numerics/feastol", 1e-4)  # Aumento la tolleranza
    set_attribute(model, "numerics/infinity", 1e10)
    set_attribute(model, "limits/time", 300)  # 5 minutes time limit
    
    # Parametri aggiuntivi per migliorare la convergenza dei problemi non lineari
    set_attribute(model, "numerics/epsilon", 1e-6)  # Aumento la tolleranza
    
    # Enable more detailed output from SCIP
    set_attribute(model, "display/verblevel", 5)
    
    # 1. Declare variables
    # Number of cores for each station (continuous)
    # Note: nc[1] is fixed at 1000, not a variable
    @variable(model, nc[i=2:M] >= MinCores[i])
    
    # Binary variables for service rate selection
    # z[1,1] è sempre 1 dato che c'è una sola opzione per la stazione 1
    @variable(model, z[i=2:M, k=1:length(PossibleServiceRates[i])], Bin)
    
    # Effective service rate for each station (continuous)
    # mu[1] è costante
    @variable(model, mu[i=2:M] >= 0)
    
    # Transition rates (T vector)
    @variable(model, T[i=1:N] >= 0)
    
    # Population derivative vector (dX)
    @variable(model, dX[i=1:M])
    
    # Steady-state queue length vector (X)
    @variable(model, X[i=1:M] >= 0)
    
    # Local response time for each station
    #@variable(model, Tl[i=1:M] >= 0)
    
    # Total response time for each station
    #@variable(model, T_total[i=1:M] >= 0)
    
    # Definisci la costante per il tasso di servizio della stazione 1
    mu_1 = PossibleServiceRates[1][1]
    
    # Imposta punti iniziali per le variabili chiave
    for i in 2:M
        # Imposta un punto iniziale ragionevole per il numero di core
        set_start_value(nc[i], max(MinCores[i], 1.0))
        
        # Imposta il punto iniziale per la selezione del tasso di servizio
        # Seleziona il tasso di servizio medio se disponibile, altrimenti il primo
        mid_idx = div(length(PossibleServiceRates[i]), 2) + 1
        for k in 1:length(PossibleServiceRates[i])
            set_start_value(z[i,k], k == mid_idx ? 1.0 : 0.0)
        end
        
        # Imposta un punto iniziale per il tasso di servizio effettivo
        set_start_value(mu[i], PossibleServiceRates[i][mid_idx])
    end
    
    # Imposta punti iniziali per le lunghezze delle code
    if !isnothing(TotalQueueLength)
        for i in 1:M
            set_start_value(X[i], TotalQueueLength / M)
        end
    else
        for i in 1:M
            set_start_value(X[i], 10.0)
        end
    end
    
    # 2. Service rate selection constraints
    # Per le stazioni 2 a M
    for i in 2:M
        # Exactly one service rate must be selected per station
        @constraint(model, sum(z[i,k] for k in 1:length(PossibleServiceRates[i])) == 1)
        
        # Determine the effective service rate based on selection
        @constraint(model, mu[i] == sum(z[i,k] * PossibleServiceRates[i][k] 
                                      for k in 1:length(PossibleServiceRates[i])))
    end
    
    # 3. Transition rate constraints
    # For each transition, the rate is proportional to the service rate and cores of the source station
    for i in 1:N
        # Find which station is the source (has -1 in the jump matrix)
        source_idx = findfirst(j -> J[i,j] == -1, 1:M)
        dest_idx = findfirst(j -> J[i,j] == 1, 1:M)
        if !isnothing(source_idx)
            if source_idx == 1
                # For station 1, use the fixed value of 1000 cores and constant service rate mu_1
                @constraint(model, T[i] == mu_1 * X[source_idx]*P_matrix[source_idx,dest_idx])
            else
                # For other stations, use the variable nc and mu
                @constraint(model, T[i] == mu[source_idx] * nc[source_idx]*P_matrix[source_idx,dest_idx])
            end
        end
    end
    
    # 4. Population dynamics constraints
    # dX = J' * T (derivative of population)
    @constraint(model,dX.==J'*T)
    
    # 5. Steady-state constraints
    # Ensure that the population derivatives are close to zero (steady-state condition)
    for i in 1:M
        @constraint(model, -SteadyStateTolerance <= dX[i] <= SteadyStateTolerance)
    end
    
    # 6. Queue length constraints
    for i in 2:M
            # For other stations, use the variable nc
            @constraint(model, X[i] >= nc[i])
    end
    
    # If TotalQueueLength is specified, add constraint on sum of queue lengths
    if !isnothing(TotalQueueLength)
        @constraint(model, sum(X[i] for i in 1:M) == TotalQueueLength)
    end
    
    # 9. Resource constraints (if MaxTotalCores is specified)
    if !isnothing(MaxTotalCores)
        # MaxTotalCores applies only to stations 2 to M (excluding station 1)
        @constraint(model, sum(nc[i] for i in 2:M) <= MaxTotalCores)
    end
    
    # 10. Define cost expressions symbolically
    # For other stations, use the variable nc
    station_costs = Dict{Int, Any}()
    
    for i in 2:M
        station_costs[i] = sum(z[i,k] * ServiceRateCosts[i][k] for k in 1:length(PossibleServiceRates[i])) * nc[i]
    end
    
    # Total cost is the sum of all station costs
    total_cost = sum(station_costs[i] for i in 2:M)
    
    # 11. Set objective: minimize weighted sum of cost and maximize throughput
    # We use a linear objective function that SCIP can handle
    # We minimize system response time and cost
    @objective(model, Min, 
        (1-PerformanceWeight) * total_cost + PerformanceWeight * (mu_1*X[1]- TotalQueueLength)^2
    )
    
    # 12. Solve the model
    optimize!(model)
    
    # 13. Return results
    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        # Create a complete cores vector including the fixed value for station 1
        # Convert DenseAxisArray to standard array before concatenation
        nc_values = [value(nc[i]) for i in 2:M]
        cores = [1000.0; nc_values]
        
        # Recupera i valori di mu per le stazioni 2 a M
        mu_values = [value(mu[i]) for i in 2:M]
        # Aggiungi il valore costante per la stazione 1
        service_rates = [mu_1; mu_values]
        
        # Calculate performance metrics for reporting
        throughput = value(mu_1*X[1])
        
        # Calculate the actual costs after optimization
        actual_station_costs = Dict{Int, Float64}()
        
        # For station 1
        selected_rate_idx_1 = 0
        actual_station_costs[1] = 0
        
        # For other stations
        for i in 2:M
            selected_rate_idx = [k for k in 1:length(PossibleServiceRates[i]) 
                                if value(z[i,k]) > 0.5][1]
            actual_station_costs[i] = ServiceRateCosts[i][selected_rate_idx] * value(nc[i])
        end
        
        actual_total_cost = sum(actual_station_costs[i] for i in 1:M)
        
        return Dict(
            "status" => status,
            "objective_value" => objective_value(model),
            "cores" => cores,
            "service_rates" => service_rates,
            "transition_rates" => value.(T),
            "population_derivative" => value.(dX),
            "queue_lengths" => value.(X),
            #"local_response_times" => value.(Tl),
           # "total_response_times" => value.(T_total),
            "selected_rates" => Dict(
                i => (i == 1 ? 1 : # Stazione 1 ha sempre il primo tasso
                     [k for k in 1:length(PossibleServiceRates[i]) 
                     if value(z[i,k]) > 0.5][1])
                for i in 1:M
            ),
            "total_cost" => actual_total_cost,
            "station_costs" => actual_station_costs,
            "throughput" => throughput
        )
    else
        return Dict(
            "status" => status,
            "error" => "Optimization did not converge to an optimal solution"
        )
    end
end

# Jump matrix J
# Negative values (-1) indicate source stations, positive values (1) indicate destination stations
J = [
    -1  1  0;
    -1  0  1;
     1 -1  0;
     1  0 -1;
]

# Example data
M = size(J,2)  # Number of stations
N = size(J,1)  # Number of transitions

P_matrix = [
    0.0 0.5 0.5;
    1.0 0.0 0.0;
    1.0 0.0 0.0;
]
PossibleServiceRates = Dict(
    1 => [1.0],
    2 => [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0],  # 8 opzioni, raddoppiando
    3 => [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0]   # 8 opzioni, raddoppiando
)
# Cost for each service rate option for each station
ServiceRateCosts = Dict(
    1 => [1.0],
    2 => [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0], # 8 opzioni, costo = 10^(indice-1)
    3 => [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0]  # 8 opzioni, costo = 10^(indice-1)
)

MinCores = [0.0, 0.0, 0.0]
CoreCosts = [1.0, 1.0, 1.0]  # Cost per core for each station
MaxTotalCores = 100000.0
TotalQueueLength = 1.0  # Total steady-state queue length across all stations
SteadyStateTolerance = 1e-6  # Tolerance for steady-state condition
CostWeight = 0.5  # Bilancio il peso dei costi nell'obiettivo
PerformanceWeight = 0.8  # Bilancio il peso delle prestazioni nell'obiettivo

# Solve the optimization problem
result = solve_lqn_optimization(
    M, N, J, P_matrix, PossibleServiceRates, ServiceRateCosts,
    MinCores, CoreCosts, MaxTotalCores, TotalQueueLength, SteadyStateTolerance,
    CostWeight, PerformanceWeight
) 

println("--- Single Run Result ---")
println(result)
println("-------------------------")

# --- Variable Load Analysis ---
println("\n--- Variable Load Analysis ---")
# Genera valori di carico con andamento sinusoidale da 1 a 800
num_points = 60
min_load = 10.0
max_load = 1000.0
amplitude = (max_load - min_load) / 2.0
midpoint = (max_load + min_load) / 2.0
x_vals = range(-π/2, 3π/2, length=num_points) # Un ciclo completo di sin partendo dal minimo
sin_values = sin.(x_vals)
load_values_float = midpoint .+ amplitude .* sin_values
load_values = round.(load_values_float, digits=1) # Arrotonda a 1 cifra decimale
println("Generated sinusoidal load values: ", load_values)

#load_values = [1,5.0, 10.0, 20.0, 50.0, 100.0,200.0,500.0] # Esempio di carichi variabili
# Inizializza un DataFrame vuoto per i risultati
results_df = DataFrame(Users=Float64[], TotalCost=Float64[], Throughput=Float64[], 
                       TotalCores=Float64[], ResponseTime=Float64[],
                       SelectedRate_Station2=Int[], SelectedRate_Station3=Int[], Status=String[])

for current_load in load_values
    println("\nSolving for TotalQueueLength = ", current_load)
    current_result = solve_lqn_optimization(
        M, N, J, P_matrix, PossibleServiceRates, ServiceRateCosts,
        MinCores, CoreCosts, MaxTotalCores, current_load, SteadyStateTolerance,
        CostWeight, PerformanceWeight
    ) 
    
    status_str = string(current_result["status"])
    if current_result["status"] == MOI.OPTIMAL || current_result["status"] == MOI.LOCALLY_SOLVED
        # Estrai i tassi selezionati per le stazioni 2 e 3
        rate_idx_2 = current_result["selected_rates"][2]
        rate_idx_3 = current_result["selected_rates"][3]
        
        # Calcola il numero totale di core (stazioni 2 e 3)
        total_cores_used = sum(current_result["cores"][2:end])
        
        # Calcola il tempo di risposta
        throughput = current_result["throughput"]
        response_time = throughput > 1e-9 ? current_load / throughput : Inf # Evita divisione per zero
        
        # Arrotonda i valori a 3 cifre decimali per il CSV
        total_cost_rounded = round(current_result["total_cost"], digits=3)
        throughput_rounded = round(throughput, digits=3)
        total_cores_rounded = round(total_cores_used, digits=3)
        response_time_rounded = isinf(response_time) ? Inf : round(response_time, digits=3)

        # Aggiungi una riga al DataFrame
        push!(results_df, (
            Users=current_load,
            TotalCost=total_cost_rounded,
            Throughput=throughput_rounded,
            TotalCores=total_cores_rounded,
            ResponseTime=response_time_rounded,
            SelectedRate_Station2=rate_idx_2,
            SelectedRate_Station3=rate_idx_3,
            Status=status_str
        ))
        println("  Status: ", status_str, " -> Results recorded.")
    else
        println("  Optimization failed for load ", current_load, " with status: ", status_str)
        # Opzionale: registra anche i fallimenti
        push!(results_df, (
            Users=current_load,
            TotalCost=NaN,
            Throughput=NaN,
            TotalCores=NaN,
            ResponseTime=NaN,
            SelectedRate_Station2=0,
            SelectedRate_Station3=0,
            Status=status_str
        ))
    end
end

# Salva il DataFrame in un file CSV
# Formatta il PerformanceWeight per il nome del file (es. 0.8 -> "0p8")
pw_str = replace(string(PerformanceWeight), "." => "p") 
output_csv_file = "variable_load_results_PW$(pw_str).csv"
try
    CSV.write(output_csv_file, results_df)
    println("\nResults saved to '$output_csv_file'")
catch e
    println("\nError saving results to CSV: ", e)
end

println("-----------------------------------") 