# beyondcore

Questa repository contiene strumenti e script per l'ottimizzazione di risorse in sistemi modellati tramite Layered Queueing Networks (LQN), con particolare attenzione al trade-off tra costi operativi e prestazioni. Il progetto include codice per la risoluzione di problemi di Mixed-Integer Nonlinear Programming (MINLP), analisi di risultati e generazione di grafici.

## Contenuto della repository

- **lqn_optimization.jl**: Script principale in Julia che implementa la risoluzione del problema di ottimizzazione delle risorse in una LQN tramite JuMP e il solver SCIP. Permette di variare i parametri di carico e di pesi nell'obiettivo, salvando i risultati in file CSV. _Vedi la sezione dedicata sotto._
- **matplot.m**: Script MATLAB che legge i risultati CSV generati da `lqn_optimization.jl` e produce grafici di analisi, sia per singole simulazioni che per confronti tra diversi pesi di performance. I grafici vengono esportati in PDF nella cartella `plots/`.
- **epew_paper_section.tex**: Sezione di articolo scientifico (in LaTeX) che descrive il modello LQN, la formulazione del problema di ottimizzazione e le principali equazioni e vincoli. Utile come riferimento teorico e per la scrittura di pubblicazioni.
- **data/**: Cartella che contiene i file CSV di output prodotti da `lqn_optimization.jl` per diversi valori del parametro PerformanceWeight. Ogni file (es. `variable_load_results_PW0p1.csv`) riporta, per una serie di carichi, i risultati dell'ottimizzazione: numero di utenti, costo totale, throughput, numero totale di core, tempo di risposta, livelli di servizio selezionati e stato della soluzione.
- **plots/**: Cartella che raccoglie i grafici PDF generati da `matplot.m`, sia come subplots per singole simulazioni (es. `variable_load_results_PW0p8_subplot_Cost.pdf`) sia come confronti aggregati tra diversi pesi (es. `comparison_AvgCost.pdf`).

## Descrizione dettagliata: `lqn_optimization.jl`

Questo script Julia implementa un framework per l'ottimizzazione della configurazione di risorse in sistemi modellati tramite Layered Queueing Networks (LQN). Le principali funzionalità sono:

- **Definizione del modello**: Il sistema è rappresentato tramite una matrice di transizione (Jump Matrix), una matrice di routing, insiemi discreti di possibili service rate e costi associati per ciascuna stazione.
- **Formulazione MINLP**: Il problema viene formulato come un Mixed-Integer Nonlinear Programming, dove le variabili decisionali includono il numero di core per stazione e la selezione del service rate. L'obiettivo è minimizzare una funzione pesata tra costo operativo e penalità sulle prestazioni (scostamento dal throughput desiderato).
- **Risoluzione tramite JuMP + SCIP**: Utilizza JuMP come linguaggio di modellazione e SCIP come solver per problemi non lineari misti interi.
- **Analisi su carico variabile**: Lo script esegue una serie di ottimizzazioni variando il carico (numero di utenti) secondo una legge sinusoidale, salvando i risultati in un DataFrame e poi in file CSV.
- **Output**: I risultati includono, per ogni carico, il costo totale, throughput, tempo di risposta, configurazione ottimale di core e service rate, e lo stato della soluzione.

## Utilizzo tipico

1. **Esecuzione dell'ottimizzazione**: Lanciare `lqn_optimization.jl` in Julia per generare i risultati su diversi carichi e pesi di performance.
2. **Analisi e visualizzazione**: Utilizzare `matplot.m` in MATLAB per generare e salvare i grafici a partire dai file CSV prodotti.
3. **Riferimento teorico**: Consultare `epew_paper_section.tex` per dettagli sul modello matematico e la formulazione del problema.

## Note aggiuntive
- I file nella cartella `data/` sono copie dei risultati CSV per diversi pesi di performance.
- I file nella cartella `plots/` sono esclusivamente grafici esportati in PDF e non contengono dati testuali.
- Il file `.DS_Store` nella cartella `plots/` può essere ignorato.