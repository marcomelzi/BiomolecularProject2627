# BiomolecularProject2627
run example: namd3 +p8 prot3.conf > prot3.log <br>
check example for creation of MUT environoment : python3 namd_mutations.py AA95M --build-only <br>
run example for MUT MD : python3 namd_mutations.py AA95M NA52M
python3 compare_q_per_residue.py in order to define the q_value change
1. Scelta degli iperparametri di simulazione: temperatura (510K), dimensioni della box (dimensioni: 43,2 × 43,9 × 51,1; centro: 9,9; 26,5; 19,6), durata della simulazione (1ns), con l’obiettivo di trovare un valore di q sufficientemente basso, così da avere margine di miglioramento: non era necessario arrivare a convergenza ma trovare un benchmark facilmente replicabile in tempi umani 
2. Simulazioni prelimiari del wild type, l’output è q_res, rmsf_core
3. Foldx su tutte le mutazioni (19x106), l’output è DeltaDeltaG
4. Excel con q_res, rmsf_core, DeltaDeltaG: le mutazioni worth it sono quelle con DeltaDeltaG molto negativo, q_res non troppo alto, e numero di contatti non troppo basso (perché numero contatti basso rende noiosi la misurazione di q_res)
5. Effettiva automatizzazione delle mutazioni (run : caffeinate python3 namd_mutations.py AA95M NA52M)
