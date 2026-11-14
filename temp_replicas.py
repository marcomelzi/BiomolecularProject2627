# to run it, just type:
# python3 temp_replicas.py
"""Lancia N repliche NAMD identiche tranne il seed, una per core, in sottocartelle separate."""

import os
import subprocess
import time
from pathlib import Path

WORK = Path("/Users/User/Desktop/BiomolecularProject2627")
CONF = "sample.conf"
INPUTS = ["ionized.psf", "ionized.pdb", "par_all27_prot_lipid.inp", "analysisMD.tcl", "nativecontacts.tcl"] #file ricopiati nelle cartellette
OUTDIR = WORK / "temperature_replicas"
NRUNS = os.cpu_count()
STEPS = 100000          # come il 'run' nel conf: serve solo per la percentuale
VMD = "/Applications/VMD2b1.app/Contents/MacOS/startup.command"

def prepara(i):
    #costruire una cartella, restituisce il percorso ad essa assoociato
    d = OUTDIR / f"run{i:02d}"
    d.mkdir(parents=True, exist_ok=True)
    for name in INPUTS:
        if not (d / name).exists():
            os.symlink(WORK / name, d / name)
    (d / "run.conf").write_text((WORK / CONF).read_text().replace("set temperature   310", f"set temperature {330 + i}"))
    # qui setto la temperatura, in particolare: scrive un file nuovo, replicas/run01/run.conf, 
    # il cui contenuto è la riga del seed seguita da tutta la copia del conf originale:
    return d


def passo(log):     #necessario solo per i print
    # d / "run.log" è documento di testo che contiene ciò che normalmente ogni processore stamperebbe
    # a terminale, noi consultiamo quelle che iniziano con la ENERGY
    if not log.exists():
        return 0
    righe_con_energy = [r for r in log.read_text().splitlines() if r.startswith("ENERGY:")]
    if righe_con_energy:
        ultima_riga_con_energy = righe_con_energy[-1]
        seconda_parola = ultima_riga_con_energy.split()[1]
        return int(seconda_parola) #l'iterazione è la seconda parola nella riga che inizia con ENERGY nel run.log
    else:
        return 0



# CODE TO ESECUTE
runs = [prepara(i) for i in range(NRUNS)]
# print(runs) stampa un vettore di percorsi, che possiamo pensare come un vettore di stringhe, del tipo:
# [PosixPath('/Users/User/Desktop/BiomolecularProject2627/replicas/run00'),  ... , PosixPath('.../run07')]

# ogni replica su un core: NAMD con +p1, il sistema le distribuisce sui core liberi
proc = [subprocess.Popen(
            ["bash", "-c", f"namd3 +p1 run.conf && {VMD} -dispdev text -e analysisMD.tcl"], #input che ogni core sta dando
            cwd=d,                                 #working directory di ogni core
            stdout=open(d / "run.log", "w"),       #output file di ogni core
            stderr=subprocess.STDOUT)              #all'atto pratico printa nello stdout 
        for d in runs]

t0 = time.time()
while any(p.poll() is None for p in proc):         
#p.poll()=None se il processo non è finito, riga significa "finchè non hanno finito tutti i core"
    time.sleep(60)
    fatti = [passo(d / "run.log") for d in runs]
    print(f"{(time.time()-t0)/60:5.0f} min  {sum(fatti)/NRUNS/STEPS:5.1%}  "
          + " ".join(f"{s//1000}k" for s in fatti), flush=True)

print("\n" + "\n".join(f"{d.name}: exit {p.returncode}" for d, p in zip(runs, proc)))