#!/usr/bin/env python3
"""Saturation scan FoldX: ogni residuo mutato nei 20 amminoacidi -> Excel ordinato per ddG."""

import multiprocessing as mp
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from openpyxl import Workbook

WORKDIR = Path("/Users/User/Desktop/Università/BML/foldx5_MacSilicon_0")
PDB = "1EW4_Repair.pdb"
FOLDX = "foldx5.1_20261231"
OUT = "1EW4_saturation_scan.xlsx"
NCORES = os.cpu_count()

# Nome della colonna -> sua posizione nel file energies_<n>_<pdb>.txt.
# L'ordine NON e' quello in cui FoldX stampa a schermo: Entropy_sidec sta in 8
# ed energy_torsion in 13. Verificato confrontando una riga del file con lo stdout.
TERMS = [
    ("Total DeltaDeltaG",  0), ("BackHbond",        1), ("SideHbond",         2),
    ("Energy_VdW",         3), ("Electro",          4), ("Energy_SolvP",      5),
    ("Energy_SolvH",       6), ("Energy_vdwclash",  7), ("energy_torsion",   13),
    ("backbone_vdwclash", 14), ("Entropy_sidec",    8), ("Entropy_mainc",     9),
    ("Energy_Ionisation", 20),
]

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def residues():
    """Tutti i residui standard del PDB: (catena, numero, amminoacido)."""
    out = []
    for line in open(WORKDIR / PDB):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in " A":
            aa = AA3.get(line[17:20])
            if aa:
                out.append((line[21], int(line[22:26]), aa))
    return out


def read_energies(path):
    """energies_<n>_<pdb>.txt -> (valori del wild-type, {amminoacido: valori})."""
    wtref, muts = None, {}
    for line in open(path):
        name, *fields = line.split("\t")
        if len(fields) < 21:
            continue
        vals = [float(x) for x in fields]
        tag = name.split("_")[0].upper()
        if tag == "WTREF":
            wtref = vals
        else:
            muts[tag] = vals
    return wtref, muts


def scan(chain, num, wt):
    """PositionScan su una posizione. Una riga per mutante, gia' come ddG."""
    box = Path(tempfile.mkdtemp())
    for name in (PDB, "molecules", "rotabase.txt"):
        os.symlink(WORKDIR / name, box / name)

    subprocess.run(
        [str(WORKDIR / FOLDX), "--command=PositionScan",
         f"--pdb={PDB}", f"--positions={wt}{chain}{num}a"],
        cwd=box, capture_output=True,
    )

    wtref, muts = read_energies(next(box.glob("energies_*.txt")))
    rows = [
        [num, chain, wt, AA3[aa]] + [round(vals[i] - wtref[i], 4) for _, i in TERMS]
        for aa, vals in muts.items() if AA3[aa] != wt
    ]

    shutil.rmtree(box)                      # niente file lasciati sul disco
    print(f"{wt}{chain}{num}: {len(rows)} mutanti", flush=True)
    return rows


def work(my_residues):
    """Un processo: scorre le posizioni che gli sono state assegnate."""
    return [row for r in my_residues for row in scan(*r)]


if __name__ == "__main__":
    res = residues()

    # round-robin: il processo k prende i residui k, k+NCORES, k+2*NCORES, ...
    with mp.Pool(NCORES) as pool:
        parts = pool.map(work, [res[k::NCORES] for k in range(NCORES)])

    rows = [row for part in parts for row in part]
    rows.sort(key=lambda r: r[4])           # ddG crescente: stabilizzanti in cima

    wb = Workbook()
    wb.active.append(["RESIDUE N", "CHAIN", "WILD TYPE RESIDUE", "MUTATION RESIDUE"]
                     + [name for name, _ in TERMS])
    for row in rows:
        wb.active.append(row)
    wb.save(WORKDIR / OUT)

    print(f"\n{len(rows)} mutanti -> {WORKDIR / OUT}")