#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simula mutanti puntiformi di 1EW4: FoldX -> sistema NAMD -> repliche -> analisi.

Uso (dalla cartella del progetto, accanto a 1EW4.pdb e sample.conf):
    python3 namd_mutations.py AA95M NA52M
    python3 namd_mutations.py AA95M --build-only     # solo FoldX + costruzione sistema
    python3 namd_mutations.py WT                     # CONTROLLO WT con la stessa pipeline

Il token WT lancia il wild type attraverso la STESSA pipeline dei mutanti: struttura
WT_1EW4_1.pdb prodotta da FoldX (BuildModel) in una mutazione gia' costruita, stesso
build_system.tcl, stessa cella misurata, stesso conf, stessi seed, stessa analisi.
Serve ad avere un confronto mutante/WT senza differenze di protocollo. Il risultato
va in mutations/WT/. La struttura WT_ si prende da mutations/<SIGLA>/foldx/ (la prima
trovata, oppure quella indicata con --wt-from SIGLA).

Le sigle sono nel formato FoldX: aa WT + catena + residuo + aa mutante (es. AA95M).

Per ogni mutazione crea:
    mutations/AA95M/
        foldx/       BuildModel (PDB mutato)
        system/      psfgen + solvatazione + ionizzazione (ionized.psf/pdb, box.txt)
        replicas/    run00..run04: NAMD + analisi (rmsd.dat, rmsfCA.dat,
                     q_value.txt, q_split.txt)
        info.txt     parametri usati (utile per il report)

Ordine di esecuzione:
    1) controllo di tutte le sigle (prima di lanciare qualunque cosa)
    2) FoldX + costruzione del sistema per TUTTE le mutazioni (fasi veloci)
    3) MD: una mutazione alla volta, con le repliche in parallelo (come il WT)
Se rilanci, le fasi gia' completate vengono saltate; una replica e' considerata
completata se esiste q_split.txt.
"""

import argparse
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# ----------------------------- CONFIGURAZIONE -----------------------------
WORK = Path(__file__).resolve().parent          # cartella del progetto
PDB_WT = WORK / "1EW4.pdb"
CONF = "sample.conf"
PAR = "par_all27_prot_lipid.inp"
BUILD_TCL = WORK / "build_system.tcl"
MUTATIONS_TXT = WORK / "mutations.txt"
OUTROOT = WORK / "mutations"
TOPOLOGY = Path.home() / "vmd-2.0.1a1/plugins/noarch/tcl/membrane1.2/top_all27_prot_lipid.inp"
FOLDX_MOLECULES = Path.home() / "foldx" / "molecules"   # linkata se esiste (innocua)

PADDING = 5.0           # A, come il WT (misurato su ionized.pdb)
CONC = 0.15             # mol/L di sale, come il WT
SEGNAME = "PROT"

# Identico al WT: stessa temperatura e stessi seed (1000 + indice replica)
TEMPERATURE_OLD = "set temperature   310"       # deve comparire ESATTAMENTE nel conf
TEMPERATURE_NEW = "set temperature 510"
SEED0 = 1000
NRUNS = 5
STEPS_DEFAULT = 500000                          # usato solo se 'run N' non e' nel conf

FOLDX, VMD, NAMD = "foldx", "vmd", "namd3"

# File linkati da WORK in ogni cartella runNN (ionized.psf/pdb arrivano da system/)
INPUT_LINKS = [PAR, "analysisMD.tcl", "nativecontacts.tcl", "nativecontacts_split.tcl"]

# PUNTO DI AGGANCIO per l'analisi per residuo (da definire con lo script del compagno):
# elenca qui gli script .tcl (presenti in WORK) da lanciare con VMD dopo run_split.tcl.
# Vengono linkati nelle cartelle runNN e lanciati in ordine.
EXTRA_TCL = ["q_value_per_residue.tcl"]
# ---------------------------------------------------------------------------

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "HSE": "H", "HSD": "H", "HSP": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA_VALIDI = set("ACDEFGHIKLMNPQRSTVWY")
RUN_SPLIT_TCL = (
    "# generato da namd_mutations.py: Q per classe di separazione (percorsi relativi)\n"
    "mol new ionized.psf type psf waitfor all\n"
    "mol addfile prot3.dcd type dcd waitfor all\n"
    "source nativecontacts_split.tcl\n"
    "exit\n"
)
OUTPUT_ATTESI = ["rmsd.dat", "rmsfCA.dat", "q_value.txt", "q_split.txt", "q_per_residue.txt"]


class PipelineError(Exception):
    pass


# ------------------------------ sigle e PDB -------------------------------
def parse_sigla(s):
    s = s.strip().upper()
    m = re.fullmatch(r"([A-Z])([A-Z])(\d+)([A-Z])", s)
    if not m:
        raise PipelineError(f"'{s}': formato non valido (atteso WT+catena+residuo+MUT, es. AA95M)")
    wt, chain, pos, mut = m.group(1), m.group(2), int(m.group(3)), m.group(4)
    for aa in (wt, mut):
        if aa not in AA_VALIDI:
            raise PipelineError(f"'{s}': '{aa}' non e' un aminoacido valido")
    return {"sigla": s, "wt": wt, "chain": chain, "pos": pos, "mut": mut}


def residui_pdb(pdb):
    """{(catena, numero): lettera} leggendo i CA dei record ATOM (prima altloc)."""
    res = {}
    for r in Path(pdb).read_text(errors="replace").splitlines():
        if r[:6] != "ATOM  " or r[12:16].strip() != "CA":
            continue
        if r[16] not in (" ", "A"):
            continue
        try:
            key = (r[21], int(r[22:26]))
        except ValueError:
            continue
        res.setdefault(key, AA3TO1.get(r[17:20].strip(), "X"))
    return res


def valida(m, wt_res):
    key = (m["chain"], m["pos"])
    if key not in wt_res:
        pos = sorted(p for c, p in wt_res if c == m["chain"])
        rng = f"{pos[0]}-{pos[-1]}" if pos else "nessun residuo"
        raise PipelineError(f"{m['sigla']}: residuo {m['chain']}{m['pos']} non trovato in "
                            f"{PDB_WT.name} (catena {m['chain']}: {rng})")
    if wt_res[key] != m["wt"]:
        raise PipelineError(f"{m['sigla']}: in {PDB_WT.name} il residuo {m['chain']}{m['pos']} "
                            f"e' {wt_res[key]}, non {m['wt']}")
    if m["wt"] == m["mut"]:
        raise PipelineError(f"{m['sigla']}: aminoacido WT e mutante coincidono")


def carica_tag():
    tags = {}
    if MUTATIONS_TXT.exists():
        for r in MUTATIONS_TXT.read_text().splitlines():
            parti = r.split()
            if parti and not parti[0].startswith("#"):
                tags[parti[0].upper()] = parti[1] if len(parti) > 1 else ""
    return tags


# ------------------------------ esecuzione -------------------------------
def esegui(cmd, cwd, log):
    with open(log, "w") as f:
        p = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                           stdin=subprocess.DEVNULL)
    if p.returncode != 0:
        coda = "\n".join(Path(log).read_text(errors="replace").splitlines()[-15:])
        raise PipelineError(f"comando fallito (exit {p.returncode}): {' '.join(map(str, cmd))}\n"
                            f"--- ultime righe di {log} ---\n{coda}")


def fase_foldx(m, mdir):
    d = mdir / "foldx"
    d.mkdir(parents=True, exist_ok=True)

    def trova():
        atteso = d / f"{PDB_WT.stem}_1.pdb"
        if atteso.exists():
            return atteso
        altri = sorted(p for p in d.glob(f"{PDB_WT.stem}_*.pdb") if not p.name.startswith("WT_"))
        return altri[0] if altri else None

    mutante = trova()
    if mutante is None:
        shutil.copy(PDB_WT, d / PDB_WT.name)
        (d / "individual_list.txt").write_text(f"{m['sigla']};\n")
        if FOLDX_MOLECULES.is_dir() and not (d / "molecules").exists():
            os.symlink(FOLDX_MOLECULES, d / "molecules")
        esegui([FOLDX, "--command=BuildModel", f"--pdb={PDB_WT.name}",
                "--mutant-file=individual_list.txt", "--numberOfRuns=1"], d, d / "foldx.log")
        mutante = trova()
        if mutante is None:
            raise PipelineError(f"{m['sigla']}: FoldX non ha prodotto {PDB_WT.stem}_1.pdb "
                                f"(vedi {d / 'foldx.log'})")
    # controllo che la mutazione sia davvero stata applicata
    res = residui_pdb(mutante).get((m["chain"], m["pos"]))
    if res != m["mut"]:
        raise PipelineError(f"{m['sigla']}: in {mutante.name} il residuo {m['chain']}{m['pos']} "
                            f"e' {res}, atteso {m['mut']}")
    return mutante


def fase_wt_controllo(mdir, wt_from):
    """Struttura WT prodotta da FoldX (WT_1EW4_1.pdb) copiata in mutations/WT/foldx/."""
    d = mdir / "foldx"
    d.mkdir(parents=True, exist_ok=True)
    nome = f"WT_{PDB_WT.stem}_1.pdb"
    if (d / nome).exists():
        return d / nome
    if wt_from:
        cand = [OUTROOT / wt_from / "foldx" / nome]
    else:
        cand = sorted(p for p in OUTROOT.glob(f"*/foldx/{nome}") if p.parent.parent.name != "WT")
    cand = [p for p in cand if p.exists()]
    if not cand:
        raise PipelineError(
            f"non trovo {nome} in mutations/{wt_from or '<SIGLA>'}/foldx/. Il controllo WT usa la "
            f"struttura che FoldX produce insieme a un mutante: lancia prima almeno una "
            f"mutazione (anche con --build-only)")
    shutil.copy(cand[0], d / nome)
    (d / "origine.txt").write_text(f"{nome} copiato da {cand[0]}\n")
    print(f"[WT] struttura di controllo: {cand[0]}", flush=True)
    return d / nome


def fase_sistema(m, mdir, mutante):
    sd = mdir / "system"
    sd.mkdir(parents=True, exist_ok=True)
    attesi = [sd / n for n in ("ionized.psf", "ionized.pdb", "box.txt")]
    if not all(p.exists() for p in attesi):
        esegui([VMD, "-dispdev", "text", "-eofexit", "-e", str(BUILD_TCL), "-args",
                str(mutante), str(sd), str(TOPOLOGY), str(PADDING), str(CONC), SEGNAME],
               mdir, sd / "build_system.log")
        mancanti = [p.name for p in attesi if not p.exists()]
        if mancanti:
            raise PipelineError(f"{m['sigla']}: build_system.tcl non ha prodotto {mancanti} "
                                f"(vedi {sd / 'build_system.log'})")
    return sd


def conf_mutante(sd):
    """sample.conf con T=510 K e cella periodica presa da box.txt (non quella a mano)."""
    testo = (WORK / CONF).read_text()
    if TEMPERATURE_OLD not in testo:
        raise PipelineError(f"in {CONF} non trovo '{TEMPERATURE_OLD}' (spaziatura esatta): "
                            f"la temperatura non verrebbe cambiata")
    testo = testo.replace(TEMPERATURE_OLD, TEMPERATURE_NEW)
    box = {r.split()[0]: r for r in (sd / "box.txt").read_text().splitlines() if r.strip()}
    for chiave in ("cellBasisVector1", "cellBasisVector2", "cellBasisVector3", "cellOrigin"):
        if chiave not in box:
            raise PipelineError(f"{sd / 'box.txt'}: manca {chiave}")
        riga = box[chiave]
        testo, n = re.subn(rf"^[ \t]*{chiave}\b.*$", lambda _: riga, testo, flags=re.M)
        if n != 1:
            raise PipelineError(f"in {CONF} trovo {n} righe '{chiave}' invece di 1")
    return testo


def scrivi_info(m, mdir, sd, tag, nruns):
    box = (sd / "box.txt").read_text().strip()
    if m.get("controllo"):
        (mdir / "info.txt").write_text(
            f"CONTROLLO WT con la stessa pipeline dei mutanti\n"
            f"creato il      : {datetime.now():%Y-%m-%d %H:%M}\n"
            f"struttura      : WT_{PDB_WT.stem}_1.pdb di FoldX (vedi foldx/origine.txt)\n"
            f"sistema        : padding {PADDING} A, sale {CONC} M, topologia {TOPOLOGY}\n"
            f"MD             : {nruns} repliche, seed {SEED0}-{SEED0 + nruns - 1}, T=510 K\n"
            f"cella periodica (da box.txt):\n{box}\n")
        return
    (mdir / "info.txt").write_text(
        f"mutazione      : {m['sigla']}  ({m['wt']} -> {m['mut']}, catena {m['chain']}, "
        f"residuo {m['pos']})\n"
        f"tag mutations.txt: {tag if tag is not None else 'non presente'}\n"
        f"creato il      : {datetime.now():%Y-%m-%d %H:%M}\n"
        f"FoldX          : BuildModel diretto su {PDB_WT.name}, numberOfRuns=1\n"
        f"sistema        : padding {PADDING} A, sale {CONC} M, topologia {TOPOLOGY}\n"
        f"MD             : {nruns} repliche, seed {SEED0}-{SEED0 + nruns - 1}, T=510 K\n"
        f"cella periodica (da box.txt, ricalcolata per questo sistema):\n{box}\n\n"
        f"NOTA: il WT e' stato simulato con la cella scritta a mano in {CONF} "
        f"(40x41x50, origine 20 20 25),\n"
        f"      diversa dalla cella misurata sul sistema. Da dichiarare nel confronto WT/mutante.\n")


# ------------------------------- MD -------------------------------------
def passo(log):
    """Ultimo step letto dalle righe ENERGY: del log (serve solo per la percentuale)."""
    if not log.exists():
        return 0
    ultimo = 0
    for r in log.read_text(errors="replace").splitlines():
        if r.startswith("ENERGY:"):
            try:
                ultimo = int(r.split()[1])
            except (IndexError, ValueError):
                pass
    return ultimo


def prepara_run(i, mdir, sd, conf):
    d = mdir / "replicas" / f"run{i:02d}"
    d.mkdir(parents=True, exist_ok=True)
    sorgenti = {"ionized.psf": sd / "ionized.psf", "ionized.pdb": sd / "ionized.pdb"}
    for nome in INPUT_LINKS + EXTRA_TCL:
        sorgenti[nome] = WORK / nome
    for nome, src in sorgenti.items():
        dst = d / nome
        if not (dst.is_symlink() or dst.exists()):
            os.symlink(src, dst)
    (d / "run.conf").write_text(f"seed {SEED0 + i}\n" + conf)
    (d / "run_split.tcl").write_text(RUN_SPLIT_TCL)
    return d


def comando_run(cpus):
    passi = [f"{NAMD} +p{cpus} run.conf",
             f"{VMD} -dispdev text -eofexit -e analysisMD.tcl",
             f"{VMD} -dispdev text -eofexit -e run_split.tcl"]
    passi += [f"{VMD} -dispdev text -eofexit -e {t}" for t in EXTRA_TCL]
    return " && ".join(passi)


def fase_md(m, mdir, sd, nruns, poll):
    conf = conf_mutante(sd)
    steps = int(next(iter(re.findall(r"^\s*run\s+(\d+)", conf, flags=re.M)), STEPS_DEFAULT))
    runs = [prepara_run(i, mdir, sd, conf) for i in range(nruns)]
    da_fare = [d for d in runs if not (d / "q_split.txt").exists()]
    if not da_fare:
        print(f"[{m['sigla']}] repliche gia' complete: salto l'MD", flush=True)
        return True
    ncpus = os.cpu_count() or 1
    cpus = max(1, ncpus // len(da_fare))
    if ncpus < len(da_fare):
        print(f"ATTENZIONE: {ncpus} core per {len(da_fare)} repliche (1 core ciascuna)", flush=True)
    print(f"[{m['sigla']}] lancio {len(da_fare)} repliche, {cpus} core ciascuna", flush=True)

    files = [open(d / "run.log", "w") for d in da_fare]
    procs = [subprocess.Popen(["bash", "-c", comando_run(cpus)], cwd=d, stdout=f,
                              stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                              start_new_session=True)
             for d, f in zip(da_fare, files)]
    t0 = time.time()
    try:
        while any(p.poll() is None for p in procs):
            time.sleep(poll)
            fatti = [passo(d / "run.log") for d in da_fare]
            print(f"[{m['sigla']}] {(time.time() - t0) / 60:5.0f} min  "
                  f"{sum(fatti) / len(da_fare) / steps:5.1%}  "
                  + " ".join(f"{s // 1000}k" for s in fatti), flush=True)
    except KeyboardInterrupt:
        print("\nInterrotto: termino le repliche in corso...", flush=True)
        for p in procs:
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGTERM)
        for p in procs:
            p.wait()
        raise
    finally:
        for f in files:
            f.close()

    ok = True
    for d, p in zip(da_fare, procs):
        mancanti = [n for n in OUTPUT_ATTESI if not (d / n).exists()]
        stato = "ok" if p.returncode == 0 and not mancanti else "PROBLEMA"
        extra = f"  (mancano: {', '.join(mancanti)})" if mancanti else ""
        print(f"[{m['sigla']}] {d.name}: exit {p.returncode}  {stato}{extra}", flush=True)
        ok = ok and stato == "ok"
    return ok


# ------------------------------- main -----------------------------------
def controlla_ambiente(build_only):
    mancanti = [str(p) for p in (PDB_WT, WORK / CONF, BUILD_TCL, TOPOLOGY) if not p.exists()]
    if mancanti:
        raise PipelineError("file non trovati: " + ", ".join(mancanti))
    if not build_only:
        for nome in INPUT_LINKS + EXTRA_TCL:
            if not (WORK / nome).exists():
                raise PipelineError(f"file non trovato: {WORK / nome}")
    programmi = [FOLDX, VMD] + ([] if build_only else [NAMD])
    assenti = [p for p in programmi if shutil.which(p) is None]
    if assenti:
        raise PipelineError("programmi non trovati nel PATH: " + ", ".join(assenti))


def main():
    ap = argparse.ArgumentParser(description="Pipeline MD per mutanti puntiformi di 1EW4")
    ap.add_argument("mutazioni", nargs="+", help="sigle FoldX (es. AA95M NA52M) e/o il token WT per il controllo")
    ap.add_argument("--wt-from", default=None, metavar="SIGLA",
                    help="per il token WT: mutazione da cui prendere WT_1EW4_1.pdb (default: la prima trovata)")
    ap.add_argument("--nruns", type=int, default=NRUNS, help=f"repliche per mutazione ({NRUNS})")
    ap.add_argument("--build-only", action="store_true", help="solo FoldX + sistema, niente MD")
    ap.add_argument("--poll", type=int, default=60, help="secondi tra due stampe di avanzamento")
    ap.add_argument("--topology", type=Path, default=None, help="file di topologia CHARMM (.inp)")
    args = ap.parse_args()

    global TOPOLOGY
    if args.topology:
        TOPOLOGY = args.topology.resolve()

    try:
        controlla_ambiente(args.build_only)

        # 1) controllo di tutte le sigle prima di fare qualunque cosa
        wt_res = residui_pdb(PDB_WT)
        tags = carica_tag()
        muts, viste = [], set()
        for s in args.mutazioni:
            if s.strip().upper() == "WT":
                m = {"sigla": "WT", "controllo": True}
            else:
                m = parse_sigla(s)
                valida(m, wt_res)
            if m["sigla"] in viste:
                raise PipelineError(f"{m['sigla']}: indicata piu' volte")
            viste.add(m["sigla"])
            muts.append(m)
        print("Mutazioni: " + ", ".join(
            f"{m['sigla']} (tag {tags.get(m['sigla'], '-')})" for m in muts), flush=True)

        # 2) FoldX + sistema per tutte le mutazioni (fasi veloci)
        info = {}
        for m in muts:
            mdir = OUTROOT / m["sigla"]
            mdir.mkdir(parents=True, exist_ok=True)
            if m.get("controllo"):
                mutante = fase_wt_controllo(mdir, args.wt_from)
            else:
                print(f"[{m['sigla']}] FoldX BuildModel (saltato se gia' fatto)...", flush=True)
                mutante = fase_foldx(m, mdir)
            print(f"[{m['sigla']}] sistema: padding {PADDING} A, {CONC} M (saltato se gia' fatto)...", flush=True)
            sd = fase_sistema(m, mdir, mutante)
            scrivi_info(m, mdir, sd, tags.get(m["sigla"]), args.nruns)
            conf_mutante(sd)                       # fallisce subito se il conf non e' patchabile
            info[m["sigla"]] = (mdir, sd)
        print("Fasi di preparazione completate.", flush=True)
        if args.build_only:
            return 0

        # 3) MD: una mutazione alla volta
        esiti = {}
        for m in muts:
            mdir, sd = info[m["sigla"]]
            esiti[m["sigla"]] = fase_md(m, mdir, sd, args.nruns, args.poll)
        print("\nRiepilogo:")
        for s, ok in esiti.items():
            print(f"  {s}: {'ok' if ok else 'con problemi (vedi run.log nelle cartelle replicas)'}")
        return 0 if all(esiti.values()) else 2

    except PipelineError as e:
        print(f"\nERRORE: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
