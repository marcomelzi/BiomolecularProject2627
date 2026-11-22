#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Confronta il Q per residuo dei mutanti con quello del wild type (WT).

Uso (dalla cartella del progetto):
    python3 compare_q_per_residue.py                # tutte le mutazioni con dati completi
    python3 compare_q_per_residue.py AA95M NA52M    # solo queste
    python3 compare_q_per_residue.py --alpha 0.1
    python3 compare_q_per_residue.py --wt-dir replicas    # forza il WT originale

Dati letti:
    mutante  mutations/<SIGLA>/replicas/run*/q_per_residue.txt
    WT       di default il CONTROLLO mutations/WT/replicas/run*/ (stessa pipeline dei
             mutanti, lanciato con "namd_mutations.py WT"); se non esiste, il WT
             originale replicas/run*/ (con un avviso). Con --wt-dir si sceglie a mano.

Output in mutations/evaluations/:
    <SIGLA>_vs_WTcontrollo.txt     report leggibile con il verdetto (WT = controllo)
    <SIGLA>_vs_WT.txt              idem se il WT e' quello originale
    <SIGLA>_vs_*_per_residue.tsv   tabella per residuo (WT vs mutante)
    summary.tsv                    una riga per (mutazione, WT usato), aggiornato a ogni lancio

Come si decide il verdetto (stabilizzante / destabilizzante / nessun effetto):
  - Per ogni replica si calcola UN Q globale: media dei q_res pesata sui n_nativi
    (equivale al Q globale medio sulla finestra di frame).
  - Si confrontano i Q globali delle repliche mutante contro quelle del WT con un
    test di permutazione esatto sulla differenza delle medie (due code).
  - Le repliche sono l'unita' statistica: i residui NON sono indipendenti, quindi
    la tabella per residuo e' solo descrittiva (si segna dove la differenza supera
    2 errori standard).
"""

import argparse
import itertools
import math
import random
import re
import statistics
import sys
from pathlib import Path

WORK = Path(__file__).resolve().parent
WT_REPLICAS = WORK / "replicas"
MUT_ROOT = WORK / "mutations"
EVAL_DIR = MUT_ROOT / "evaluations"
MUTATIONS_TXT = WORK / "mutations.txt"
FILE_Q = "q_per_residue.txt"
ALPHA = 0.05
INTORNO = 5                 # residui da ciascun lato del sito mutato per la media locale
MAX_PERM_ESATTE = 200000    # oltre, permutazioni casuali (seed fisso)
SIGLA_RE = re.compile(r"[A-Z][A-Z]\d+[A-Z]")


class ConfrontoError(Exception):
    pass


# ------------------------------- lettura ---------------------------------
def leggi_q(path):
    finestra, res = None, {}
    for r in Path(path).read_text(errors="replace").splitlines():
        if r.startswith("#"):
            if "finestra" in r:
                finestra = r.lstrip("# ").split("(")[0].replace("finestra", "", 1).strip()
            continue
        p = r.split()
        if len(p) != 4:
            continue
        q = None if p[3].lower() == "nan" else float(p[3])
        res[int(p[0])] = (p[1], int(p[2]), q)
    if not res:
        raise ConfrontoError(f"{path}: nessun dato letto")
    return {"finestra": finestra, "res": res}


def raccogli(cartella_replicas):
    runs = []
    for d in sorted(Path(cartella_replicas).glob("run*")):
        if (d / FILE_Q).exists():
            runs.append((d.name, leggi_q(d / FILE_Q)))
    return runs


def carica_tag():
    tags = {}
    if MUTATIONS_TXT.exists():
        for r in MUTATIONS_TXT.read_text().splitlines():
            p = r.split()
            if p and not p[0].startswith("#"):
                tags[p[0].upper()] = p[1] if len(p) > 1 else ""
    return tags


# ------------------------------- statistica -------------------------------
def media(x):
    return sum(x) / len(x) if x else math.nan


def sem(x):
    return statistics.stdev(x) / math.sqrt(len(x)) if len(x) > 1 else math.nan


def q_globale(dati):
    """Media dei q_res pesata sui n_nativi (= Q globale medio sulla finestra)."""
    v = [(n, q) for _, n, q in dati["res"].values() if q is not None and n > 0]
    return sum(n * q for n, q in v) / sum(n for n, _ in v)


def test_permutazione(a, b):
    """p-value a due code per la differenza media(b) - media(a)."""
    tutti, na, nb = list(a) + list(b), len(a), len(b)
    obs = abs(media(b) - media(a))
    tot = math.comb(na + nb, nb)
    conta = 0
    if tot <= MAX_PERM_ESATTE:
        for idx in itertools.combinations(range(na + nb), nb):
            s = set(idx)
            gb = [tutti[i] for i in idx]
            ga = [tutti[i] for i in range(na + nb) if i not in s]
            if abs(media(gb) - media(ga)) >= obs - 1e-12:
                conta += 1
        return conta / tot, "esatto", tot
    rng = random.Random(12345)
    n_perm = 100000
    for _ in range(n_perm):
        rng.shuffle(tutti)
        if abs(media(tutti[na:]) - media(tutti[:na])) >= obs - 1e-12:
            conta += 1
    return (conta + 1) / (n_perm + 1), "Monte Carlo", n_perm


# ------------------------------- confronto -------------------------------
def confronta(sigla, wt_runs, mut_runs, tag, alpha, wt_tipo, wt_label):
    if len(wt_runs) < 2 or len(mut_runs) < 2:
        raise ConfrontoError(f"servono almeno 2 repliche per gruppo (WT: {len(wt_runs)}, "
                             f"{sigla}: {len(mut_runs)})")
    finestre = {d["finestra"] for _, d in wt_runs + mut_runs}
    if len(finestre) != 1:
        raise ConfrontoError(f"finestre di frame diverse tra le repliche: {sorted(map(str, finestre))}. "
                             f"Rilancia q_value_per_residue.tcl con la stessa finestra ovunque")
    finestra = finestre.pop()

    tutte = [set(d["res"]) for _, d in wt_runs + mut_runs]
    resids = sorted(set.intersection(*tutte))
    avviso_res = None
    if len(resids) != len(set.union(*tutte)):
        avviso_res = (f"i residui non coincidono tra tutte le repliche: uso i {len(resids)} "
                      f"in comune su {len(set.union(*tutte))}")

    # --- verdetto globale
    gw = [q_globale(d) for _, d in wt_runs]
    gm = [q_globale(d) for _, d in mut_runs]
    delta = media(gm) - media(gw)
    p, metodo, n_perm = test_permutazione(gw, gm)
    if p < alpha:
        esito = "STABILIZZANTE" if delta > 0 else "DESTABILIZZANTE"
    else:
        esito = "NESSUN EFFETTO SIGNIFICATIVO"

    # --- dettaglio per residuo
    pos = int(re.search(r"\d+", sigla).group())
    righe = []
    for r in resids:
        qw = [d["res"][r][2] for _, d in wt_runs if d["res"][r][2] is not None]
        qm = [d["res"][r][2] for _, d in mut_runs if d["res"][r][2] is not None]
        nw = media([d["res"][r][1] for _, d in wt_runs])
        nm = media([d["res"][r][1] for _, d in mut_runs])
        d_r = media(qm) - media(qw) if qw and qm else math.nan
        se = math.sqrt(sem(qw) ** 2 + sem(qm) ** 2) if len(qw) > 1 and len(qm) > 1 else math.nan
        flag = ""
        if not math.isnan(d_r) and not math.isnan(se) and abs(d_r) > 2 * se and abs(d_r) > 0:
            flag = "UP" if d_r > 0 else "DOWN"
        righe.append({"resid": r, "resname": wt_runs[0][1]["res"][r][0], "sito": r == pos,
                      "n_wt": nw, "n_mut": nm, "q_wt": media(qw), "sem_wt": sem(qw),
                      "q_mut": media(qm), "sem_mut": sem(qm), "delta": d_r, "flag": flag})
    return {"wt_tipo": wt_tipo, "wt_label": wt_label,
            "sigla": sigla, "tag": tag, "finestra": finestra, "avviso_res": avviso_res,
            "gw": gw, "gm": gm, "delta": delta, "p": p, "metodo": metodo, "n_perm": n_perm,
            "alpha": alpha, "esito": esito, "righe": righe, "pos": pos,
            "wt_run_names": [n for n, _ in wt_runs], "mut_run_names": [n for n, _ in mut_runs]}


# ------------------------------- output -------------------------------
def f(x, nd=2):
    return "nan" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}f}"


def testo_report(c):
    gw, gm = c["gw"], c["gm"]
    righe = c["righe"]
    validi = [r for r in righe if not math.isnan(r["delta"])]
    su = [r for r in validi if r["flag"] == "UP"]
    giu = [r for r in validi if r["flag"] == "DOWN"]
    sito = next((r for r in righe if r["sito"]), None)
    loc = [r for r in validi if abs(r["resid"] - c["pos"]) <= INTORNO and not r["sito"]]
    top = sorted([r for r in validi if r["flag"]], key=lambda r: -abs(r["delta"]))[:8]

    L = []
    L.append(f"CONFRONTO Q PER RESIDUO: {c['sigla']} vs WT ({c['wt_tipo']})"
             + (f"   (tag mutations.txt: {c['tag']})" if c["tag"] else ""))
    L.append("=" * 72)
    L.append(f"WT di riferimento : {c['wt_label']}")
    L.append(f"Finestra di frame : {c['finestra']}")
    L.append(f"Repliche          : WT {len(gw)} ({', '.join(c['wt_run_names'])}) | "
             f"{c['sigla']} {len(gm)} ({', '.join(c['mut_run_names'])})")
    if c["avviso_res"]:
        L.append(f"ATTENZIONE        : {c['avviso_res']}")
    L.append("")
    L.append("Q globale per replica (media pesata dei q_res sui contatti nativi, %)")
    L.append("  WT          : " + "  ".join(f(x) for x in gw))
    L.append(f"  {c['sigla']:<11} : " + "  ".join(f(x) for x in gm))
    L.append(f"  Media +- SEM: WT {f(media(gw))} +- {f(sem(gw))}  |  "
             f"{c['sigla']} {f(media(gm))} +- {f(sem(gm))}")
    L.append(f"  Differenza  : {c['delta']:+.2f} punti di Q (mutante - WT)")
    L.append(f"  p-value     : {c['p']:.4f}  (test di permutazione {c['metodo']}, "
             f"{c['n_perm']} permutazioni, due code)")
    L.append("")
    L.append(f"ESITO: {c['esito']}  (soglia p < {c['alpha']})")
    if c["esito"] == "STABILIZZANTE":
        L.append(f"  -> {c['sigla']} conserva piu' contatti nativi del WT: il Q indica un "
                 f"AUMENTO di stabilita'.")
    elif c["esito"] == "DESTABILIZZANTE":
        L.append(f"  -> {c['sigla']} perde piu' contatti nativi del WT: il Q indica una "
                 f"DIMINUZIONE di stabilita'.")
    else:
        tendenza = "in aumento" if c["delta"] > 0 else "in diminuzione"
        L.append(f"  -> Differenza non distinguibile dalla variabilita' tra repliche "
                 f"(tendenza del Q: {tendenza}, ma non supportata).")
    L.append("")
    L.append("Dettaglio per residuo (descrittivo, nessuna correzione per test multipli)")
    L.append(f"  Residui con Q piu' alto del WT di oltre 2 errori standard : {len(su)}")
    L.append(f"  Residui con Q piu' basso del WT di oltre 2 errori standard: {len(giu)}")
    if sito and not math.isnan(sito["delta"]):
        L.append(f"  Sito mutato ({sito['resname']}{sito['resid']}): Q WT {f(sito['q_wt'])} -> "
                 f"{f(sito['q_mut'])}  (delta {sito['delta']:+.2f}; contatti nativi "
                 f"{f(sito['n_wt'], 1)} -> {f(sito['n_mut'], 1)})")
    if loc:
        L.append(f"  Intorno +-{INTORNO} residui (sito escluso): delta medio "
                 f"{media([r['delta'] for r in loc]):+.2f}")
    if top:
        L.append("  Residui con differenza maggiore (flag):")
        for r in top:
            L.append(f"    {r['resname']}{r['resid']:<4} delta {r['delta']:+7.2f}  "
                     f"(WT {f(r['q_wt'])} -> {f(r['q_mut'])})  {r['flag']}")
    L.append("")
    L.append("Note per l'interpretazione")
    L.append("  - Q e' relativo al frame 0 di OGNI traiettoria: misura quanto ciascun sistema "
             "conserva la propria")
    L.append("    struttura iniziale. L'insieme dei contatti nativi del mutante puo' differire da "
             "quello del WT.")
    if c["wt_tipo"] == "controllo":
        L.append("  - Il WT e' il controllo passato nella STESSA pipeline dei mutanti (struttura WT_ di "
                 "FoldX, stessa")
        L.append("    costruzione del sistema, cella misurata, conf, seed e analisi).")
    else:
        L.append("  - ATTENZIONE: WT originale. E' stato simulato con la cella scritta a mano in "
                 "sample.conf e senza")
        L.append("    passare da FoldX, mentre i mutanti usano la cella misurata e la pipeline FoldX: "
                 "una differenza")
        L.append("    sistematica di protocollo puo' spostare il Q di tutte le mutazioni nello stesso "
                 "verso.")
    L.append("  - Q a 510 K su 1 ns e' un indicatore della resistenza al disordine termico, "
             "non una misura")
    L.append("    di stabilita' termodinamica (ddG). Con poche repliche la potenza del test e' limitata.")
    return "\n".join(L) + "\n"


def scrivi_tsv(c, path):
    cols = ["resid", "resname", "sito_mutato", "n_nativi_wt", "n_nativi_mut", "q_wt", "sem_wt",
            "q_mut", "sem_mut", "delta", "flag"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in c["righe"]:
            fh.write("\t".join([str(r["resid"]), r["resname"], "yes" if r["sito"] else "",
                                f(r["n_wt"], 1), f(r["n_mut"], 1), f(r["q_wt"]), f(r["sem_wt"]),
                                f(r["q_mut"]), f(r["sem_mut"]), f(r["delta"]), r["flag"]]) + "\n")


SUMMARY_COLS = ["mutazione", "wt", "tag", "n_wt", "n_mut", "q_wt", "q_mut", "delta", "p_value",
                "esito", "res_up", "res_down"]


def aggiorna_summary(risultati):
    path = EVAL_DIR / "summary.tsv"
    righe = {}
    if path.exists():
        vecchie = path.read_text().splitlines()
        if vecchie and vecchie[0].split("\t") == SUMMARY_COLS:
            for r in vecchie[1:]:
                p = r.split("\t")
                if len(p) > 1:
                    righe[(p[0], p[1])] = r
        else:
            print("ATTENZIONE: summary.tsv in formato precedente, rigenerato con i soli "
                  "risultati di questo lancio", flush=True)
    for c in risultati:
        v = [r for r in c["righe"] if not math.isnan(r["delta"])]
        righe[(c["sigla"], c["wt_tipo"])] = "\t".join([
            c["sigla"], c["wt_tipo"], c["tag"] or "-", str(len(c["gw"])), str(len(c["gm"])),
            f(media(c["gw"])), f(media(c["gm"])), f"{c['delta']:+.2f}", f"{c['p']:.4f}",
            c["esito"].replace(" ", "_"),
            str(sum(1 for r in v if r["flag"] == "UP")),
            str(sum(1 for r in v if r["flag"] == "DOWN"))])
    with open(path, "w") as fh:
        fh.write("\t".join(SUMMARY_COLS) + "\n")
        for k in sorted(righe):
            fh.write(righe[k] + "\n")
    return path


# ------------------------------- main -------------------------------
def trova_mutazioni():
    if not MUT_ROOT.is_dir():
        return []
    return sorted(d.name for d in MUT_ROOT.iterdir()
                  if d.is_dir() and SIGLA_RE.fullmatch(d.name)
                  and any((d / "replicas").glob(f"run*/{FILE_Q}")))


def main():
    ap = argparse.ArgumentParser(description="Confronto Q per residuo mutanti vs WT")
    ap.add_argument("mutazioni", nargs="*", help="sigle (default: tutte quelle con dati)")
    ap.add_argument("--wt-dir", type=Path, default=None,
                    help="cartella del WT da usare (con run*/ oppure con replicas/run*/); "
                         "default: mutations/WT se esiste, altrimenti replicas/")
    ap.add_argument("--alpha", type=float, default=ALPHA, help=f"soglia p-value ({ALPHA})")
    args = ap.parse_args()

    try:
        if args.wt_dir:
            wt_dir = args.wt_dir.resolve()
            if (wt_dir / "replicas").is_dir():
                wt_dir = wt_dir / "replicas"
        elif raccogli(MUT_ROOT / "WT" / "replicas"):
            wt_dir = MUT_ROOT / "WT" / "replicas"
        else:
            wt_dir = WT_REPLICAS
            print("ATTENZIONE: controllo mutations/WT non trovato, uso il WT originale (replicas/). "
                  "Lancia 'python3 namd_mutations.py WT' per il confronto a pari protocollo.",
                  flush=True)
        wt_tipo = "originale" if wt_dir.resolve() == WT_REPLICAS.resolve() else "controllo"
        wt_label = f"{wt_tipo} ({wt_dir})"
        print(f"WT di riferimento: {wt_label}", flush=True)
        wt_runs = raccogli(wt_dir)
        if not wt_runs:
            raise ConfrontoError(
                f"nessun {FILE_Q} in {wt_dir}/run*/. Calcolalo prima con "
                f"q_value_per_residue.tcl su ogni replica del WT")
        sigle = [s.upper() for s in args.mutazioni] or trova_mutazioni()
        if not sigle:
            raise ConfrontoError(f"nessuna mutazione con {FILE_Q} in {MUT_ROOT}/<SIGLA>/replicas/run*/")
    except ConfrontoError as e:
        print(f"ERRORE: {e}", file=sys.stderr)
        return 1

    tags = carica_tag()
    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    risultati, problemi = [], []
    n_run_wt = len(list(wt_dir.glob("run*")))
    if len(wt_runs) < n_run_wt:
        print(f"ATTENZIONE: WT con {FILE_Q} in {len(wt_runs)} repliche su {n_run_wt}", flush=True)

    for s in sigle:
        try:
            if not SIGLA_RE.fullmatch(s):
                raise ConfrontoError("formato sigla non valido (es. AA95M)")
            rep = MUT_ROOT / s / "replicas"
            mut_runs = raccogli(rep) if rep.is_dir() else []
            if not mut_runs:
                raise ConfrontoError(f"nessun {FILE_Q} in {rep}/run*/")
            n_run_mut = len(list(rep.glob("run*")))
            if len(mut_runs) < n_run_mut:
                print(f"ATTENZIONE [{s}]: {FILE_Q} in {len(mut_runs)} repliche su {n_run_mut}", flush=True)
            c = confronta(s, wt_runs, mut_runs, tags.get(s), args.alpha, wt_tipo, wt_label)
        except ConfrontoError as e:
            print(f"[{s}] saltata: {e}", file=sys.stderr)
            problemi.append(s)
            continue
        rep_txt = testo_report(c)
        suffisso = "WTcontrollo" if wt_tipo == "controllo" else "WT"
        (EVAL_DIR / f"{s}_vs_{suffisso}.txt").write_text(rep_txt)
        scrivi_tsv(c, EVAL_DIR / f"{s}_vs_{suffisso}_per_residue.tsv")
        risultati.append(c)
        print(f"[{s}] Q globale WT {f(media(c['gw']))} -> {f(media(c['gm']))} "
              f"(delta {c['delta']:+.2f}, p = {c['p']:.4f})  ESITO: {c['esito']}", flush=True)

    if risultati:
        sp = aggiorna_summary(risultati)
        print(f"\nOutput in {EVAL_DIR}/  (report per mutazione, tabelle per residuo, {sp.name})")
    return 0 if not problemi else 2


if __name__ == "__main__":
    sys.exit(main())
