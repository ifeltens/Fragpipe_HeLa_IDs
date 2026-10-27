#!/usr/bin/env python3
"""
Convert new files in a folder to mzML with MSConvert (Peak Picking, Vendor),
run FragPipe on them, log the number of identifications, and plot them over
time against all previously processed files.

Usage:
    1. Edit the CONFIG section below (once).
    2. Drop your new .raw / .d / .wiff (or .mzML) file into WATCH_DIR.
    3. Run:  python run_fragpipe_qc.py

Requires: Python 3.8+, matplotlib  (pip install matplotlib)
"""

import csv
import datetime as dt
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt

# ============================== CONFIG ======================================
FRAGPIPE_BIN = r"C:\FragPipe\FragPipe-24.0\bin\fragpipe.bat"  # path to fragpipe launcher
WORKFLOW     = r"C:\Users\Izzy\Temp_Data\FragPipe_Methods\HeLaDigest_BasicSearch.workflow"  # workflow saved from the FragPipe GUI
WATCH_DIR    = Path(r"C:\Users\Izzy\Temp_Data\HeLa_digest")  # folder you drop files into
RESULTS_DIR  = Path(r"C:\Users\Izzy\Temp_Data\HeLaDigest_Output")  # one subfolder per file

# MSConvert (ProteoWizard) settings
MSCONVERT_BIN = r"C:\Users\Izzy\ProteoWizard\msconvert.exe"
MZML_DIR      = Path(r"C:\Users\Izzy\Temp_Data\HeLa_digest_mzML")  # converted files go here (NOT inside WATCH_DIR)
PEAK_PICKING_FILTER = "peakPicking vendor msLevel=1-"  # Peak Picking, Vendor, MS levels 1-
HISTORY_CSV  = RESULTS_DIR / "id_history.csv"          # running log of all runs
PLOT_PNG     = RESULTS_DIR / "ids_over_time.png"

DATA_TYPE    = "DDA"            # DDA, DIA, DIA-Quant, DDA+ ...
EXTENSIONS   = {".raw", ".d", ".wiff", ".mzml", ".mzxml"}   # files picked up from WATCH_DIR
NEEDS_CONVERSION = {".raw", ".d", ".wiff"}   # these get converted to mzML first; mzML/mzXML are used as-is

# Process ONLY this file (name without the extension). It is re-run every time you
# run the script, and its row in id_history.csv is replaced.
# Set to None to go back to processing every new file in WATCH_DIR.
TARGET_STEM  = "HeLa_Digest_test"
THREADS      = 4
RAM_GB       = 10               # fixed limit in GB (0 = let FragPipe decide)
EXCLUDE_CONTAMINANTS = True     # drop contam_/Cont_ proteins from protein count
# Use "acquired" date from file modified time. Change to "now" to use run date.
DATE_SOURCE  = "mtime"          # "mtime" or "now"
# ============================================================================

FIELDS = ["file", "date", "psms", "peptides", "proteins"]


def read_history():
    if not HISTORY_CSV.exists():
        return []
    with open(HISTORY_CSV, newline="") as f:
        return list(csv.DictReader(f))


def write_history(rows):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(HISTORY_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)


def find_new_files(done):
    candidates = [p for p in WATCH_DIR.iterdir() if p.suffix.lower() in EXTENSIONS]
    if TARGET_STEM:
        matches = [p for p in candidates if p.stem.lower() == TARGET_STEM.lower()]
        if not matches:
            print(f"No file named {TARGET_STEM} (with a supported extension) in {WATCH_DIR}")
        return sorted(matches)
    return sorted(p for p in candidates if p.name not in done)


def convert_to_mzml(src: Path) -> Path:
    """Convert a vendor file to mzML with MSConvert (Peak Picking, Vendor).
    Returns the path to the mzML file. mzML/mzXML inputs are returned unchanged."""
    if src.suffix.lower() not in NEEDS_CONVERSION:
        return src

    MZML_DIR.mkdir(parents=True, exist_ok=True)
    out = MZML_DIR / f"{src.stem}.mzML"
    if out.exists():
        print(f"\n>>> {out.name} already converted, skipping MSConvert")
        return out

    cmd = [
        MSCONVERT_BIN, str(src),
        "--mzML",
        "--filter", PEAK_PICKING_FILTER,
        "-o", str(MZML_DIR),
    ]
    print(f"\n>>> Converting {src.name} to mzML with MSConvert")
    try:
        subprocess.run(cmd, check=True)
        if not out.exists():
            raise FileNotFoundError(f"MSConvert finished but {out} was not created")
    except Exception:
        out.unlink(missing_ok=True)  # don't leave a partial file that would be skipped next time
        raise
    return out


def run_fragpipe(data_file: Path) -> Path:
    stem = data_file.stem
    workdir = RESULTS_DIR / stem
    workdir.mkdir(parents=True, exist_ok=True)

    manifest = workdir / f"{stem}.fp-manifest"
    # manifest columns: path <tab> experiment <tab> bioreplicate <tab> data type
    manifest.write_text(f"{data_file.resolve()}\t{stem}\t\t{DATA_TYPE}\n")

    cmd = [
        FRAGPIPE_BIN, "--headless",
        "--workflow", WORKFLOW,
        "--manifest", str(manifest),
        "--workdir", str(workdir),
        "--threads", str(THREADS),
        "--ram", str(RAM_GB),
    ]
    print(f"\n>>> Running FragPipe on {data_file.name}")
    subprocess.run(cmd, check=True)
    return workdir


def count_rows(tsv: Path, skip_contaminants=False):
    n = 0
    with open(tsv, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if skip_contaminants:
                pid = (row.get("Protein") or row.get("Protein ID") or "").lower()
                if pid.startswith(("contam_", "cont_")):
                    continue
            n += 1
    return n


def count_ids(workdir: Path):
    def find(name):
        hits = sorted(workdir.rglob(name))
        if not hits:
            raise FileNotFoundError(f"{name} not found under {workdir}")
        return hits[0]

    return (
        count_rows(find("psm.tsv")),
        count_rows(find("peptide.tsv")),
        count_rows(find("protein.tsv"), EXCLUDE_CONTAMINANTS),
    )


def file_date(p: Path):
    if DATE_SOURCE == "now":
        return dt.datetime.now()
    return dt.datetime.fromtimestamp(p.stat().st_mtime)


def make_plot(rows):
    if not rows:
        print("Nothing to plot yet.")
        return
    rows = sorted(rows, key=lambda r: r["date"])
    dates = [dt.datetime.fromisoformat(r["date"]) for r in rows]
    series = [("psms", "PSMs"), ("peptides", "Peptides"), ("proteins", "Proteins")]

    fig, axes = plt.subplots(len(series), 1, figsize=(10, 8), sharex=True)
    for ax, (key, label) in zip(axes, series):
        vals = [int(r[key]) for r in rows]
        ax.plot(dates, vals, "o-", color="tab:blue", alpha=0.6)
        ax.plot(dates[-1], vals[-1], "o", color="tab:red", markersize=10,
                label=f"Latest: {vals[-1]:,}")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.legend(loc="best")
    axes[-1].set_xlabel("Date")
    fig.autofmt_xdate()
    fig.suptitle("FragPipe identifications over time")
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    print(f"Plot saved to {PLOT_PNG}")
    try:
        plt.show()
    except Exception:
        pass


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    history = read_history()
    done = {r["file"] for r in history}

    new_files = find_new_files(done)
    if not new_files:
        print("No new files found in", WATCH_DIR)
    for f in new_files:
        try:
            mzml = convert_to_mzml(f)
            workdir = run_fragpipe(mzml)
            psms, peps, prots = count_ids(workdir)
        except Exception as e:
            print(f"!!! Failed on {f.name}: {e}", file=sys.stderr)
            continue
        history = [r for r in history if r["file"] != f.name]  # replace any earlier row for this file
        history.append({
            "file": f.name,
            "date": file_date(f).isoformat(timespec="seconds"),
            "psms": psms, "peptides": peps, "proteins": prots,
        })
        write_history(history)  # save after each file in case a later one fails
        print(f"    {f.name}: {psms:,} PSMs, {peps:,} peptides, {prots:,} proteins")

    make_plot(history)


if __name__ == "__main__":
    main()
