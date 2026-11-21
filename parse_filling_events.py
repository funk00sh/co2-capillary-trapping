"""Parse the PNflow .prt logs of every run into one table of filling events and trapped-cluster statistics.

Usage: python parse_filling_events.py <netB_nw_work_dir>
Writes filling_events_all_runs.csv (A: raw_runs_462.zip, B: old_network_DongBlunt/raw_runs_306.zip + new 95-115 block).
"""
import re
import sys
import zipfile
from pathlib import Path
from typing import Dict, Iterator, Tuple

import pandas as pd

HERE = Path(__file__).resolve().parent
KEYS = {
    "p_snap": r"Pore Filling Process[\s\S]*?Snap off:\s+(\d+)",
    "p_piston": r"Pore Filling Process[\s\S]*?Piston type displacement:\s+(\d+)",
    "p_I2": r"Pore body filling, I2:\s+(\d+)",
    "p_I3": r"Pore body filling, I3:\s+(\d+)",
    "p_I4": r"Pore body filling, I4\+:\s+(\d+)",
    "p_uninv": r"Pore Filling Process[\s\S]*?Uninvaded:\s+(\d+)",
    "t_snap": r"Throat Filling Process[\s\S]*?Snap off:\s+(\d+)",
    "t_piston": r"Throat Filling Process[\s\S]*?Piston type displacement:\s+(\d+)",
    "t_uninv": r"Throat Filling Process[\s\S]*?Uninvaded:\s+(\d+)",
    "trapped_oil_regions": r"Number of trapped oil regions:\s+(\d+)",
    "trapped_oil_elems": r"Total number of trapped oil elements:\s+(\d+)",
    "invaded_imb": r"Number of elements invaded:\s+(\d+)",
}


def parse(prt: str) -> Dict[str, int]:
    """Counts from the imbibition section of one log."""
    sec = prt.split("Imbibition  =====")[-1]
    return {k: (int(m.group(1)) if (m := re.search(p, sec)) else -1) for k, p in KEYS.items()}


def case_info(name: str) -> Tuple[str, float, int]:
    """Angle label, maximum Pc (kPa), seed from a case name."""
    m = re.match(r"(\w+_\d+-\d+)_ift(\d+)_pc([\d.]+)_s(\d+)", name)
    return m.group(1), float(m.group(3)), int(m.group(4))


def logs() -> Iterator[Tuple[str, str, str]]:
    """Yield (network, case, prt text) for every run."""
    for net, zp in (("A", HERE / "raw_runs_462.zip"), ("B", HERE / "old_network_DongBlunt/raw_runs_306.zip")):
        z = zipfile.ZipFile(zp)
        for n in z.namelist():
            if n.endswith("run_pnflow.prt") and "_ift37_" in n:
                yield net, n.split("/")[0], z.read(n).decode(errors="replace")
    for d in Path(sys.argv[1]).glob("NW_95-115_*"):
        yield "B", d.name, (d / "run_pnflow.prt").read_text(errors="replace")


def main() -> None:
    """Build the table."""
    rows = []
    for net, case, prt in logs():
        w, pc, seed = case_info(case)
        rows.append(dict(network=net, case=case, wett=w, pc_max_kPa=pc, seed=seed, **parse(prt)))
    df = pd.DataFrame(rows)
    sm = pd.concat([pd.read_csv(HERE / "summary_all_runs.csv").assign(network="A"),
                    pd.read_csv(HERE / "old_network_DongBlunt/summary_all_runs.csv").assign(network="B"),
                    pd.read_csv(Path(sys.argv[1]) / "summary.csv").assign(network="B")])
    df = df.merge(sm[["network", "case", "Sgi", "Sgr"]], on=["network", "case"], how="left")
    df.to_csv(HERE / "filling_events_all_runs.csv", index=False)
    print(len(df), df.groupby("network").size().to_dict(), df.Sgr.isna().sum(), (df == -1).sum().sum())


if __name__ == "__main__":
    main()
