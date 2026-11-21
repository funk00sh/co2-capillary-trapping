"""Extract imbibition event counts from a PNflow log; -1 indicates a missing field."""
import re
from typing import Dict

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

