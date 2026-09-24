"""
lookup.py - Replicon-specific risk priors derived from PIPdb.

The lookup table provides empirical median values for ALL ten risk
dimensions (S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE, S_BMG, S_GEO,
S_HAB, S_GROW) for known replicon types, based on analysis of
792,964 PSCs in PIPdb.

When sequence-derived dimensions cannot be computed (e.g., abricate not
available, no mobility genes detected), these replicon-specific medians
are used as fallback priors.
"""

from __future__ import annotations

import os
from typing import Optional

import pandas as pd

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
_LOOKUP_FILE = os.path.join(_DATA_DIR, "replicon_lookup.csv")

# All ten dimensions that may have replicon priors
ALL_DIMENSIONS = (
    "S_ARG", "S_VF", "S_MOB", "S_HOST", "S_REP",
    "S_SIZE", "S_BMG", "S_GEO", "S_HAB", "S_GROW",
)


def load_replicon_lookup(path: Optional[str] = None) -> pd.DataFrame:
    """
    Load the replicon prior lookup table.

    Parameters
    ----------
    path : str, optional
        Path to a custom CSV file. If None, uses the bundled table.

    Returns
    -------
    pd.DataFrame with columns:
        replicon_primary, n_PSC, and all ten S_* dimension columns
        containing empirical median values from PIPdb.
    """
    filepath = path or _LOOKUP_FILE
    df = pd.read_csv(filepath)
    required = {"replicon_primary", "n_PSC"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Lookup table missing columns: {missing}")
    # Ensure all dimension columns exist; fill missing with 0
    for dim in ALL_DIMENSIONS:
        if dim not in df.columns:
            df[dim] = 0.0
    return df
