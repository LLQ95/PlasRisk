"""
lookup.py - Replicon-specific risk priors derived from PIPdb.

The lookup table provides S_REP, S_GEO, S_HAB, S_GROW, and S_HOST
prior values for known replicon types, based on empirical analysis
of 792,964 PSCs in PIPdb.
"""

from __future__ import annotations

import os
from typing import Optional

import pandas as pd

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
_LOOKUP_FILE = os.path.join(_DATA_DIR, "replicon_lookup.csv")


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
        replicon_primary, S_REP, S_GEO, S_HAB, S_GROW, S_HOST,
        n_PSC, notes
    """
    filepath = path or _LOOKUP_FILE
    df = pd.read_csv(filepath)
    required = {"replicon_primary", "S_REP", "S_GEO", "S_HAB", "S_GROW", "S_HOST"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Lookup table missing columns: {missing}")
    return df
