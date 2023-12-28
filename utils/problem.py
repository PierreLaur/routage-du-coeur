from dataclasses import dataclass
import pandas as pd


@dataclass
class Problem:
    matrix: pd.DataFrame
    n: int
    n_pdr: int
    m: int
    n_days: int
    demands: int
    freqs_pdr: int
    capacities: list[int]
    sizes: list[int]
    frais: list[bool]
    max_palette_capacity: int
