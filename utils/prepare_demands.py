from collections import defaultdict
from math import ceil
import pandas as pd
from utils.datatypes import Demand, LoadType, Params, ProductType

frais_sur_demi_palettes = [
    "L ISLE EN DODON",
    "RIEUMES",
    "LUCHON",
    "CAZERES",
    "SALIES DU SALAT",
    "CARBONNE",
    "LEGUEVIN",
    "LEVIGNAC",
    "REVEL",
    "PLAISANCE",
]


def prepare_demands(demands_file: str, params: Params):
    demands_df = pd.read_excel(demands_file, index_col=0)
    demands: dict[int, list[Demand]] = defaultdict(list)
    for c, (index, row) in enumerate(demands_df.iterrows()):
        weight_A = ceil(row["A"])
        while weight_A > 0:
            new_weight = min(weight_A, params.max_palette_capacity)
            weight_A -= new_weight
            new_demand = Demand(new_weight, LoadType.PALETTE, ProductType.A)
            demands[c].append(new_demand)

        weight_F = ceil(row["F"])
        while weight_F > 0:
            if weight_F > params.demi_palette_capacity:
                new_weight = min(weight_F, ceil(params.max_palette_capacity / 2))
                load_type = LoadType.PALETTE
                weight_F -= new_weight
            else:
                new_weight = min(weight_F, ceil(params.demi_palette_capacity / 2))
                load_type = LoadType.DEMI_PALETTE
                weight_F -= new_weight

            new_demand = Demand(new_weight, load_type, ProductType.F)
            demands[c].append(new_demand)

        weight_S = ceil(row["S"])
        if weight_S <= 120 or index == "BESSIERES":
            # use norvegiennes
            while weight_S > 0:
                new_weight = min(weight_S, params.norvegienne_capacity)
                weight_S -= new_weight
                new_demand = Demand(new_weight, LoadType.NORVEGIENNE, ProductType.S)
                demands[c].append(new_demand)
        else:
            # use demi-palettes
            while weight_S > 0:
                new_weight = min(weight_S, params.demi_palette_capacity)
                weight_S -= new_weight
                new_demand = Demand(new_weight, LoadType.DEMI_PALETTE, ProductType.S)
                demands[c].append(new_demand)

    return demands
