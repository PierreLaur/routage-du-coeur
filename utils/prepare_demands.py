from collections import defaultdict
from math import ceil
import pandas as pd
from utils.datatypes import Demand, Params, ProductType

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
            new_demand = Demand(new_weight, 1, 0, ProductType.A)
            demands[c].append(new_demand)

        weight_F = ceil(row["F"])
        while weight_F > 0:
            if weight_F > params.demi_palette_capacity:
                new_weight = min(weight_F, ceil(params.max_palette_capacity / 2))
                palettes = 1
                weight_F -= new_weight
            else:
                new_weight = min(weight_F, ceil(params.demi_palette_capacity / 2))
                palettes = 0.5
                weight_F -= new_weight

            new_demand = Demand(new_weight, palettes, 0, ProductType.F)
            demands[c].append(new_demand)

        weight_S = ceil(row["S"])
        if weight_S <= 120 or index == "BESSIERES":
            # use norvegiennes
            while weight_S > 0:
                new_weight = min(weight_S, params.norvegienne_capacity)
                weight_S -= new_weight
                new_demand = Demand(new_weight, 0, 1, ProductType.S)
                demands[c].append(new_demand)
        else:
            # use demi-palettes
            while weight_S > 0:
                new_weight = min(weight_S, params.demi_palette_capacity)
                weight_S -= new_weight
                new_demand = Demand(new_weight, 0.5, 0, ProductType.S)
                demands[c].append(new_demand)

    return demands


def prepare_demands_no_split(demands_file: str, params: Params):
    demands_df = pd.read_excel(demands_file, index_col=0)
    demands: dict[int, list[Demand]] = defaultdict(list)
    for c, (index, row) in enumerate(demands_df.iterrows()):
        weight = ceil(row["A"])
        palettes = ceil(weight / params.max_palette_capacity)
        new_demand = Demand(weight, palettes, 0, ProductType.A)

        if weight > 0:
            demands[c].append(new_demand)

        weight = ceil(row["F"])
        if index in frais_sur_demi_palettes:
            palettes = ceil(2 * weight / params.demi_palette_capacity) / 2
        else:
            palettes = ceil(2 * weight / params.max_palette_capacity)
        new_demand = Demand(weight, palettes, 0, ProductType.F)
        if weight > 0:
            if weight > 800:
                new_demand2 = Demand(700, 2, 0, ProductType.F)
                new_demand.weight -= 700
                new_demand.palettes -= 1
                demands[c].append(new_demand)
                demands[c].append(new_demand2)
            else:
                demands[c].append(new_demand)

        weight = ceil(row["S"])
        if weight <= 120 or index == "BESSIERES":
            # use norvegiennes
            norvegiennes = ceil(weight / params.norvegienne_capacity)
            new_demand = Demand(weight, 0, norvegiennes, ProductType.S)
        else:
            # use demi-palettes
            palettes = ceil(2 * weight / params.demi_palette_capacity) / 2
            new_demand = Demand(weight, palettes, 0, ProductType.S)
        if weight > 0:
            demands[c].append(new_demand)

    return demands
