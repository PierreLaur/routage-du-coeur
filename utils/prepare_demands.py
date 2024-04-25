from collections import defaultdict
import json
from math import ceil
import sys
import numpy as np
import pandas as pd
from utils.datatypes import Demand, Params, ProductType
from numpy.random import Generator, PCG64
from enum import Enum

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


class PalettePolicy(Enum):
    MAX = 0
    RANDOM = 1


def generate_weights(n_instances: int):
    distributions = pd.read_csv("data/demands/distributions.csv", index_col=0)
    rng = Generator(PCG64())
    n_centres = distributions.index.unique().shape[0]

    instances = np.zeros((n_instances, n_centres, 3))
    for centre in range(n_centres):
        for product in range(3):
            a, b, c = distributions.iloc[centre * 3 + product, 1:]
            if a == c:
                instances[:, centre, product] = a
            else:
                instances[:, centre, product] = rng.triangular(
                    a,
                    b,
                    c,
                    n_instances,
                )
    instances = np.ceil(instances).astype(int)
    return [
        pd.DataFrame(
            instances[i], index=distributions.index.unique(), columns=["A", "F", "S"]
        )
        for i in range(n_instances)
    ]


def palette_weight(
    quantity: int, capacity: int, policy: PalettePolicy = PalettePolicy.RANDOM
) -> int:
    if policy == PalettePolicy.RANDOM:
        capacity = ceil(np.random.uniform(capacity * 0.8, capacity))

    weight = min(quantity, capacity)
    return weight


def prepare_demands(demands_df: pd.DataFrame, params: Params, max=False):
    policy = PalettePolicy.MAX if max else PalettePolicy.RANDOM
    demands: dict[int, list[Demand]] = defaultdict(list)
    for c, (index, row) in enumerate(demands_df.iterrows()):
        weight_A = ceil(row["A"])
        while weight_A > 0:
            new_weight = palette_weight(weight_A, params.max_palette_capacity, policy)
            weight_A -= new_weight
            new_demand = Demand(new_weight, 1, 0, ProductType.A)
            demands[c].append(new_demand)

        weight_F = ceil(row["F"])
        while weight_F > 0:
            if weight_F > 0.5 * params.max_palette_capacity:
                new_weight = palette_weight(
                    weight_F, params.max_palette_capacity // 2, policy
                )
                palettes = 1
                weight_F -= new_weight
            else:
                new_weight = palette_weight(
                    weight_F, params.max_palette_capacity // 4, policy
                )
                palettes = 0.5
                weight_F -= new_weight

            new_demand = Demand(new_weight, palettes, 0, ProductType.F)
            demands[c].append(new_demand)

        weight_S = ceil(row["S"])
        if weight_S <= 120 or index == "BESSIERES":
            # use norvegiennes
            while weight_S > 0:
                new_weight = palette_weight(
                    weight_S, params.norvegienne_capacity, policy
                )
                weight_S -= new_weight
                new_demand = Demand(new_weight, 0, 1, ProductType.S)
                demands[c].append(new_demand)
        else:
            # use demi-palettes
            while weight_S > 0:
                new_weight = palette_weight(
                    weight_S, params.max_palette_capacity // 2, policy
                )
                weight_S -= new_weight
                new_demand = Demand(new_weight, 0.5, 0, ProductType.S)
                demands[c].append(new_demand)

    return demands


def prepare_demands_max_no_split(demands_df: pd.DataFrame, params: Params):
    demands: dict[int, list[Demand]] = defaultdict(list)
    for c, (index, row) in enumerate(demands_df.iterrows()):
        weight = ceil(row["A"])
        palettes = ceil(weight / params.max_palette_capacity)
        new_demand = Demand(weight, palettes, 0, ProductType.A)

        if weight > 0:
            demands[c].append(new_demand)

        weight = ceil(row["F"])
        if index in frais_sur_demi_palettes:
            palettes = ceil(4 * weight / params.max_palette_capacity) / 2
        else:
            palettes = ceil(2 * weight / params.max_palette_capacity)
        new_demand = Demand(weight, palettes, 0, ProductType.F)
        if weight > 0:
            if weight > 800:
                new_demand2 = Demand(700, 1, 0, ProductType.F)
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
            palettes = ceil(4 * weight / params.max_palette_capacity) / 2
            new_demand = Demand(weight, palettes, 0, ProductType.S)
        if weight > 0:
            demands[c].append(new_demand)

    return demands


if __name__ == "__main__":
    n_instances = int(sys.argv[1])

    params = Params.from_dict(json.load(open("data/params.json")))

    scenarios = {}

    weights = generate_weights(n_instances)
    for i in range(n_instances):
        demands = prepare_demands(weights[i], params)
        scenarios[i] = {k: [d.to_dict() for d in demands[k]] for k in demands}

    json.dump(scenarios, open("problems/demands/test_demands.json", "w"))
