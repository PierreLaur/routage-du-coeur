from dataclasses import dataclass
import pandas as pd
from math import ceil


@dataclass
class Problem:
    matrix: pd.DataFrame
    n: int
    n_pdr: int
    m: int
    n_days: int
    demands: list[int]
    weights: list[int]
    j_de_ramasse: list[list[int]]
    use_pl: set[tuple[int, int]]
    capacities: list[int]
    sizes: list[int]
    frais: list[bool]
    max_palette_capacity: int


def read_problem(
    centres_file, points_de_ramasse_file, vehicles_file, matrix_file, week
):
    centres = pd.read_excel(centres_file, index_col=0)
    points_de_ramasse = pd.read_excel(points_de_ramasse_file, index_col=0)
    vehicles = pd.read_excel(vehicles_file, index_col=0)
    matrix = pd.read_excel(matrix_file, index_col=0)

    n = len(centres.index)
    n_pdr = len(points_de_ramasse.index)
    n_days = 5  # Single week scheduling
    m = len(vehicles.index)

    capacities = vehicles["Capacité (kg)"].astype(int).tolist()
    sizes = vehicles["Taille(Palettes)"].astype(int).tolist()
    max_palette_capacity = 800

    demands = {
        "a": centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist(),
        "f": centres["Tonnage Frais (kg)"].fillna(0).astype(int).tolist(),
        "s": centres["Tonnage Surgelé (kg)"].fillna(0).astype(int).tolist(),
    }

    weights = points_de_ramasse["Poids par ramasse(kg)"].fillna(0).astype(int).tolist()

    # Ignore centres that are not delivered this week :
    for i in range(n):
        # 0 = every week
        if centres["Semaine"][i] not in [week, 0]:
            demands["a"][i] = 0
            demands["f"][i] = 0
            demands["s"][i] = 0

    # Add 15% for robustness
    for c in range(1, n):
        demands["a"][c] = ceil(demands["a"][c] * 1.15)
        demands["f"][c] = ceil(demands["f"][c] * 1.15)
        demands["s"][c] = ceil(demands["s"][c] * 1.15)
    for p in range(n_pdr):
        weights[p] = ceil(weights[p] * 1.15)

    # Demand for depot is set to 0
    for d in demands.values():
        d[0] = 0

    jours_de_ramasse = points_de_ramasse[f"Jours de Ramasse S{week}"].tolist()
    j_de_ramasse = []
    use_pl = set()
    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4}
    for pdr in range(n_pdr):
        jours = jours_de_ramasse[pdr].split(", ")
        for i, j in enumerate(jours):
            if "(PL)" in j:
                jours[i] = jours[i].replace("(PL)", "")
                use_pl.add((jours_map[jours[i]], pdr))
            jours[i] = jours_map[jours[i]]
        j_de_ramasse.append(jours)

    frais = [f == "Oui" for f in vehicles["Frais"]]

    problem = Problem(
        matrix,
        n,
        n_pdr,
        m,
        n_days,
        demands,
        weights,
        j_de_ramasse,
        use_pl,
        capacities,
        sizes,
        frais,
        max_palette_capacity,
    )

    return problem
