from dataclasses import dataclass
import pandas as pd
from math import ceil
import json
import random


@dataclass
class Problem:
    distance_matrix: pd.DataFrame
    duration_matrix: pd.DataFrame
    n: int
    n_pdr: int
    m: int
    n_days: int
    demands: list[int]
    weights: list[int]
    j_de_livraison_possibles: list[list[int]]
    j_de_ramasse: list[list[int]]
    capacities: list[int]
    consumptions: list[int]
    sizes: list[int]
    frais: list[bool]
    max_palette_capacity: int
    demi_palette_capacity: int
    n_norvegiennes: int
    norvegienne_capacity: int
    max_stops: int
    max_tour_duration: int
    max_first_pickup_time: int
    wait_at_centres: int
    wait_at_pdrs: int
    duration_coefficients: list[int]
    vehicle_allowed: list[bool]
    disallow_norvegiennes_in_PL: bool
    weekly_fixed_cost: float
    fuel_cost: float


def read_problem(
    centres_file,
    points_de_ramasse_file,
    vehicles_file,
    distance_matrix_file,
    duration_matrix_file,
    params_file,
    week,
    custom_assignment=None,
):
    """Reads the problem from the usual input files and returns a Problem object"""
    centres = pd.read_excel(centres_file, index_col=0)
    points_de_ramasse = pd.read_excel(points_de_ramasse_file, index_col=0)
    vehicles = pd.read_excel(vehicles_file, index_col=0)
    distance_matrix = pd.read_excel(distance_matrix_file, index_col=0)
    duration_matrix = pd.read_excel(duration_matrix_file, index_col=0)
    params = json.load(open(params_file, "r"))

    n = len(centres.index)
    n_pdr = len(points_de_ramasse.index)
    n_days = 5  # Single week scheduling
    m = len(vehicles.index)

    capacities = vehicles["Capacité (kg)"].astype(int).tolist()
    consumptions = vehicles["Consommation (L/100km)"].astype(int).tolist()
    sizes = vehicles["Taille(Palettes)"].astype(int).tolist()

    max_palette_capacity = params["max_palette_capacity"]
    demi_palette_capacity = params["demi_palette_capacity"]
    max_stops = params["max_stops"]
    norvegienne_capacity = params["norvegienne_capacity"]
    n_norvegiennes = params["n_norvegiennes"]

    demands = {
        "a": centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist(),
        "f": centres["Tonnage Frais (kg)"].fillna(0).astype(int).tolist(),
        "s": centres["Tonnage Surgelé (kg)"].fillna(0).astype(int).tolist(),
    }

    weights = points_de_ramasse["Poids par ramasse(kg)"].fillna(0).astype(int).tolist()

    if custom_assignment:
        params["week_assignment"] = custom_assignment

    # Ignore centres that are not delivered this week :
    i = 0
    for c in range(n):
        # 0 = every week
        if centres["Semaine"][c] == 0:
            continue

        if params["week_assignments"][i] != week:
            demands["a"][c] = 0
            demands["f"][c] = 0
            demands["s"][c] = 0
        i += 1

    # Make the demands higher for robustness
    for c in range(1, n):
        demands["a"][c] = ceil(demands["a"][c] * params["robustness_factor"])
        demands["f"][c] = ceil(demands["f"][c] * params["robustness_factor"])
        demands["s"][c] = ceil(demands["s"][c] * params["robustness_factor"])
    for p in range(n_pdr):
        weights[p] = ceil(weights[p])

    # Demand for depot is set to 0
    for d in demands.values():
        d[0] = 0

    # Allowed delivery days
    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4}
    jours_de_livraison = centres["Jours de Livraison possibles"].tolist()
    jours_de_livraison = [[]] + [
        j.replace(" ", "").split(",") for j in jours_de_livraison[1:]
    ]
    j_de_livraison_possibles = [
        list(map(lambda x: jours_map[x], j)) for j in jours_de_livraison
    ]

    # Mandatory pickup days
    jours_de_ramasse = points_de_ramasse[f"Jours de Ramasse"].tolist()
    j_de_ramasse = []
    for pdr in range(n_pdr):
        jours = jours_de_ramasse[pdr].split(", ")
        for i, j in enumerate(jours):
            jours[i] = jours_map[jours[i]]
        j_de_ramasse.append(jours)

    frais = [f == "Oui" for f in vehicles["Frais"]]

    problem = Problem(
        distance_matrix,
        duration_matrix,
        n,
        n_pdr,
        m,
        n_days,
        demands,
        weights,
        j_de_livraison_possibles,
        j_de_ramasse,
        capacities,
        consumptions,
        sizes,
        frais,
        max_palette_capacity,
        demi_palette_capacity,
        n_norvegiennes,
        norvegienne_capacity,
        max_stops,
        params["max_tour_duration"],
        params["max_first_pickup_time"],
        params["wait_at_centres"],
        params["wait_at_pdrs"],
        params["duration_coefficients"],
        [bool(x) for x in params["vehicle_allowed"]],
        params["disallow_norvegiennes_in_PL"],
        params["weekly_fixed_cost"],
        params["fuel_cost"],
    )

    return problem
