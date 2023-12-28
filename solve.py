import pandas as pd
from utils.plots import plot_tours
from utils.problem import Problem
from models.cp_solver import solve_vrp, solve_vrp_single_serve
from models.routing_solver import route_vrp
import json
from math import ceil


def solve_with_routing(
    problem: Problem,
    current_tours,
    current_arcs,
):
    obj, tours = route_vrp(
        problem,
        hint=(current_tours, current_arcs),
    )

    for d in range(problem.n_days):
        print(f"-- JOUR {d} --")
        for v in range(problem.m):
            if (d, v) in tours and len(tours[d, v]) > 2:
                print(f"Vehicule {v} tour : \n", end="")
                for t in tours[d, v][1:-1]:
                    print("\t", problem.matrix.index[t])
    exit()


def tour_length(matrix, tour):
    length = 0
    for i in range(len(tour) - 1):
        length += matrix.iloc[tour[i], tour[i + 1]]
    length += matrix.iloc[tour[-1], tour[0]]
    return length


def print_tour_to_console(file):
    # TODO
    # for d in range(n_days):
    #     print(f"-- JOUR {d} --")
    #     for v in range(m):
    #         if len(tours[d, v]) == 1:
    #             continue
    #         print(f"Vehicule {v} ({vehicles.index[v]}) tour : \n", end="")

    #         for t in tours[d, v]:
    #             delivery = None
    #             pals = 0
    #             name = matrix.index[t]
    #             if t == 0 or t >= n:
    #                 continue
    #             else:
    #                 delivery = tuple(deliveries[i][d, v, t] for i in range(3))
    #                 pals = palettes[d, v, t]
    #             print(
    #                 f"\t{name:20} {str(delivery) if delivery else ''} {' - ' + str(pals)+' palettes'}"
    #             )
    pass


def solve_with_cp(ignore_centres):
    centres = pd.read_excel("data/centres.xlsx", index_col=0)
    points_de_ramasse = pd.read_excel("data/points_de_ramasse.xlsx", index_col=0)
    vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
    matrix = pd.read_excel("data/euclidean_matrix.xlsx", index_col=0)

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

    # Ignore some centres (deliver them the other week) :
    indexes = [centres.loc[centres["Nom"] == c].index.values[0] for c in ignore_centres]
    for i in indexes:
        demands["a"][i] = 0
        demands["f"][i] = 0
        demands["s"][i] = 0

    # Add 15% for robustness
    for c in range(1, n):
        demands["a"][c] = ceil(demands["a"][c] * 1.15)
        demands["f"][c] = ceil(demands["f"][c] * 1.15)
        demands["s"][c] = ceil(demands["s"][c] * 1.15)

    for d in demands.values():
        d[0] = 0

    freqs_pdr = (
        points_de_ramasse["Fréquence de Ramasse(/w)"].fillna(0).astype(int).tolist()
    )

    frais = [f == "Oui" for f in vehicles["Frais"]]

    problem = Problem(
        matrix,
        n,
        n_pdr,
        m,
        n_days,
        demands,
        freqs_pdr,
        capacities,
        sizes,
        frais,
        max_palette_capacity,
    )

    # current_tours = open("data/tours_tournees_actuelles_w1.json", "r")
    # current_tours = json.load(current_tours)
    # current_tours = str_to_tuple(current_tours)

    # current_arcs = json.load(open("data/arcs_tournees_actuelles.json", "r"))
    # current_arcs = str_to_tuple(current_arcs)

    tours, obj = solve_vrp(
        problem,
    )

    return obj, tours


def str_to_tuple(my_dict):
    new_dict = {}
    for k, v in my_dict.items():
        key = k.strip("()").split(", ")
        key = tuple(map(int, key))
        new_dict[key] = v
    return new_dict


if __name__ == "__main__":
    free = ["Bagneres de Luchon", "L Isle en Dodon"]
    centres_not_delivered_w1 = ["Cugnaux", "Levignac", "Rieumes"]
    centres_not_delivered_w2 = ["Fonsorbes", "Bessieres", "Fronton"]
    centres_semi_hebdo = free + centres_not_delivered_w1 + centres_not_delivered_w2

    obj, tours = solve_with_cp(ignore_centres=centres_not_delivered_w1)
