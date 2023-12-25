import pandas as pd
from plots import plot_tours
from cp_solver import solve_vrp
from routing_solver import route_vrp

import json


def str_to_tuple(my_dict):
    new_dict = {}
    for k, v in my_dict.items():
        key = k.strip("()").split(", ")
        key = tuple(map(int, key))
        new_dict[key] = v
    return new_dict


if __name__ == "__main__":
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

    # Add 15% for robustness
    for c in range(1, n):
        demands["a"][c] = int(demands["a"][c] * 1.15)
        demands["f"][c] = int(demands["f"][c] * 1.15)
        demands["s"][c] = int(demands["s"][c] * 1.15)

    for d in demands.values():
        d[0] = 0

    freqs_pdr = (
        points_de_ramasse["Fréquence de Ramasse(/w)"].fillna(0).astype(int).tolist()
    )

    frais = [f == "Oui" for f in vehicles["Frais"]]

    current_tours = open("data/tours_tournees_actuelles.json", "r")
    current_tours = json.load(current_tours)
    current_tours = str_to_tuple(current_tours)

    current_arcs = json.load(open("data/arcs_tournees_actuelles.json", "r"))
    current_arcs = str_to_tuple(current_arcs)

    tours, obj, deliveries, visits, arcs = solve_vrp(
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
        hint=(current_tours, current_arcs),
    )

    # obj, tours = route_vrp(
    #     matrix,
    #     n,
    #     n_pdr,
    #     m,
    #     n_days,
    #     demands,
    #     freqs_pdr,
    #     capacities,
    #     sizes,
    #     frais,
    #     max_palette_capacity,
    #     hint=(current_tours, current_arcs),
    # )
    # for t in tours.values():
    #     if len(t) > 2:
    #         print(t)
    # exit()

    # check_solution(centres, vehicles, tours, obj, deliveries, visits, arcs)

    for d in range(n_days):
        print(f"-- JOUR {d} --")
        for v in range(m):
            if len(tours[v, d]) == 1:
                continue
            print(f"Vehicule {v} ({vehicles.index[v]}) tour : \n", end="")
            for t in tours[v, d]:
                delivery = None
                if t >= n:
                    name = points_de_ramasse.index[t - n]
                elif t == 0:
                    continue
                else:
                    name = centres["Nom"][t]
                    delivery = tuple(deliveries[i][v, d, t] for i in range(3))
                print(f"\t{name:20} {str(delivery) if delivery else ''}")

    plot_tours([tour for tour in tours.values() if len(tour) > 1])
