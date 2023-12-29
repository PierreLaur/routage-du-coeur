from utils.plots import plot_tours
from utils.problem import Problem, read_problem
from models.cp_solver import solve_vrp, solve_vrp_single_serve
from models.routing_solver import route_vrp
import json
import argparse


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


def str_to_tuple(my_dict):
    new_dict = {}
    for k, v in my_dict.items():
        key = k.strip("()").split(", ")
        key = tuple(map(int, key))
        new_dict[key] = v
    return new_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("week", type=int, help="week number (1 or 2)")
    parser.add_argument("--infile", type=str, help="init solution")
    parser.add_argument("--outfile", type=str, help="desired output file path")

    args = parser.parse_args()

    if not args.outfile:
        print("Warning : No output file specified")

    problem = read_problem(
        "data/centres.xlsx",
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/euclidean_matrix.xlsx",
        week=2,
    )

    tours, obj = solve_vrp(problem, hint=args.infile, outfile=args.outfile)
