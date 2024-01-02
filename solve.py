from utils.problem import Problem, read_problem
from models.cp_solver import solve_vrp, solve_vrp_single_serve
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("week", type=int, help="week number (1 or 2)")
    parser.add_argument("--infile", type=str, help="init solution")
    parser.add_argument("--outfile", type=str, help="desired output file path")
    parser.add_argument(
        "--solver", type=str, help="solver to use (hexaly/ortools)", default="hexaly"
    )

    args = parser.parse_args()

    if not args.outfile:
        print("Please specify an output file")
        exit()

    problem = read_problem(
        "data/centres.xlsx",
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/euclidean_matrix.xlsx",
        week=args.week,
    )

    tours, obj = solve_vrp(problem, hint=args.infile, outfile=args.outfile)
