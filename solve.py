from utils.problem import Problem, read_problem
from models.cp_solver import solve_vrp, cluster_nodes
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("week", type=int, help="week number (1 or 2)")
    parser.add_argument("--infile", type=str, help="init solution")
    parser.add_argument("--outfile", type=str, help="desired output file path")
    parser.add_argument("--improve", type=str, help="desired file to improve")

    args = parser.parse_args()
    assert args.week in [1, 2]

    if args.improve:
        outfile = args.improve
        infile = args.improve
    else:
        outfile = args.outfile
        infile = args.infile

    if not outfile:
        print("Please specify an output file")
        exit()

    centres_file = "data/centres.xlsx"

    problem = read_problem(
        centres_file,
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/euclidean_matrix.xlsx",
        "data/duration_matrix_w_traffic.xlsx",
        "data/params.json",
        week=args.week,
    )

    # visits = cluster_nodes(problem)

    tours, obj = solve_vrp(
        problem,
        hint=infile,
        outfile=outfile,
        violation_cost=None,
        predefined_visits=None,
    )
