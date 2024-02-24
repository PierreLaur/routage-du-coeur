from utils.problem import Problem, Solution
from models.cp_solver import solve_vrp
from models.routing_solver import solve
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
        init_sol = Solution.read_from_json(args.improve)
    else:
        outfile = args.outfile
        if args.infile:
            init_sol = Solution.read_from_json(args.infile)
        else:
            init_sol = None

    centres_file = "data/centres_variations/centres_keep.xlsx"

    problem = Problem.from_files(
        centres_file,
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/distance_matrix.xlsx",
        "data/traffic_duration_matrix.xlsx",
        "data/no_traffic_duration_matrix.xlsx",
        "data/params.json",
        week=args.week,
    )

    solve_vrp(
        problem,
        hint=init_sol,
        outfile=outfile,
        violation_cost=None,
    )
