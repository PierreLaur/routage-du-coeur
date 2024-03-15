from utils.problem import Problem, Solution
from models.cp_solver import solve_vrp
import argparse
import os


def compute_scenarios(centre_files, scenario_names, time_limit=300):
    for file, name in zip(centre_files, scenario_names):
        for week in [1, 2]:
            print("- - - - - Computing scenario", name, "for week", week)

            outfile_name = f"solutions/scenarios/{name}_{week}.json"

            problem = Problem.from_files(
                file,
                "data/points_de_ramasse.xlsx",
                "data/vehicules.xlsx",
                "data/distance_matrix.xlsx",
                "data/traffic_duration_matrix.xlsx",
                "data/no_traffic_duration_matrix.xlsx",
                "data/params.json",
                week=week,
            )

            if os.path.exists(outfile_name):
                init_sol = Solution.read_from_json(outfile_name)
            else:
                init_sol = None

            solve_vrp(
                problem,
                hint=init_sol,
                outfile=outfile_name,
                violation_cost=None,
                time_limit=time_limit,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "problem_file", type=str, help="problem file (json, generated with problem.py)"
    )
    parser.add_argument(
        "--initsol", type=str, help="init solution (json solution file)"
    )
    parser.add_argument("--outfile", type=str, help="desired output file path")
    parser.add_argument("--time_limit", "-t", type=int, help="time limit")

    args = parser.parse_args()

    problem = Problem.from_json(args.problem_file)

    if args.initsol:
        init_sol = Solution.read_from_json(args.initsol)
    else:
        init_sol = None

    solve_vrp(
        problem,
        hint=init_sol,
        outfile=args.outfile,
        time_limit=args.time_limit,
    )
