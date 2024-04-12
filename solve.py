from utils.problem import Problem, Solution
from models.cp_solver import solve_vrp
import argparse
import os


def recompute_scenarios():
    print("todo")
    exit()
    for file in os.listdir("problems"):
        if "flex" not in file:
            continue
        if file.endswith(".json"):
            problem = Problem.from_json(f"problems/{file}")
            initsol = Solution.read_from_json(f"solutions/scenarios/{file}")
            solve_vrp(
                problem, initsol, outfile=f"solutions/scenarios/{file}", time_limit=10
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
