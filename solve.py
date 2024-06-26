from utils.check_solution import check_solution
from utils.plots import print_to_txt
from utils.problem import DeliveryWeek, Problem, Solution
from models.cp_solver import solve_day_per_day, solve_vrp  # noqa: F401
import argparse
import subprocess
from utils.evaluate_flexibility import evaluate_flexibility
from utils.prepare_demands import make_test_instances


def call_ortools(
    problem: Problem, week: DeliveryWeek, hint=None, time_limit=None, outfile=None
):
    status, solution = solve_vrp(
        problem,
        DeliveryWeek(args.week),
        hint=init_sol,
        time_limit=args.time_limit,
    )

    return solution


def call_hexaly(problem: str, week: int, hint=None, time_limit=None, outfile=None):
    command = [
        "localsolver",
        "models/ls_solver.lsp",
        problem,
        "-w",
        str(week),
    ]

    if hint:
        command += ["-i", hint]
    if time_limit:
        command += ["-t", str(time_limit)]
    if outfile:
        command += ["-o", outfile]

    try:
        subprocess.run(command)
    except KeyboardInterrupt:
        print("Solve interrupted")

    if outfile:
        try:
            solution = Solution.from_json(outfile)
        except FileNotFoundError:
            solution = None
        return solution


def validate_solution(problem_file: str, solution_file: str):
    problem = Problem.from_json(problem_file)
    solution = Solution.from_json(solution_file)
    check_solution(problem, solution)
    make_test_instances(100)
    evaluate_flexibility(problem_file, solution_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "problem_file", type=str, help="problem file (json, generated with problem.py)"
    )
    parser.add_argument(
        "--solver", "-s", type=str, help="solver (hexaly/ortools)", default="ortools"
    )
    parser.add_argument(
        "--week", "-w", type=int, help="week (1 or 2 for ODD/EVEN)", default=1
    )
    parser.add_argument(
        "--initsol", "-i", type=str, help="init solution (json solution file)"
    )
    parser.add_argument("--outfile", "-o", type=str, help="desired output file path")
    parser.add_argument("--time_limit", "-t", type=int, help="time limit")

    parser.add_argument(
        "--validate", "-v", type=str, help="a solution file to validate"
    )

    args = parser.parse_args()

    if args.validate:
        validate_solution(args.problem_file, args.validate)
        exit()

    problem = Problem.from_json(args.problem_file)

    if args.solver == "hexaly":
        solution = call_hexaly(
            args.problem_file,
            args.week,
            hint=args.initsol,
            time_limit=args.time_limit,
            outfile=args.outfile,
        )
    elif args.solver == "ortools":
        if args.initsol:
            init_sol = Solution.from_json(args.initsol)
        else:
            init_sol = None
        status, solution = solve_vrp(
            problem,
            DeliveryWeek(args.week),
            hint=init_sol,
            time_limit=args.time_limit,
        )
        # solution = solve_day_per_day(problem, DeliveryWeek(args.week))
    else:
        raise ValueError(f"Unknown solver: {args.solver}")

    if solution:
        solution.adjust_durations(problem)
        check_solution(problem, solution)
        print("[TEST] No constraint violations detected - the solution is valid.")
    else:
        print("No solution found")
        exit()

    if args.outfile:
        solution.to_json(args.outfile)

    if args.outfile:
        print_to_txt(solution, args.outfile.split(".json")[0] + ".txt")
