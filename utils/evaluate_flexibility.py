import argparse
import subprocess


def evaluate_flexibility(problem: str, solution):
    command = [
        "localsolver",
        "models/ls_solver.lsp",
        problem,
        "-e",
        solution,
    ]

    try:
        subprocess.run(command)
    except KeyboardInterrupt:
        print("Solve interrupted")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "problem_file", type=str, help="problem file (json, generated with problem.py)"
    )

    parser.add_argument("solution", type=str, help="the solution to evaluate")

    args = parser.parse_args()

    evaluate_flexibility(args.problem_file, args.solution)
