import argparse
from utils.check_solution import check_solution
from utils.evaluate_flexibility import evaluate_flexibility
from utils.prepare_demands import make_test_instances
from utils.problem import Problem, Solution


def validate_solution(problem_file: str, solution_file: str):
    problem = Problem.from_json(problem_file)
    solution = Solution.from_json(solution_file)
    check_solution(problem, solution)
    print("[TEST] No constraint violations detected - the solution is valid.")
    make_test_instances(100)
    evaluate_flexibility(problem_file, solution_file, solution.week.value)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "problem_file", type=str, help="problem file (json, generated with problem.py)"
    )

    parser.add_argument("solution_file", type=str, help="a solution file to validate")

    args = parser.parse_args()
    validate_solution(args.problem_file, args.solution_file)
