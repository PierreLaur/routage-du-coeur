import os
import pytest

from utils.problem import Problem
from models.cp_solver import solve_vrp

problem_files = [
    os.path.join("problems", f)
    for f in [
        # "1s_111_w1.json",
        # "1s_111_w2.json",
        # "1s_113_w1.json",
        # "1s_113_w2.json",
        "1s_114_w1.json",
        "1s_114_w2.json",
        # "1s_116_w1.json",
        # "1s_116_w2.json",
    ]
]


@pytest.fixture(scope="package")
def make_test_problem():
    path = "problems/example.json"
    problem = Problem.from_json(path)
    return problem


@pytest.fixture(scope="package")
def make_test_solution(make_test_problem):
    path = "solutions/example.json"
    status, solution = solve_vrp(make_test_problem, hint=None, time_limit=3)
    if not solution:
        raise Exception("No solution found to example problem")
    solution.to_json(path)
    return solution
