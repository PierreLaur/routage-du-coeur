import os
import pytest
from utils.problem import DeliveryWeek, Problem
from models.cp_solver import solve_vrp

problem_files = [
    os.path.join("problems", f)
    for f in [
        "fixe_usual_original.json",
        "fixe_usual_0226.json",
        "fixe_usual_median.json",
        # "fixe_0226_original.json",
        # "fixe_0226_0226.json",
        # "fixe_0226_median.json",
    ]
]

current_problem = "problems/current_problem.json"
current_solution = "solutions/current_tours/current_problem.json"


@pytest.fixture(scope="package")
def make_test_problem():
    path = "problems/example.json"
    problem = Problem.from_json(path)
    return problem


@pytest.fixture(params=(DeliveryWeek.ODD, DeliveryWeek.EVEN), scope="package")
def get_week(request):
    return request.param


@pytest.fixture(scope="package")
def make_test_solution(make_test_problem):
    path = "solutions/example.json"
    status, solution = solve_vrp(
        make_test_problem, week=DeliveryWeek.ODD, hint=None, time_limit=3
    )
    if not solution:
        raise Exception("No solution found to example problem")
    solution.to_json(path)
    return solution


def pytest_report_teststatus(report):
    if report.when == "call" and report.passed:
        if "margin" in dict(report.user_properties):
            margin = dict(report.user_properties)["margin"]
            short_outcome = "."
            long_outcome = f"PASSED with {margin * 100:.2f}% improvement"
            return report.outcome, short_outcome, long_outcome

        if "score" in dict(report.user_properties):
            score = dict(report.user_properties)["score"]
            short_outcome = "."
            long_outcome = f"PASSED with cost={score:.2f}€"
            return report.outcome, short_outcome, long_outcome
