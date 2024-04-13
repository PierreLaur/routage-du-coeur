import pytest
from ortools.sat.python import cp_model
from models.cp_solver import solve_vrp
from utils.problem import Problem
from utils.check_solution import check_solution
from tests.conftest import problem_files


@pytest.fixture(params=problem_files, scope="module")
def get_problem(request):
    problem = Problem.from_json(request.param)
    return problem


@pytest.fixture(scope="module")
def solve_instance(get_problem):
    time_limit = 30
    status, solution = solve_vrp(
        get_problem, time_limit=time_limit, stop_at_first_solution=True
    )
    return status, solution


def test_hint(make_test_problem, make_test_solution):
    time_limit = 5
    solution = make_test_solution
    status, new_solution = solve_vrp(
        make_test_problem,
        hint=solution,
        time_limit=time_limit,
        stop_at_first_solution=True,
        fix_hint=True,
    )
    assert status in [cp_model.FEASIBLE, cp_model.OPTIMAL]
    assert new_solution is not None
    assert new_solution.total_costs <= solution.total_costs


def test_finds_feasible(solve_instance):
    status, solution = solve_instance
    assert status in [cp_model.FEASIBLE, cp_model.OPTIMAL]
    assert solution is not None


def test_valid_solution(get_problem, solve_instance):
    status, solution = solve_instance
    if not solution:
        pytest.skip("no solution")
    check_solution(get_problem, solution)
