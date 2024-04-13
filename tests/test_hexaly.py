import pytest
from utils.problem import Problem, Solution
from utils.check_solution import check_solution
from tests.conftest import problem_files
import subprocess
import os


@pytest.fixture(params=problem_files, scope="module")
def solve_instance(request):
    time_limit = 2
    tmp_file = "solutions/tmp.json"
    try:
        os.remove(tmp_file)
    except OSError:
        pass
    subprocess.run(
        [
            "localsolver",
            "models/ls_solver.lsp",
            request.param,
            "nil",
            tmp_file,
            f"{time_limit}",
        ]
    )
    yield request.param, tmp_file
    try:
        os.remove(tmp_file)
    except OSError:
        pass


def test_finds_feasible(solve_instance):
    _, solution_file = solve_instance
    assert os.path.isfile(solution_file)


def test_valid_solution(solve_instance):
    problem_file, solution_file = solve_instance

    if not os.path.isfile(solution_file):
        pytest.skip("no solution")

    problem = Problem.from_json(problem_file)
    solution = Solution.from_json(solution_file)
    check_solution(problem, solution)


# @pytest.mark.parametrize("total_costs", [2000, 1500, 1000, 900, 800, 700, 600])
# def test_solution_quality(solve_instance, total_costs):
#     problem_file, solution_file = solve_instance
#     solution = Solution.from_json(solution_file)
#     quality = solution.total_costs
#     assert quality <= total_costs
