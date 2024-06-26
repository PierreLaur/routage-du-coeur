import pytest
from utils.problem import Problem, Solution
from utils.check_solution import check_solution
from tests.conftest import problem_files, current_problem, current_solution
import subprocess
import os


@pytest.fixture(params=problem_files, scope="module")
def solve_instance(request, get_week):
    time_limit = 3
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
            "-w",
            str(get_week.value),
            "-o",
            tmp_file,
            "-t",
            str(time_limit),
        ]
    )
    return request.param, tmp_file


def test_finds_feasible(solve_instance):
    _, solution_file = solve_instance
    assert os.path.isfile(solution_file)


def test_valid_solution(solve_instance, record_property):
    problem_file, solution_file = solve_instance

    if not os.path.isfile(solution_file):
        pytest.skip("no solution")

    problem = Problem.from_json(problem_file)
    solution = Solution.from_json(solution_file)
    check_solution(problem, solution)
    record_property("score", solution.total_costs)

    try:
        os.remove(solution_file)
    except OSError:
        pass


def test_improves_over_current_solution(record_property):
    time_limit = 20
    tmp_file = "solutions/tmp.json"
    if not os.path.isfile(current_solution):
        pytest.skip("current solution not found")
    cur_solution = Solution.from_json(current_solution)
    try:
        os.remove(tmp_file)
    except OSError:
        pass
    subprocess.run(
        [
            "localsolver",
            "models/ls_solver.lsp",
            current_problem,
            "-o",
            tmp_file,
            "-t",
            str(time_limit),
        ]
    )

    if not os.path.isfile(tmp_file):
        pytest.skip("no solution")
    solution = Solution.from_json(tmp_file)
    assert solution.total_costs < cur_solution.total_costs
    record_property(
        "margin",
        (cur_solution.total_costs - solution.total_costs) / cur_solution.total_costs,
    )


# @pytest.mark.parametrize("total_costs", [2000, 1500, 1000, 900, 800, 700, 600])
# def test_solution_quality(solve_instance, total_costs):
#     problem_file, solution_file = solve_instance
#     solution = Solution.from_json(solution_file)
#     quality = solution.total_costs
#     assert quality <= total_costs
