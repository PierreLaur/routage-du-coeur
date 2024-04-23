import os
from utils.problem import Problem, Solution
from utils.plots import print_to_txt, plot_solution


def test_write_read_problem(make_test_problem):
    tmp_file = "problems/tmp.json"
    problem = make_test_problem
    problem.to_json(tmp_file)
    assert os.path.isfile(tmp_file)
    new_problem = Problem.from_json(tmp_file)
    assert new_problem == problem
    os.remove(tmp_file)


def test_write_read_solution(make_test_solution):
    solution = make_test_solution
    tmp_file = "solutions/tmp.json"
    solution.to_json(tmp_file)
    assert os.path.isfile(tmp_file)
    new_solution = Solution.from_json(tmp_file)
    assert new_solution == solution
    os.remove(tmp_file)


def test_to_txt(make_test_solution):
    solution = make_test_solution
    path = "solutions/example.txt"
    try:
        os.remove(path)
    except OSError:
        pass
    print_to_txt(solution, path)
    assert os.path.isfile(path)


def test_plot_solution(make_test_solution):
    solution = make_test_solution
    path = "solutions/example.html"
    try:
        os.remove(path)
    except OSError:
        pass
    plot_solution(solution, path)
    assert os.path.isfile(path)
