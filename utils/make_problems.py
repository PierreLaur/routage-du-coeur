import argparse
from utils.problem import Problem
import os
from itertools import combinations


def make_problem(demand_files, allowed_days_file, week_assignments_file, week, outfile):
    problem = Problem.from_files(
        demand_files,
        allowed_days_file,
        week_assignments_file,
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/matrices/distance_matrix.xlsx",
        "data/matrices/traffic_duration_matrix.xlsx",
        "data/matrices/no_traffic_duration_matrix.xlsx",
        "data/params.json",
        week,
    )

    problem.write_as_json(outfile)


def make_all(n_scenarios):
    allowed_days_root = "problems/allowed_days"
    week_assignments_root = "problems/week_assignments"
    demands_root = "problems/demands"
    for i, allowed_days_file in enumerate(os.listdir(allowed_days_root)):
        if not allowed_days_file.endswith(".xlsx"):
            continue

        for j, week_assignments_file in enumerate(os.listdir(week_assignments_root)):
            if not week_assignments_file.endswith(".xlsx"):
                continue

            demand_files = [
                os.path.join(demands_root, demand_file)
                for demand_file in os.listdir(demands_root)
                if demand_file.endswith(".xlsx")
            ]

            for k, files in enumerate(combinations(demand_files, n_scenarios)):
                for week in range(2):
                    outfile = f"problems/{n_scenarios}s_{i+1}{j+1}{k+1}_w{week+1}.json"
                    print("Making problem", outfile)
                    make_problem(
                        files,
                        os.path.join(allowed_days_root, allowed_days_file),
                        os.path.join(week_assignments_root, week_assignments_file),
                        week,
                        outfile,
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("n_scenarios", type=int)

    args = parser.parse_args()

    make_all(args.n_scenarios)
