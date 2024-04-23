from utils.problem import Problem
import os


def make_problem(demand_files, allowed_days_file, week_assignments_file):
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
    )
    return problem


def make_all():
    allowed_days_root = "problems/allowed_days"
    week_assignments_root = "problems/week_assignments"
    demands_root = "problems/demands"
    for i, allowed_days_file in enumerate(os.listdir(allowed_days_root)):
        if not allowed_days_file.endswith(".xlsx"):
            continue
        allowed_days_name = allowed_days_file.split(".")[0]

        for j, week_assignments_file in enumerate(os.listdir(week_assignments_root)):
            if not week_assignments_file.endswith(".xlsx"):
                continue
            week_assignments_name = week_assignments_file.split(".")[0]

            demand_files = [
                os.path.join(demands_root, demand_file)
                for demand_file in os.listdir(demands_root)
                if demand_file.endswith(".xlsx")
            ]

            for file in demand_files:
                name = os.path.basename(file).split(".")[0]
                outfile = (
                    f"problems/{allowed_days_name}_{week_assignments_name}_{name}.json"
                )
                print("Making problem", outfile)
                problem = make_problem(
                    file,
                    os.path.join(allowed_days_root, allowed_days_file),
                    os.path.join(week_assignments_root, week_assignments_file),
                )
                problem.to_json(outfile)


if __name__ == "__main__":
    make_all()
