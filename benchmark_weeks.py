import pandas as pd
from itertools import product
from subprocess import Popen, PIPE
from time import time
from random import shuffle


def bench_assignments():
    """Centres that are delivered every other week must be assigned to either week 1 or week 2
    This function generates all possible assignments of these centres, runs the solver for each one, and generates a
    assignments.csv file with the results"""

    centres = pd.read_csv("data/centres_original.csv")

    semi_hebdo = list()
    for i in range(len(centres)):
        if centres["Semaine"][i] != 0:
            semi_hebdo.append(i)

    assignments = list(product([1, 2], repeat=len(semi_hebdo)))[:128]
    shuffle(assignments)

    res_week1 = []
    res_week2 = []
    total = []
    for index, a in enumerate(assignments):
        centres_new = pd.read_csv("data/centres_original.csv")
        week_new = [0] * len(centres)

        for i, j in zip(semi_hebdo, a):
            week_new[i] = j

        centres_new["Semaine"] = week_new
        centres_new.to_csv("data/centres.csv", index=False)

        res = {}
        print(f"\n- - - - EVALUATING ASSIGNMENT {index+1}/{len(assignments)} - - - -\n")

        print(
            f"Semaine 1 : ",
            ", ".join(
                [centres["Nom"][i] for ind, i in enumerate(semi_hebdo) if a[ind] == 1]
            ),
        )
        print(
            f"Semaine 2 : ",
            ", ".join(
                [centres["Nom"][i] for ind, i in enumerate(semi_hebdo) if a[ind] == 2]
            ),
        )

        for week in [1, 2]:
            print("\n- - - - - - - - WEEK ", week, "\n")

            file_name = f"overnight/assignment{index}week{week}.json"
            process = Popen(
                [
                    "localsolver",
                    "models/ls_solver.lsp",
                    # f"solutions/week_{week}.json",
                    "nil",
                    file_name,
                    f"{week}",
                    "18",
                ],
                stderr=PIPE,
                stdout=PIPE,
            )

            st = time()
            while time() - st < 20:
                output = process.stdout.readline()
                if output == "" and process.poll() is not None:
                    break
                if output:
                    line = output.decode("utf-8").strip()
                    if "%" not in line:
                        print(line)

            res[week] = float(line.split("Total fuel consumption : ")[1])

        res_week1.append(res[1])
        res_week2.append(res[2])
        total.append(sum(res.values()))

    assignments_df = pd.DataFrame(
        pd.Series(name="assignments", data=assignments[: len(res_week1)])
    )
    assignments_df["week1"] = res_week1
    assignments_df["week2"] = res_week2
    assignments_df["total"] = total

    assignments_df.to_csv("data/assignments.csv")

    centres.to_csv("data/centres.csv", index=False)


if __name__ == "__main__":
    bench_assignments()
