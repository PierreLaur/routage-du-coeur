import pandas as pd
import json
import sys
from problem import Problem, read_problem


def check_solution_file(
    centres_file,
    points_de_ramasse_file,
    vehicles_file,
    matrix_file,
    week,
    solution_file,
):
    sol = json.load(open(solution_file))
    total_distance = sol["total_distance"]
    tours = sol["tours"]

    pb = read_problem(
        centres_file,
        points_de_ramasse_file,
        vehicles_file,
        matrix_file,
        week,
    )

    obj = 0
    tours_flat = {}
    deliveries = {}
    palettes = {}
    for d in range(pb.n_days):
        for v in range(pb.m):
            key = str(d) + ", " + str(v)
            if key in tours:
                a = 0
                tours_flat[d, v] = [0]
                for place in tours[key]:
                    b = place["index"]
                    obj += pb.matrix.iloc[a, b]
                    tours_flat[d, v].append(b)
                    a = b

                    if b < pb.n:
                        deliveries[d, v, b] = place["delivery"]
                        palettes[d, v, b] = place["palettes"]

                obj += pb.matrix.iloc[a, 0]
                tours_flat[d, v].append(0)

                if not pb.frais[v]:
                    assert (
                        sum(
                            deliveries[d, v, c][1]
                            for c in tours_flat[d, v][1:-1]
                            if c < pb.n
                        )
                        == 0
                    )

                assert (
                    sum(
                        sum(deliveries[d, v, c])
                        for c in tours_flat[d, v][1:-1]
                        if c < pb.n
                    )
                    <= pb.capacities[v] + 0.001
                ), [
                    sum(
                        sum(deliveries[d, v, c])
                        for c in tours_flat[d, v][1:-1]
                        if c < pb.n
                    ),
                    pb.capacities[v],
                ]

                for c in tours_flat[d, v][1:-1]:
                    if c < pb.n:
                        assert (
                            deliveries[d, v, c][0] + deliveries[d, v, c][1]
                        ) <= pb.max_palette_capacity * palettes[d, v, c], [
                            deliveries[d, v, c][0] + deliveries[d, v, c][1],
                            palettes[d, v, c],
                        ]

                assert (
                    sum(palettes[d, v, c] for c in tours_flat[d, v][1:-1] if c < pb.n)
                    <= pb.sizes[v]
                )

    assert obj == int(total_distance), [obj, total_distance]

    for p in range(pb.n_pdr):
        for d in range(pb.n_days):
            node = pb.n + p
            n_visits = sum(
                [
                    tours_flat[d, v].count(node)
                    for v in range(pb.m)
                    if (d, v) in tours_flat
                ]
            )
            assert n_visits == pb.j_de_ramasse[p].count(d)
            if (d, p) in pb.use_pl:
                assert node in tours_flat[d, 0]

    for d in range(pb.n_days):
        for v in range(pb.m):
            if not (d, v) in tours_flat:
                continue
            assert len([p for p in tours_flat[d, v] if p > pb.n]) <= pb.sizes[v]
            assert (
                sum([pb.weights[p - pb.n] for p in tours_flat[d, v] if p > pb.n])
                <= pb.capacities[v]
            )

    for c in range(pb.n):
        assert (
            sum(
                deliveries[d, v, c][0]
                for v in range(pb.m)
                for d in range(pb.n_days)
                if (d, v, c) in deliveries
            )
            >= pb.demands["a"][c]
        )
        assert (
            sum(
                deliveries[d, v, c][1]
                for v in range(pb.m)
                for d in range(pb.n_days)
                if (d, v, c) in deliveries
            )
            >= pb.demands["f"][c]
        )
        assert (
            sum(
                deliveries[d, v, c][2]
                for v in range(pb.m)
                for d in range(pb.n_days)
                if (d, v, c) in deliveries
            )
            >= pb.demands["s"][c]
        )

    print("Assertions ok")


if __name__ == "__main__":
    file = sys.argv[1]
    week = int(sys.argv[2])

    check_solution_file(
        "data/centres.xlsx",
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/euclidean_matrix.xlsx",
        week,
        file,
    )
