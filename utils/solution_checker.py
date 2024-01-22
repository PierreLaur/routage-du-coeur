import json
import sys
from problem import Problem, read_problem
import pandas as pd


def check_solution_file(
    centres_file,
    points_de_ramasse_file,
    vehicles_file,
    matrix_file,
    params_file,
    week,
    solution_file,
):
    """Checks whether a solution file satisfies all constraints"""
    sol = json.load(open(solution_file))
    total_distance = sol["total_distance"]
    fuel_consumption = sol["fuel_consumption"]
    tours = sol["tours"]

    pb = read_problem(
        centres_file,
        points_de_ramasse_file,
        vehicles_file,
        matrix_file,
        params_file,
        week,
    )

    duration_matrix = pd.read_excel("data/duration_matrix.xlsx", index_col=0)
    warning_duration_threshold = (
        pb.max_tour_duration
        - 30  # Show a warning if estimated duration for a tour is higher than this
    )

    obj = 0
    distance = 0
    tours_flat = {}
    deliveries = {}
    palettes = {}
    norvegiennes = {}
    for d in range(pb.n_days):
        for v in range(pb.m):
            key = str(d) + ", " + str(v)
            if key in tours:
                a = 0
                tours_flat[d, v] = [0]
                for place in tours[key]:
                    b = place["index"]
                    distance += pb.matrix.iloc[a, b]
                    obj += pb.matrix.iloc[a, b] * pb.consumptions[v]
                    tours_flat[d, v].append(b)
                    a = b

                    if b < pb.n:
                        deliveries[d, v, b] = place["delivery"]
                        palettes[d, v, b] = place["palettes"]
                        norvegiennes[d, v, b] = place["norvegiennes"]

                distance += pb.matrix.iloc[a, 0]
                obj += pb.matrix.iloc[a, 0] * pb.consumptions[v]
                tours_flat[d, v].append(0)

                # Only deliver frais with camions frigos
                if not pb.frais[v]:
                    assert (
                        sum(
                            deliveries[d, v, c][1]
                            for c in tours_flat[d, v][1:-1]
                            if c < pb.n
                        )
                        == 0
                    )

                # Only deliver surgelé on pallets with camions frigos
                if not pb.frais[v]:
                    assert (
                        sum(
                            palettes[d, v, c][2]
                            for c in tours_flat[d, v][1:-1]
                            if c < pb.n
                        )
                        == 0
                    )

                # Check the vehicle capacities
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
                        # Use enough palettes for each centre
                        assert (
                            deliveries[d, v, c][0]
                            <= pb.max_palette_capacity * palettes[d, v, c][0]
                        )
                        assert (
                            deliveries[d, v, c][1]
                            <= pb.demi_palette_capacity * palettes[d, v, c][1]
                        )

                        # use either palettes, norvegiennes or both but enough of them
                        assert (
                            deliveries[d, v, c][2]
                            <= pb.demi_palette_capacity * palettes[d, v, c][2]
                            + norvegiennes[d, v, c] * pb.norvegienne_capacity
                        ), [
                            deliveries[d, v, c][2],
                            pb.demi_palette_capacity * palettes[d, v, c][2],
                            norvegiennes[d, v, c] * pb.norvegienne_capacity,
                        ]

                # Size limit for each vehicle
                assert (
                    sum(
                        palettes[d, v, c][0]
                        + 0.5 * (palettes[d, v, c][1] + palettes[d, v, c][2])
                        for c in tours_flat[d, v][1:-1]
                        if c < pb.n
                    )
                    <= pb.sizes[v]
                )

                estimated_duration = sum(
                    duration_matrix.iloc[a, b]
                    for a, b in zip(tours_flat[d, v][:-2], tours_flat[d, v][1:-1])
                ) / 60 + sum(
                    pb.wait_at_centres if c < pb.n else pb.wait_at_pdrs
                    for c in tours_flat[d, v][1:-1]
                )

                if estimated_duration > warning_duration_threshold:
                    print(
                        f"Warning : tour {key} is ~{estimated_duration//60:.0f}h{estimated_duration%60:.0f}min long :",
                        [duration_matrix.index[c] for c in tours_flat[d, v][1:-1]],
                    )

    for d in range(pb.n_days):
        # Don't use too many norvegiennes
        assert (
            sum(
                norvegiennes[d, v, c]
                for v in range(pb.m)
                if (d, v) in tours_flat
                for c in tours_flat[d, v][1:-1]
                if c < pb.n
            )
            <= pb.n_norvegiennes
        )

    # Check that the objective value is correct
    assert obj / 100000 == fuel_consumption
    assert distance == int(total_distance), [distance, total_distance]

    # Check that the number of pickups is correct for each day
    for p in range(pb.n_pdr):
        node = pb.n + p
        for d in range(pb.n_days):
            n_visits = sum(
                [
                    tours_flat[d, v].count(node)
                    for v in range(pb.m)
                    if (d, v) in tours_flat and pb.frais[v]
                ]
            )
            assert n_visits == pb.j_de_ramasse[p].count(d), [
                p,
                d,
                n_visits,
                pb.j_de_ramasse[p].count(d),
            ]

    # Specific requirements
    assert (2, 0) in tours_flat
    assert (2, 2) in tours_flat
    assert (4, 2) in tours_flat
    carrefour_centrale_index = pb.n + pb.n_pdr - 1

    ### Carrefour centrale must be visited by these vehicles on these days
    assert carrefour_centrale_index in tours_flat[2, 0][1:-1]
    assert carrefour_centrale_index in tours_flat[2, 2][1:-1]
    assert carrefour_centrale_index in tours_flat[4, 2][1:-1]

    # Check time window constraints
    for d in range(pb.n_days):
        for v in range(pb.m):
            if not (d, v) in tours_flat:
                continue
            for node in tours_flat[d, v][1:-1]:
                if node < pb.n:
                    assert d in pb.j_de_livraison_possibles[node]

    for d in range(pb.n_days):
        for v in range(pb.m):
            if not (d, v) in tours_flat:
                continue
            # Check that not too many pickups are done - each takes 2 palettes
            assert len([p for p in tours_flat[d, v] if p > pb.n]) * 2 <= pb.sizes[v]
            # Capacity constraints on vehicles for the pickups
            assert (
                sum([pb.weights[p - pb.n] for p in tours_flat[d, v] if p > pb.n])
                <= pb.capacities[v]
            )

    # Check that demands are met
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

    print("No constraint violations detected - the solution is valid.")


if __name__ == "__main__":
    file = sys.argv[1]
    week = int(sys.argv[2])

    check_solution_file(
        "data/centres.xlsx",
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/euclidean_matrix.xlsx",
        "data/params.json",
        week,
        file,
    )
