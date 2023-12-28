import pandas as pd
import json


def check_solution_file(
    matrix_file, centres_file, points_de_ramasse_file, vehicles_file, solution_file
):
    centres = pd.read_excel(centres_file, index_col=0)
    points_de_ramasse = pd.read_excel(points_de_ramasse_file, index_col=0)
    vehicles = pd.read_excel(vehicles_file, index_col=0)
    matrix = pd.read_excel(matrix_file, index_col=0)

    sol = json.load(open(solution_file))
    total_distance = sol["total_distance"]
    tours = sol["tours"]

    n = len(centres.index)
    n_pdr = len(points_de_ramasse.index)
    n_days = 5  # Single week scheduling
    m = len(vehicles.index)

    capacities = vehicles["Capacité (kg)"].astype(int).tolist()
    sizes = vehicles["Taille(Palettes)"].astype(int).tolist()
    max_palette_capacity = 800

    demands = {
        "a": centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist(),
        "f": centres["Tonnage Frais (kg)"].fillna(0).astype(int).tolist(),
        "s": centres["Tonnage Surgelé (kg)"].fillna(0).astype(int).tolist(),
    }

    ignore_centres = ["Cugnaux", "Levignac", "Rieumes"]

    # Ignore some centres (deliver them the other week) :
    indexes = [centres.loc[centres["Nom"] == c].index.values[0] for c in ignore_centres]
    for i in indexes:
        demands["a"][i] = 0
        demands["f"][i] = 0
        demands["s"][i] = 0

    # Add 15% for robustness
    for c in range(1, n):
        demands["a"][c] = int(demands["a"][c] * 1.15)
        demands["f"][c] = int(demands["f"][c] * 1.15)
        demands["s"][c] = int(demands["s"][c] * 1.15)

    for d in demands.values():
        d[0] = 0

    freqs_pdr = (
        points_de_ramasse["Fréquence de Ramasse(/w)"].fillna(0).astype(int).tolist()
    )

    frais = [f == "Oui" for f in vehicles["Frais"]]

    obj = 0
    tours_flat = {}
    deliveries = {}
    palettes = {}
    for d in range(n_days):
        for v in range(m):
            key = str(d) + ", " + str(v)
            if key in tours:
                a = 0
                tours_flat[d, v] = [0]
                for place in tours[key]:
                    b = place["index"]
                    obj += matrix.iloc[a, b]
                    tours_flat[d, v].append(b)
                    a = b

                    if b < n:
                        deliveries[d, v, b] = place["delivery"]
                        palettes[d, v, b] = place["palettes"]

                obj += matrix.iloc[a, 0]
                tours_flat[d, v].append(0)

                if not frais[v]:
                    assert (
                        sum(
                            deliveries[d, v, c][1]
                            for c in tours_flat[d, v][1:-1]
                            if c < n
                        )
                        == 0
                    )

                assert (
                    sum(
                        sum(deliveries[d, v, c])
                        for c in tours_flat[d, v][1:-1]
                        if c < n
                    )
                    <= capacities[v] + 0.001
                ), [
                    sum(
                        sum(deliveries[d, v, c])
                        for c in tours_flat[d, v][1:-1]
                        if c < n
                    ),
                    capacities[v],
                ]

                for c in tours_flat[d, v][1:-1]:
                    if c < n:
                        assert (
                            deliveries[d, v, c][0] + deliveries[d, v, c][1]
                        ) <= max_palette_capacity * palettes[d, v, c], [
                            deliveries[d, v, c][0] + deliveries[d, v, c][1],
                            palettes[d, v, c],
                        ]

                assert (
                    sum(palettes[d, v, c] for c in tours_flat[d, v][1:-1] if c < n)
                    <= sizes[v]
                )

    assert obj == int(total_distance), [obj, total_distance]

    for p in range(n_pdr):
        node = n + p
        n_visits = sum([tour.count(node) for tour in tours_flat.values()])
        assert n_visits >= freqs_pdr[p]

    for c in range(n):
        assert (
            sum(
                deliveries[d, v, c][0]
                for v in range(m)
                for d in range(n_days)
                if (d, v, c) in deliveries
            )
            >= demands["a"][c]
        )
        assert (
            sum(
                deliveries[d, v, c][1]
                for v in range(m)
                for d in range(n_days)
                if (d, v, c) in deliveries
            )
            >= demands["f"][c]
        )
        assert (
            sum(
                deliveries[d, v, c][2]
                for v in range(m)
                for d in range(n_days)
                if (d, v, c) in deliveries
            )
            >= demands["s"][c]
        )

    print("Assertions ok")


check_solution_file(
    "data/euclidean_matrix.xlsx",
    "data/centres.xlsx",
    "data/points_de_ramasse.xlsx",
    "data/vehicules.xlsx",
    "solutions/test2.json",
)
