def check_solution(centres, vehicles, tours, obj, deliveries, visits, arcs):
    """TODO : update"""
    n = len(centres.index)
    m = len(vehicles.index)
    n_days = 5

    demands = {}
    demands["a"] = centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist()
    demands["f"] = centres["Tonnage Frais (kg)"].fillna(0).astype(int).tolist()
    demands["s"] = centres["Tonnage Surgelé (kg)"].fillna(0).astype(int).tolist()

    for c in range(n):
        demands["a"][c] = int(demands["a"][c] * 1.15)
        demands["f"][c] = int(demands["f"][c] * 1.15)
        demands["s"][c] = int(demands["s"][c] * 1.15)

    capacities = vehicles["Capacité (kg)"].astype(int).tolist()
    sizes = vehicles["Taille(Palettes)"].astype(int).tolist()

    da, df, ds = deliveries

    for v in range(m):
        for d in range(n_days):
            for c in range(n):
                vd = v + d * m

                if visits[v, d, c]:
                    assert sum([arcs[vd, c2, c] for c2 in range(n) if c2 != c]) == 1
                    assert arcs[vd, c, c] == 0
                    if da[v, d, c] + df[v, d, c] + ds[v, d, c] == 0:
                        f"Warning : day {d} vehicle {v} visits centre {c} with no deliveries"
                else:
                    assert sum([arcs[vd, c2, c] for c2 in range(n) if c2 != c]) == 0
                    assert arcs[vd, c, c] == 1
                    assert da[v, d, c] + df[v, d, c] + ds[v, d, c] == 0

            assert sum([da[v, d, c] + df[v, d, c] for c in range(n)]) <= capacities[v]

            if vehicles["Frais"].iloc[v] == "Non":
                assert sum([df[v, d, c] for c in range(n)]) == 0

    # TODO : check palettes

    for c in range(1, n):
        assert (
            sum([da[v, d, c] for v in range(m) for d in range(n_days)])
            == demands["a"][c]
        )
        assert (
            sum([df[v, d, c] for v in range(m) for d in range(n_days)])
            == demands["f"][c]
        )
        assert (
            sum([ds[v, d, c] for v in range(m) for d in range(n_days)])
            == demands["s"][c]
        )

    print("Assertions ok")
