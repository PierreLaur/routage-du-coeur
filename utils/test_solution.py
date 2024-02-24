import json
import sys
from .problem import Problem, Solution, ProductType, DeliveryWeek, Stop, StopType
import pandas as pd
from numpy import int64


def check_solution_file(
    centres_file,
    points_de_ramasse_file,
    vehicles_file,
    dist_matrix_file,
    duration_matrix_file,
    no_traffic_duration_matrix_file,
    params_file,
    solution_file_path,
):
    """Checks whether a solution file satisfies all constraints"""
    sol = Solution.read_from_json(solution_file_path)

    pb = Problem.from_files(
        centres_file,
        points_de_ramasse_file,
        vehicles_file,
        dist_matrix_file,
        duration_matrix_file,
        no_traffic_duration_matrix_file,
        params_file,
        sol.week,
    )

    no_traffic_duration_matrix = pd.read_excel(
        "data/duration_matrix.xlsx", index_col=0
    ).astype(int)

    warning_duration_threshold = (
        200  # Show a warning if estimated duration for a tour is higher than this
    )

    actual_dist = 0
    for (d, v), tour in sol.tours.items():
        tour_indices = [stop.index for stop in tour]
        for a, b in zip([0] + tour_indices, tour_indices + [0]):
            new_dist = pb.distance_matrix.iloc[a, b]
            if type(new_dist) != int64:
                continue
            actual_dist += new_dist

        estimated_duration = 0
        for a, b in zip([0] + tour_indices, tour_indices + [0]):
            if estimated_duration <= 90:
                estimated_duration += pb.duration_matrix.iloc[a, b] / 60  # type: ignore
            else:
                estimated_duration += no_traffic_duration_matrix.iloc[a, b] / 60  # type: ignore

        estimated_duration += sum(
            pb.params.wait_at_centres * (stop.type == StopType.Livraison)
            + pb.params.wait_at_pdrs * (stop.type == StopType.Ramasse)
            for stop in tour
        )

        if estimated_duration > warning_duration_threshold:
            print(
                f"Warning : tour {d, v} is ~{estimated_duration//60:.0f}h{estimated_duration%60:.0f}min long :",
                [stop.name for stop in tour],
            )

        # first_pickup_estimated_duration = (
        #     sum(
        #         pb.duration_matrix.iloc[a, b]
        #         for a, b in zip(tours_flat[d, v][:-2], tours_flat[d, v][1:-1])
        #         if a == 0 or a < pb.n_centres
        #                             )
        #     / 60
        #     + sum(
        #         pb.params.wait_at_centres for c in tours_flat[d, v][1:-1] if c < pb.n_centres
        #                             )
        #     if any(node >= pb.n for node in tours_flat[d, v])
        #     else 0
        # )

        # if first_pickup_estimated_duration > pb.max_first_pickup_time:
        #     print(
        #         f"Warning : first pickup of tour {key} is at ~{first_pickup_estimated_duration//60:.0f}h{first_pickup_estimated_duration%60:.0f}min :",
        #         [pb.duration_matrix.index[c] for c in tours_flat[d, v][1:-1]],
        #     )

    # Check that the objective values are correct
    assert actual_dist == int(sol.total_distance)

    for (d, v), tour in sol.tours.items():
        for stop in tour:
            if stop.type == StopType.Livraison:

                # Product type requirements
                assert all(
                    ProductType(t) in pb.vehicles[v].can_carry
                    for t in range(3)
                    if stop.delivery[t] > 0
                )
                assert all(
                    stop.norvegiennes == 0
                    for (_, v), tour in sol.tours.items()
                    if ProductType.S not in pb.vehicles[v].can_carry
                    for stop in tour
                )

                # Only deliver palettes of S products if the vehicle allows it
                if not pb.vehicles[v].allows_isotherm_cover:
                    assert stop.palettes[2] == 0

                # Use enough palettes for each centre
                assert (
                    stop.delivery[0]
                    <= pb.params.max_palette_capacity * stop.palettes[0]
                )

                # Number of A palettes must be int
                assert stop.palettes[0] % 1 == 0

                full, rest = divmod(stop.palettes[1], 1)
                halves = rest * 2
                assert (
                    stop.delivery[1]
                    <= pb.params.demi_palette_capacity * halves
                    + pb.params.max_palette_capacity * full
                )

                # use either palettes, norvegiennes or both but enough of them
                full, rest = divmod(stop.palettes[2], 1)
                halves = rest * 2
                assert (
                    stop.delivery[2]
                    <= pb.params.demi_palette_capacity * halves
                    + pb.params.max_palette_capacity * full
                    + stop.norvegiennes * pb.params.norvegienne_capacity
                )

        # Check the vehicle capacities
        assert (
            sum(sum(stop.delivery) for stop in tour if stop.type == StopType.Livraison)
            <= pb.vehicles[v].capacity
        )
        assert (
            sum(
                pb.pdrs[stop.index - pb.n_centres].weight
                for stop in tour
                if stop.type == StopType.Ramasse
            )
            <= pb.vehicles[v].capacity
        )

        # Size limit for each vehicle
        assert (
            sum(sum(stop.palettes) for stop in tour if stop.type == StopType.Livraison)
            <= pb.vehicles[v].size
        )
        # Each ramasse takes 2 palettes
        assert (
            sum(stop.type == StopType.Ramasse for stop in tour) * 2
            <= pb.vehicles[v].size
        )

    for d in range(pb.n_days):
        # Don't use too many norvegiennes
        assert (
            sum(
                stop.norvegiennes
                for (d2, v), tour in sol.tours.items()
                for stop in tour
                if d2 == d
            )
            <= pb.params.n_norvegiennes
        )

    # Check that the number of pickups with appropriate vehicles is correct for each day
    for p in range(pb.n_pdr):
        node = pb.n_centres + p
        for d in range(pb.n_days):
            n_visits = sum(
                stop.index == node
                for (d2, v), tour in sol.tours.items()
                if d2 == d and pb.pdrs[p].product_type in pb.vehicles[v].can_carry
                for stop in tour
            )

            assert n_visits == (d in pb.pdrs[p].required_days)

    # Specific requirements
    assert (2, 0) in sol.tours
    assert (2, 2) in sol.tours
    assert (4, 2) in sol.tours

    ### Carrefour centrale must be visited by these vehicles on these days
    assert "Carrefour Supply Chain (A)" in [node.name for node in sol.tours[2, 0]]
    assert "Carrefour Supply Chain" in [node.name for node in sol.tours[2, 2]]
    assert "Carrefour Supply Chain" in [node.name for node in sol.tours[4, 2]]

    ### No other pickups
    assert sum(stop.type == StopType.Ramasse for stop in sol.tours[2, 0]) == 1
    assert sum(stop.type == StopType.Ramasse for stop in sol.tours[4, 2]) == 1
    assert sum(stop.type == StopType.Ramasse for stop in sol.tours[2, 2]) == 1

    ### PL can't deliver Gde Bretagne
    assert not any(
        stop.name == "Toulouse/Grande-Bretagne(Casselardit)"
        for (_, v), tour in sol.tours.items()
        if v == 0
        for stop in tour
    )

    ### Revel must be delivered on tuesdays with an exclusive tour
    assert len(sol.tours[1, 2]) == 1
    assert sol.tours[1, 2][0].name == "Revel"

    assert not any(
        stop.name == "Revel"
        for (d, v), tour in sol.tours.items()
        if (d, v) != (1, 2)
        for stop in tour
    )

    ### Les arènes must be delivered on Friday along with the pickup at Leclerc Blagnac
    assert any(
        stop.name == "Toulouse/Negogousses(Les Arenes)" for stop in sol.tours[4, 3]
    )
    assert any(stop.name == "Leclerc Blagnac" for stop in sol.tours[4, 3])

    # Don't visit outside of allowed/required days
    for (d, v), tour in sol.tours.items():
        for stop in tour:
            if stop.type == StopType.Livraison:
                assert d in pb.centres[stop.index].allowed_days
            else:
                assert d in pb.pdrs[stop.index - pb.n_centres].required_days

    # Check that demands are met
    for c in range(pb.n_centres):
        if pb.centres[c].delivery_week in (DeliveryWeek.ANY, pb.params.week):
            deliv = [
                stop.delivery
                for tour in sol.tours.values()
                for stop in tour
                if stop.index == c
            ]
            for t in ProductType:
                assert sum(d[t.value] for d in deliv) >= pb.centres[c].demands[t]

    print("No constraint violations detected - the solution is valid.")


def test_file():
    file = "solutions/test.json"

    check_solution_file(
        "data/centres_variations/centres_keep.xlsx",
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/distance_matrix.xlsx",
        "data/traffic_duration_matrix.xlsx",
        "data/no_traffic_duration_matrix.xlsx",
        "data/params.json",
        file,
    )


if __name__ == "__main__":
    file = sys.argv[1]
    week = int(sys.argv[2])

    check_solution_file(
        "data/centres_keep.xlsx",
        "data/points_de_ramasse.xlsx",
        "data/vehicules.xlsx",
        "data/distance_matrix.xlsx",
        "data/traffic_duration_matrix.xlsx",
        "data/no_traffic_duration_matrix.xlsx",
        "data/params.json",
        file,
    )
