import argparse
from utils.problem import Problem, Solution, ProductType, DeliveryWeek, StopType
from models.cp_solver import solve_vrp
from numpy import int64
import pytest
from os import remove
import subprocess


def check_solution_file(
    problem_file,
    solution_file,
    warning_duration_threshold=200,  # Show a warning if estimated duration for a tour is higher than this
):
    """Checks whether a solution file satisfies all constraints specified in a problem file"""
    sol = Solution.read_from_json(solution_file)
    pb = Problem.from_json(problem_file)

    # Check the durations
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
                estimated_duration += pb.no_traffic_duration_matrix.iloc[a, b] / 60  # type: ignore

        estimated_duration += sum(
            pb.params.wait_at_centres * (stop.type == StopType.Livraison)
            + pb.params.wait_at_pdrs * (stop.type == StopType.Ramasse)
            for stop in tour
        )

        if estimated_duration > warning_duration_threshold:
            print(
                f"[TEST] Warning : tour {d, v} is ~{estimated_duration//60:.0f}h{estimated_duration%60:.0f}min long :",
                [stop.name for stop in tour],
            )

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

    print("[TEST] No constraint violations detected - the solution is valid.")


test_instances = [
    "problems/centres_fixe_1.json",
    "problems/centres_fixe_2.json",
    "problems/centres_flex_1.json",
    "problems/centres_flex_2.json",
    "problems/centres_any_1.json",
    "problems/centres_any_2.json",
    "problems/centres_1.json",
    "problems/centres_2.json",
    "problems/centres_adjust_1.json",
    "problems/centres_adjust_2.json",
]


@pytest.mark.parametrize("instance", test_instances)
def test_ortools(instance):
    tmp_file = "solutions/tmp.json"
    time_limit = 10
    problem = Problem.from_json(instance)
    solve_vrp(problem, outfile=tmp_file, time_limit=time_limit)
    check_solution_file(instance, tmp_file)
    remove(tmp_file)


@pytest.mark.parametrize("instance", test_instances)
def test_hexaly(instance):
    tmp_file = "solutions/tmp.json"
    time_limit = 10
    subprocess.run(
        [
            "localsolver",
            "models/ls_solver.lsp",
            instance,
            "nil",
            tmp_file,
            str(time_limit),
        ]
    )
    check_solution_file(instance, tmp_file)
    remove(tmp_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("instance")
    parser.add_argument("solution")

    args = parser.parse_args()
    check_solution_file(args.instance, args.solution)
