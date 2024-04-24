import argparse
from collections import defaultdict
from utils.problem import Problem, Solution, ProductType, DeliveryWeek, StopType, Stop
from numpy import int64


def get_trips(tour: list[Stop]):
    tour_indices = [stop.index for stop in tour]

    retours = [i for i in range(len(tour_indices)) if tour_indices[i] == 0]
    n_trips = len(retours) + 1

    split_indices = [-1] + retours + [len(tour_indices)]
    trips = [tour[split_indices[i] + 1 : split_indices[i + 1]] for i in range(n_trips)]
    return trips


def check_durations(pb: Problem, sol: Solution):
    actual_dist = 0
    for (d, v), tour in sol.tours.items():
        tour_indices = [stop.index for stop in tour]
        for a, b in zip([0] + tour_indices, tour_indices + [0]):
            new_dist = pb.distance_matrix.iloc[a, b]
            if type(new_dist) != int64:
                continue
            actual_dist += new_dist
    assert actual_dist == int(sol.total_distance)

    def duration(a, b, tour_duration):
        if tour_duration <= 90:
            arc_duration = pb.duration_matrix.iloc[a, b] / 60  # type: ignore
        else:
            arc_duration = pb.no_traffic_duration_matrix.iloc[a, b] / 60  # type: ignore
        return arc_duration

    for (d, v), tour in sol.tours.items():
        tour_duration = 0

        trips = get_trips(tour)
        assert len(trips) <= pb.params.max_trips

        for i, trip in enumerate(trips):
            trip_duration = 0
            a = 0
            for stop in trip:
                b = stop.index
                arc_duration = duration(a, b, tour_duration)
                tour_duration += arc_duration
                trip_duration += arc_duration
                if stop.type == StopType.Livraison:
                    tour_duration += pb.params.wait_at_centres
                    trip_duration += pb.params.wait_at_centres
                elif stop.type == StopType.Ramasse:
                    tour_duration += pb.params.wait_at_pdrs
                    trip_duration += pb.params.wait_at_pdrs
                a = b
            # Trip back to depot
            arc_duration = duration(a, 0, tour_duration)
            tour_duration += arc_duration

            if i == 0 and any(stop.type == StopType.Ramasse for stop in trip):
                assert trip_duration <= pb.params.max_tour_duration_with_pickup

                if trip_duration > 0.8 * pb.params.max_tour_duration_with_pickup:
                    print(
                        f"[TEST] Warning : trip {d, v} is ~{trip_duration//60:.0f}h{trip_duration%60:.0f}min long :",
                        [stop.name for stop in trip],
                    )

            if trip_duration >= 120 * 60:
                if any(stop.palettes[2] > 0 for stop in trip) and any(
                    stop.palettes[1] > 0 for stop in trip
                ):
                    print(
                        "[TEST] Warning : long trip with F and S on palettes :",
                        [stop.name for stop in trip],
                    )

            if i < len(trips) - 1:
                tour_duration += pb.params.wait_between_trips

        if tour_duration > 0.85 * pb.params.max_tour_duration:
            print(
                f"[TEST] Warning : tour {d, v} is ~{tour_duration//60:.0f}h{tour_duration%60:.0f}min long :",
                [stop.name for stop in tour],
            )

        if any(stop.name in ["BAGNERES DE LUCHON", "MONTREJEAU"] for stop in tour):
            assert tour_duration <= pb.params.max_tour_duration + 90 * 60
        else:
            assert tour_duration <= pb.params.max_tour_duration


def check_load_constraints(pb: Problem, sol: Solution):
    for (_, v), tour in sol.tours.items():
        for stop in tour:
            if stop.type == StopType.Ramasse:
                continue
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
            assert stop.delivery[0] <= pb.params.max_palette_capacity * stop.palettes[0]

            # Number of A palettes must be int
            assert stop.palettes[0] % 1 == 0

            full, rest = divmod(stop.palettes[1], 1)
            halves = rest * 2
            assert (
                2 * stop.delivery[1]
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


def check_capacity_constraints(pb: Problem, sol: Solution):
    for (d, v), tour in sol.tours.items():
        trips = get_trips(tour)
        for trip in trips:
            # Capacity limit for deliveries
            assert (
                sum(
                    sum(stop.delivery)
                    for stop in trip
                    if stop.type == StopType.Livraison
                )
                <= pb.vehicles[v].capacity
            )

            # Capacity limit for pickups
            assert (
                sum(
                    pb.pdrs[stop.index - pb.n_centres].weight
                    for stop in trip
                    if stop.type == StopType.Ramasse
                )
                <= pb.vehicles[v].capacity
            )

            # Size limit for deliveries
            assert (
                sum(
                    sum(stop.palettes)
                    for stop in trip
                    if stop.type == StopType.Livraison
                )
                <= pb.vehicles[v].size
            )

            # Size limit for pickups (each pickup takes 2 palettes)
            assert (
                sum(stop.type == StopType.Ramasse for stop in trip) * 2
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


def check_specific_requirements(pb: Problem, sol: Solution):
    ### Carrefour centrale must be visited by these vehicles on these days
    assert (2, 0) in sol.tours
    assert (2, 6) in sol.tours
    assert (4, 4) in sol.tours

    assert "CARREFOUR EN JACCA" in [node.name for node in sol.tours[2, 0]]
    assert "CARREFOUR LA MENUDE" in [node.name for node in sol.tours[2, 6]]
    assert "CARREFOUR LA MENUDE" in [node.name for node in sol.tours[4, 4]]

    ### No other pickups
    assert sum(stop.type == StopType.Ramasse for stop in sol.tours[2, 0]) == 1
    assert sum(stop.type == StopType.Ramasse for stop in sol.tours[2, 6]) == 1
    assert sum(stop.type == StopType.Ramasse for stop in sol.tours[4, 4]) == 1

    ### PL can't deliver Gde Bretagne
    assert not any(
        stop.name == "TOULOUSE GRANDE BRETAGNE"
        for (_, v), tour in sol.tours.items()
        if v == 0
        for stop in tour
    )

    ### PL can't deliver Escalquens
    assert not any(
        stop.name == "ESCALQUENS"
        for (_, v), tour in sol.tours.items()
        if v == 0 or v == 1
        for stop in tour
    )

    ### Bessières can't be delivered with camions frigos (or PL)
    assert not any(
        stop.name == "BESSIERES"
        for (_, v), tour in sol.tours.items()
        if ProductType.F in pb.vehicles[v].can_carry
        for stop in tour
    )

    ### St Orens has to be picked up with CF 2
    assert not any(
        stop.name == "LECLERC ST ORENS"
        for (_, v), tour in sol.tours.items()
        if v != 2
        for stop in tour
    )

    ### Fenouillet not too early
    for (d, v), tour in sol.tours.items():
        assert (
            tour[0].name != "FENOUILLET"
            or sol.tour_durations[d, v] <= (pb.params.max_tour_duration - 45) * 60
        )

    ### Fronton not too late
    assert not any(
        tour[i].name == "FRONTON"
        for tour in sol.tours.values()
        for i in range(1, len(tour))
    )

    ### Livraisons de Ramasse are done with no load
    for d, _, _, c in pb.livraisons_de_ramasses:
        done = False
        for v in pb.allowed_vehicles:
            if (d, v) not in sol.tours:
                continue

            current_trip = 0
            for stop in sol.tours[d, v]:
                if stop.index == 0:
                    current_trip += 1
                    continue

                if stop.index == c:
                    done = True
                    done &= all(deliv == 0 for deliv in stop.delivery[:2])
                    if d not in pb.centres[c].allowed_days:
                        done &= stop.delivery[2] == 0
                    if done:
                        break
            if done:
                break

        assert done


def check_time_window_constraints(pb: Problem, sol: Solution):
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

    n_visits = defaultdict(int)
    # Don't visit outside of allowed/required days
    for (d, v), tour in sol.tours.items():
        trip = 0
        trip_used = False
        for stop in tour:
            if stop.index == 0:
                # Make sure trips are done in order
                assert trip_used
                trip += 1
                trip_used = False
                continue

            trip_used = True

            if stop.type == StopType.Livraison:
                assert (
                    d in pb.centres[stop.index].allowed_days
                    or (d, v, trip, stop.index) in pb.livraisons_de_ramasses
                )
                if (d, v, trip, stop.index) not in pb.livraisons_de_ramasses:
                    n_visits[stop.index] += 1
            else:
                assert d in pb.pdrs[stop.index - pb.n_centres].required_days

                # No pickups outside of first trip
                assert trip == 0

    for c, n_vis in n_visits.items():
        if n_vis > 3:
            print(f"[TEST] Warning : {pb.centres[c].name} is visited {n_vis} times")


def check_demand_constraints(pb: Problem, sol: Solution):
    for c in range(1, pb.n_centres):
        if pb.centres[c].delivery_week in (DeliveryWeek.ANY, sol.week):
            deliv = [
                stop.delivery
                for tour in sol.tours.values()
                for stop in tour
                if stop.index == c
            ]
            for t in ProductType:
                assert sum(d[t.value] for d in deliv) >= sum(
                    dem.weight for dem in pb.demands[c] if dem.product_type == t
                ), f"{c} {t}"


def check_solution(pb: Problem, sol: Solution):
    check_durations(pb, sol)
    check_load_constraints(pb, sol)
    check_capacity_constraints(pb, sol)
    check_specific_requirements(pb, sol)
    check_time_window_constraints(pb, sol)
    check_demand_constraints(pb, sol)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("instance")
    parser.add_argument("solution")

    args = parser.parse_args()
    pb = Problem.from_json(args.instance)
    sol = Solution.from_json(args.solution)
    check_solution(pb, sol)
    print("[TEST] No constraint violations detected - the solution is valid.")
