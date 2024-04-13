import argparse
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

            if i == 0 and any(stop.type == StopType.Ramasse for stop in trip):
                assert trip_duration <= pb.params.max_tour_duration_with_pickup

                if trip_duration > 0.8 * pb.params.max_tour_duration_with_pickup:
                    print(
                        f"[TEST] Warning : trip {d, v} is ~{tour_duration//60:.0f}h{tour_duration%60:.0f}min long :",
                        [stop.name for stop in trip],
                    )

            if i < len(trips) - 1:
                tour_duration += pb.params.wait_between_trips

        assert tour_duration <= pb.params.max_tour_duration

        if tour_duration > 0.8 * pb.params.max_tour_duration:
            print(
                f"[TEST] Warning : tour {d, v} is ~{tour_duration//60:.0f}h{tour_duration%60:.0f}min long :",
                [stop.name for stop in tour],
            )


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

            assert (
                2 * stop.delivery[1]
                <= 2 * stop.palettes[1] * pb.params.demi_palette_capacity
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

    ### Livraisons de Ramasse are done with no load
    for d, v, trip, c in pb.livraisons_de_ramasses:
        done = False
        current_trip = 0
        for stop in sol.tours[d, v]:
            if stop.index == 0:
                current_trip += 1
                continue

            if stop.index == c and current_trip == trip:
                done = True
                assert all(d == 0 for d in stop.delivery)
                break
        if not done:
            print(d, v, trip, c)
            print(sol.tours[d, v])
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

    # Don't visit outside of allowed/required days
    for (d, v), tour in sol.tours.items():
        trip = 0
        for stop in tour:
            if stop.index == 0:
                trip += 1
                continue

            if stop.type == StopType.Livraison:
                assert (
                    d in pb.centres[stop.index].allowed_days
                    or (d, v, trip, stop.index) in pb.livraisons_de_ramasses
                )
            else:
                assert d in pb.pdrs[stop.index - pb.n_centres].required_days

                # No pickups outside of first trip
                assert trip == 0


def check_demand_constraints(pb: Problem, sol: Solution):
    demands_met = []

    # Check that demands are met for each scenario
    for scenario in range(len(pb.centres[0].demands)):
        met = True
        for c in range(pb.n_centres):
            if pb.centres[c].delivery_week in (DeliveryWeek.ANY, pb.params.week):
                deliv = [
                    stop.delivery
                    for tour in sol.tours.values()
                    for stop in tour
                    if stop.index == c
                ]
                for t in ProductType:
                    if (
                        sum(d[t.value] for d in deliv)
                        < pb.centres[c].demands[scenario][t]
                    ):
                        met = False
        demands_met.append(met)

    assert any(demands_met)
    print(f"[TEST] Demands met in {sum(demands_met)}/{len(demands_met)} scenarios")


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
