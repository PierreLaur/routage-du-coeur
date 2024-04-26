import json
from ortools.sat.python import cp_model
from math import inf, ceil
from utils.datatypes import (
    DeliveryWeek,
    ProductType,
    Stop,
    StopType,
)
from utils.problem import Problem, Solution
from typing import Any
import numpy as np
from copy import deepcopy


class LogPrinter(cp_model.CpSolverSolutionCallback):
    """Print the current objective values and the optimality gap as the search progresses"""

    def __init__(self, vars):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.best_solution = inf
        self.vars = vars

    def on_solution_callback(self):
        obj = self.objective_value

        if obj < self.best_solution:
            self.best_solution = obj
            if self.best_objective_bound == 0:
                return
            gap = (
                100
                * (self.best_solution - self.best_objective_bound)
                / self.best_objective_bound
            )
            n_used = self.value(sum(self.vars.used.values()))
            print(
                f"[CP-SAT] {f'[{self.wall_time:.2f}s]':8}  cost={obj:<8.2f}€ "
                f"| n_vehicles={n_used:<2} "
                f"| dist={self.value(self.vars.total_distance)/1000:.1f}km   "
                f"gap = {gap:.2f}% [{self.best_objective_bound:.2f}, {obj:.2f}]"
            )


class VarContainer:
    fixed_costs: Any
    variable_costs: Any
    total_costs: Any
    total_distance: Any

    def __init__(self) -> None:
        self.arcs: dict[tuple[int, int, int], np.ndarray] = {}
        self.visits: dict[tuple[int, int, int, int], cp_model.IntVar] = {}
        self.delivers: dict[tuple[int, int, int, int, int], cp_model.IntVar] = {}
        self.route_dist: dict[tuple[int, int], Any] = {}
        self.tour_duration: dict[tuple[int, int], Any] = {}
        self.trip_duration: dict[tuple[int, int, int], Any] = {}
        self.used: dict[int, cp_model.IntVar] = {}


def add_hint(
    model: cp_model.CpModel,
    pb: Problem,
    vars: VarContainer,
    init_sol: Solution,
) -> None:
    """
    Adds an initial solution as a hint to the model.
    """

    for d in range(pb.n_days):
        for v in range(pb.m):
            if not pb.vehicles[v].allowed:
                continue

            if (d, v) not in init_sol.tours:
                for trip in range(pb.params.max_trips):
                    model.add_hint(vars.visits[d, v, trip, 0], 0)
                continue

            trip = 0
            to_assign = set(
                (trip, start, end)
                for trip in range(pb.params.max_trips)
                for start in range(pb.n_centres + pb.n_pdr)
                for end in range(pb.n_centres + pb.n_pdr)
                if start != end
            )

            if (d, v, 0, 0) in vars.visits:
                model.add_hint(vars.visits[d, v, 0, 0], 1)

            a = 0
            trip = 0

            tour = init_sol.tours[d, v]
            for stop in tour:
                b = stop.index
                model.add_hint(vars.arcs[d, v, trip][a, b], 1)
                to_assign.remove((trip, a, b))
                a = b
                if b == 0:
                    trip += 1
                    continue

            model.add_hint(vars.arcs[d, v, trip][a, 0], 1)
            to_assign.remove((trip, a, 0))

            # for trip, start, end in to_assign:
            #     model.add_hint(vars.arcs[d, v, trip][start, end], 0)


def set_current_tours(
    pb: Problem,
    model: cp_model.CpModel,
    vars: VarContainer,
):
    """Sets the current tours as a hint to the solver"""
    # current_tours_file = (
    #     f"data/current/tours_tournees_actuelles_w{pb.params.week.value}.json"
    # )
    current_tours_file = "data/current/2602.json"
    with open(current_tours_file, "r") as f:
        tours = json.load(f)
        for d in range(pb.n_days):
            for v in range(pb.m):
                key = f"({d}, {v})"
                if key in tours:
                    trip = 0
                    n_stops = 0

                    for stop in tours[key][1:-1]:
                        if stop == 0:
                            model.Add(
                                sum(
                                    vars.visits[d, v, trip, s]
                                    for s in range(1, pb.n_centres + pb.n_pdr)
                                    if (d, v, trip, s) in vars.visits
                                )
                                == n_stops
                            )
                            trip += 1
                            n_stops = 0
                            continue
                        model.add_hint(vars.visits[d, v, trip, stop], 1)
                        n_stops += 1

                    model.Add(
                        sum(
                            vars.visits[d, v, trip, s]
                            for s in range(1, pb.n_centres + pb.n_pdr)
                            if (d, v, trip, s) in vars.visits
                        )
                        == n_stops
                    )
                    for t in range(trip + 1, pb.params.max_trips):
                        model.add_hint(vars.visits[d, v, t, 0], 0)

                else:
                    for trip in range(pb.params.max_trips):
                        if (d, v, trip, 0) in vars.visits:
                            model.add_hint(vars.visits[d, v, trip, 0], 0)


def prune_arcs(n_nodes, matrix, limit=1.0):
    """Evaluates whether each arc is longer than a given ratio of the maximum distance between nodes in the matrix, and returns the list of arcs that should be pruned"""
    max_dist = matrix.max().max()
    pruned = []
    for node1 in range(n_nodes):
        for node2 in range(n_nodes):
            if matrix.iloc[node1, node2] > limit * max_dist:
                pruned.append((node1, node2))
    return pruned


def read_tours(pb: Problem, vars: VarContainer, solver):
    """Reads tours from the solver
    Returns a dictionary of tours with keys (day, vehicle)"""

    tours: dict[tuple[int, int], list[int]] = {}
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            tours[d, v] = []
            for trip in range(pb.params.max_trips):
                if not solver.value(vars.visits[d, v, trip, 0]):
                    continue

                tour = [0]
                while True:
                    a = tour[-1]
                    goes_to = set()

                    for b in range(pb.n_centres + pb.n_pdr):
                        if a == b:
                            continue
                        arc = vars.arcs[d, v, trip][a, b]

                        if solver.value(arc):
                            goes_to.add(b)

                    if len(goes_to) > 1:
                        print("Warning : ", d, v, a, "goes to ", goes_to)
                    if len(goes_to) == 0:
                        if a != 0:
                            print("Warning : ", d, v, a, "goes nowhere ")
                        break
                    goto = goes_to.pop()
                    if goto == 0:
                        break
                    tour.append(goto)

                if len(tour) == 1:
                    continue

                tours[d, v] += tour
            tours[d, v].append(0)
    return tours


def make_solution(
    pb: Problem, week: DeliveryWeek, vars: VarContainer, solver: cp_model.CpSolver
):
    tours_flat = read_tours(pb, vars, solver)
    delivers = {k: solver.value(v) for k, v in vars.delivers.items()}

    tour_duration = {
        k: solver.value(v) for k, v in vars.tour_duration.items() if k in tours_flat
    }

    tours = {}
    for (d, v), tour in tours_flat.items():
        if len(tour) <= 1:
            continue
        t = []
        trip = 0
        for node in tour[1:-1]:
            if node >= pb.n_centres:
                name = pb.pdrs[node - pb.n_centres].name
                stop_type = StopType.Ramasse
            else:
                name = pb.centres[node].name
                stop_type = (
                    StopType.Liv_Ramasse
                    if (d, v, trip, node) in pb.livraisons_de_ramasses
                    else StopType.Livraison
                )
            stop = Stop(node, name, stop_type)

            if node == 0:
                trip += 1
            elif node < pb.n_centres:
                stop.delivery = (
                    sum(
                        delivers[d, v, trip, node, dem] * pb.demands[node][dem].weight
                        for dem in range(len(pb.demands[node]))
                        if pb.demands[node][dem].product_type == ProductType.A
                    ),
                    sum(
                        delivers[d, v, trip, node, dem] * pb.demands[node][dem].weight
                        for dem in range(len(pb.demands[node]))
                        if pb.demands[node][dem].product_type == ProductType.F
                    ),
                    sum(
                        delivers[d, v, trip, node, dem] * pb.demands[node][dem].weight
                        for dem in range(len(pb.demands[node]))
                        if pb.demands[node][dem].product_type == ProductType.S
                    ),
                )

                stop.palettes = (
                    sum(
                        delivers[d, v, trip, node, dem] * pb.demands[node][dem].palettes
                        for dem in range(len(pb.demands[node]))
                        if pb.demands[node][dem].product_type == ProductType.A
                    ),
                    sum(
                        delivers[d, v, trip, node, dem] * pb.demands[node][dem].palettes
                        for dem in range(len(pb.demands[node]))
                        if pb.demands[node][dem].product_type == ProductType.F
                    ),
                    sum(
                        delivers[d, v, trip, node, dem] * pb.demands[node][dem].palettes
                        for dem in range(len(pb.demands[node]))
                        if pb.demands[node][dem].product_type == ProductType.S
                    ),
                )

                stop.norvegiennes = sum(
                    delivers[d, v, trip, node, dem] * pb.demands[node][dem].norvegiennes
                    for dem in range(len(pb.demands[node]))
                )

            t.append(stop)
        tours[d, v] = t

    sol = Solution(
        week,
        round(solver.value(vars.total_costs), 3),  # type: ignore
        round(solver.value(vars.fixed_costs), 3),  # type: ignore
        round(solver.value(vars.variable_costs), 3),  # type: ignore
        solver.value(vars.total_distance),
        [
            bool(solver.value(vars.used[v])) if v in pb.allowed_vehicles else False
            for v in range(pb.m)
        ],
        tours,
        {k: solver.value(v) for k, v in tour_duration.items() if k in tours_flat},
    )
    sol.adjust_durations(pb)
    return sol


def create_tour_variables(
    pb: Problem,
    week: DeliveryWeek,
    model: cp_model.CpModel,
    vars: VarContainer,
    use_all_vehicles=False,
):
    """Creates arcs and visits variables and loads them into the VarContainer and the model"""
    n_nodes = pb.n_nodes
    used_nodes = pb.used_nodes(week)
    allowed_vehicles = pb.allowed_vehicles

    # Optionally prune some long arcs to make solving easier
    pruned = prune_arcs(n_nodes, pb.distance_matrix, limit=1)
    if pruned:
        print(f"    Pruned arcs : {100*len(pruned)/(n_nodes)**2:.2f}%")

    vars.used = {v: model.new_bool_var("") for v in allowed_vehicles}
    if use_all_vehicles:
        model.add_bool_and(vars.used.values())

    ########## Arcs & Visit variables ###########
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            for trip in range(pb.params.max_trips):
                vars.arcs[d, v, trip] = np.array(
                    [
                        [model.new_constant(0) for _ in range(n_nodes)]
                        for _ in range(n_nodes)
                    ]
                )
                for node in used_nodes:
                    # pickup only on the first trip
                    if trip > 0 and node >= pb.n_centres:
                        vars.visits[d, v, trip, node] = model.new_constant(0)
                    else:
                        vars.visits[d, v, trip, node] = model.new_bool_var(
                            "visits_%i_%i_%i_%i" % (d, v, trip, node)
                        )

                    # Add an arc from the node to itself (necessary for Circuit constraints)
                    # If used, the node is not visited
                    vars.arcs[d, v, trip][node, node] = vars.visits[
                        d, v, trip, node
                    ].negated()

                    for node2 in used_nodes:
                        if node == node2:
                            continue
                        if (node, node2) in pruned:
                            continue
                        # No arc from pdr to centre (except depot) - we do deliveries first then pickups
                        elif (
                            node >= pb.n_centres and 0 < node2 and node2 < pb.n_centres
                        ):
                            continue
                        else:
                            vars.arcs[d, v, trip][node, node2] = model.new_bool_var("")

                # If you don't visit the depot, visit nothing
                model.add_bool_and(
                    vars.visits[d, v, trip, node].negated()
                    for node in pb.delivered_centres(week)
                    + list(range(pb.n_centres, pb.n_centres + pb.n_pdr))
                ).only_enforce_if(vars.visits[d, v, trip, 0].negated())

                # Trips should be done in order
                for prev_trip in range(trip):
                    model.add_implication(
                        vars.visits[d, v, trip, 0], vars.visits[d, v, prev_trip, 0]
                    )

                model.add_circuit(
                    [
                        (node1, node2, vars.arcs[d, v, trip][node1, node2])
                        for node1 in used_nodes
                        for node2 in used_nodes
                    ]
                )


def create_load_variables(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    """Creates the deliveries, palettes and norvegiennes variables and loads them into the VarContainer and the model"""

    ########## Load variables & constraints ###########
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            for trip in range(pb.params.max_trips):
                for c in pb.delivered_centres(week):
                    for i, demand in enumerate(pb.demands[c]):
                        can_deliver = True
                        if demand.product_type not in pb.vehicles[v].can_carry:
                            can_deliver = False
                        if (
                            demand.product_type == ProductType.S
                            and demand.palettes > 0
                            and not pb.vehicles[v].allows_isotherm_cover
                        ):
                            can_deliver = False

                        if can_deliver:
                            vars.delivers[d, v, trip, c, i] = model.new_bool_var("")
                        else:
                            vars.delivers[d, v, trip, c, i] = model.new_constant(0)


def create_duration_variables(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    delivered_centres = pb.delivered_centres(week)
    pdr_nodes = pb.pdr_nodes
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            for trip in range(pb.params.max_trips):
                vars.trip_duration[d, v, trip] = np.sum(
                    pb.duration_matrix.values * vars.arcs[d, v, trip]
                ) + 60 * (
                    pb.params.wait_at_centres
                    * sum(vars.visits[d, v, trip, c] for c in delivered_centres)
                    + pb.params.wait_at_pdrs
                    * sum(vars.visits[d, v, trip, p] for p in pdr_nodes)
                )  # type: ignore

            vars.route_dist[d, v] = sum(
                np.sum(pb.distance_matrix.values * vars.arcs[d, v, trip])
                for trip in range(pb.params.max_trips)
            )
            vars.tour_duration[d, v] = sum(
                vars.trip_duration[d, v, trip] for trip in range(pb.params.max_trips)
            ) + sum(
                vars.visits[d, v, trip, 0] * pb.params.wait_between_trips * 60
                for trip in range(1, pb.params.max_trips)  # type: ignore
            )


def add_load_constraints(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    """Adds the deliveries, palettes and norvegiennes constraints to the model"""
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            for trip in range(pb.params.max_trips):
                for c in pb.delivered_centres(week):
                    # Not visited -> delivery is 0
                    model.add_bool_and(
                        vars.delivers[d, v, trip, c, dem].negated()
                        for dem in range(len(pb.demands[c]))
                    ).only_enforce_if(vars.visits[d, v, trip, c].negated())


def add_time_windows(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    allowed_vehicles = pb.allowed_vehicles

    for d in range(pb.n_days):
        for pdr in range(pb.n_pdr):
            # Visit each exactly once on each required days with an appropriate vehicle
            if d in pb.pdrs[pdr].required_days:
                model.AddExactlyOne(
                    vars.visits[d, v, trip, pdr + pb.n_centres]
                    for v in allowed_vehicles
                    for trip in range(pb.params.max_trips)
                    if pb.pdrs[pdr].product_type in pb.vehicles[v].can_carry
                )

                # Disallow any ramasses with vehicles that cannot carry that type of product
                model.add_bool_and(
                    vars.visits[d, v, trip, pdr + pb.n_centres].negated()
                    for v in allowed_vehicles
                    for trip in range(pb.params.max_trips)
                    if pb.pdrs[pdr].product_type not in pb.vehicles[v].can_carry
                )
            else:
                # Do not visit them if not required (this seems to improve pruning)
                model.add_bool_and(
                    [
                        vars.visits[d, v, trip, pdr + pb.n_centres].negated()
                        for v in allowed_vehicles
                        for trip in range(pb.params.max_trips)
                    ]
                )

    # Only deliver centres on allowed days
    for c in pb.delivered_centres(week):
        for d in range(pb.n_days):
            if d not in pb.centres[c].allowed_days:
                model.add_bool_and(
                    [
                        vars.visits[d, v, trip, c].negated()
                        for v in allowed_vehicles
                        for trip in range(pb.params.max_trips)
                        if (d, v, trip, c) not in pb.livraisons_de_ramasses
                    ]
                )


def add_cover_constraints(
    pb: Problem, week: DeliveryWeek, model, vars: VarContainer, loose=False
):
    fourgon_capacity = 1200
    for c in pb.delivered_centres(week):
        # Redundant Cover constraint
        sum_demands = sum(dem.weight for dem in pb.demands[c])

        if sum_demands > fourgon_capacity:
            min_visits = 2
        else:
            min_visits = 1

        # # Heuristic pruning rule (very effective)
        max_visits = ceil(sum_demands / pb.params.max_palette_capacity)

        # more than 3 visits is unreasonable in practice
        max_visits = min(max_visits, 3)

        if loose:
            max_visits = 3

        model.add_linear_constraint(
            sum(
                vars.visits[d, v, trip, c]
                for d in range(pb.n_days)
                for v in pb.allowed_vehicles
                for trip in range(pb.params.max_trips)
                if d in pb.centres[c].allowed_days
                and (d, v, trip, c) not in pb.livraisons_de_ramasses
            ),
            min_visits,
            max_visits,
        )


def add_demand_constraints(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    # Every customer must be served his demand
    for c in pb.delivered_centres(week):
        for dem in range(len(pb.demands[c])):
            model.add_exactly_one(
                vars.delivers[d, v, trip, c, dem]
                for d in range(pb.n_days)
                for v in pb.allowed_vehicles
                for trip in range(pb.params.max_trips)
            )


def add_capacity_constraints(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    delivered_centres = pb.delivered_centres(week)
    allowed_vehicles = pb.allowed_vehicles

    # Capacity constraints
    for d in range(pb.n_days):
        model.add(
            sum(
                vars.delivers[d, v, trip, c, dem] * pb.demands[c][dem].norvegiennes
                for v in allowed_vehicles
                for trip in range(pb.params.max_trips)
                for c in delivered_centres
                for dem in range(len(pb.demands[c]))
            )
            <= pb.params.n_norvegiennes
        )

        for v in allowed_vehicles:
            for trip in range(pb.params.max_trips):
                # Delivery
                model.add(
                    sum(
                        vars.delivers[d, v, trip, c, dem] * pb.demands[c][dem].weight
                        for c in delivered_centres
                        for dem in range(len(pb.demands[c]))
                    )
                    <= pb.vehicles[v].capacity
                )

                # Pickup
                model.add(
                    sum(
                        pb.pdrs[p].weight * vars.visits[d, v, trip, p + pb.n_centres]
                        for p in range(pb.n_pdr)
                    )
                    <= pb.vehicles[v].capacity
                )

                # Size limits
                model.add(
                    sum(
                        vars.delivers[d, v, trip, c, dem]
                        * int(2 * pb.demands[c][dem].palettes)
                        for c in delivered_centres
                        for dem in range(len(pb.demands[c]))
                    )
                    <= 2 * pb.vehicles[v].size
                )

                model.add(
                    sum(
                        pb.pdrs[p - pb.n_centres].palettes * vars.visits[d, v, trip, p]
                        for p in pb.pdr_nodes
                    )
                    <= pb.vehicles[v].size
                )


def add_duration_constraints(
    pb: Problem,
    week: DeliveryWeek,
    model: cp_model.CpModel,
    vars: VarContainer,
    loose=False,
):
    delivered_centres = pb.delivered_centres(week)
    pdr_nodes = pb.pdr_nodes

    # Duration constraints
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            if loose:
                if (d, v) in [(1, 7), (1, 0)]:
                    continue
            pickups = model.new_bool_var("")
            model.add_bool_and(
                vars.visits[d, v, 0, p].negated() for p in pdr_nodes
            ).only_enforce_if(pickups.negated())
            model.add(
                vars.trip_duration[d, v, 0]
                <= pb.params.max_tour_duration_with_pickup * 60
            ).only_enforce_if(pickups)

            model.add(vars.tour_duration[d, v] <= pb.params.max_tour_duration * 60)

            model.add(
                sum(
                    vars.visits[d, v, trip, node]
                    for node in delivered_centres + list(pdr_nodes)
                    for trip in range(pb.params.max_trips)
                )
                <= pb.params.max_stops
            )


def add_specific_requirements(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    """
    Adds a few "hard-coded" constraints
    These are not written in the input files but have been specified by the client
    """

    # The PL must visit Carrefour En Jacca on Wednesday
    # He must only do one pickup
    carrefour_en_jacca_index = pb.get_index_by_name("CARREFOUR EN JACCA")
    model.add(vars.visits[2, 0, 0, carrefour_en_jacca_index] == 1)
    model.add_bool_and(
        vars.visits[2, 0, 0, p].negated()
        for p in pb.pdr_nodes
        if p != carrefour_en_jacca_index
    )

    # A camion Frigo must visit Carrefour la menude on Wednesday and Friday
    # To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    carrefour_menude_index = pb.get_index_by_name("CARREFOUR LA MENUDE")
    model.add(vars.visits[2, 6, 0, carrefour_menude_index] == 1)
    model.add_bool_and(
        vars.visits[2, 6, 0, p].negated()
        for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
        if p != carrefour_menude_index
    )

    model.add(vars.visits[4, 4, 0, carrefour_menude_index] == 1)
    model.add_bool_and(
        vars.visits[4, 4, 0, p].negated()
        for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
        if p != carrefour_menude_index
    )

    # The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = pb.get_index_by_name("TOULOUSE GRANDE BRETAGNE")
    model.add_bool_and(
        vars.visits[d, 0, trip, index_gde_bretagne].negated()
        for d in range(pb.n_days)
        for trip in range(pb.params.max_trips)
    )

    # both PLs cannot deliver Escalquens
    index_escalquens = pb.get_index_by_name("ESCALQUENS")
    model.add_bool_and(
        vars.visits[d, v, trip, index_escalquens].negated()
        for d in range(pb.n_days)
        for trip in range(pb.params.max_trips)
        for v in (0, 1)
    )

    # Bessières cannot be delivered with camions frigos or PL
    index_bessieres = pb.get_index_by_name("BESSIERES")
    model.add_bool_and(
        vars.visits[d, v, trip, index_bessieres].negated()
        for d in range(pb.n_days)
        for v in pb.allowed_vehicles
        if (ProductType.F in pb.vehicles[v].can_carry or v <= 1)
        for trip in range(pb.params.max_trips)
        if (d, v, trip, index_bessieres) in vars.visits
    )

    # St Orens must be picked up with a specific CF (the smallest one, for height limit reasons)
    index_st_orens = pb.get_index_by_name("LECLERC ST ORENS")
    model.add_bool_and(
        vars.visits[d, v, 0, index_st_orens].negated()
        for d in range(pb.n_days)
        for v in pb.allowed_vehicles
        if v != 2
    )

    # Fenouillet not too early
    index_fenouillet = pb.get_index_by_name("FENOUILLET")
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            for trip in range(pb.params.max_trips):
                model.add(
                    vars.tour_duration[d, v] <= (pb.params.max_tour_duration - 45) * 60
                ).only_enforce_if(vars.visits[d, v, trip, index_fenouillet])

    # Fronton in first only (they distribute food in the morning)
    index_fronton = pb.get_index_by_name("FRONTON")
    model.add_bool_and(
        vars.arcs[d, v, trip][c, index_fronton].negated()
        for d in range(pb.n_days)
        for v in pb.allowed_vehicles
        for c in pb.delivered_centres(week)
        if c != index_fronton and c != 0
        for trip in range(pb.params.max_trips)
    )
    model.add_bool_and(
        vars.arcs[d, v, trip][0, index_fronton].negated()
        for d in range(pb.n_days)
        for v in pb.allowed_vehicles
        for trip in range(1, pb.params.max_trips)
    )

    # "Livraisons de Ramasses" these deliver the food gathered in previous pickups. Only a single delivery can be done each time.
    for d, v, trip, c in pb.livraisons_de_ramasses:
        model.add_bool_and(
            vars.visits[d, v, trip, c2]
            if c2 == c
            else vars.visits[d, v, trip, c2].negated()
            for c2 in pb.delivered_centres(week)
            if (d, v, trip, c2) in vars.visits
        )
        model.add_bool_and(
            vars.delivers[d, v, trip, c, i].negated()
            for i in range(len(pb.demands[c]))
            if pb.demands[c][i].product_type != ProductType.S
            or d not in pb.centres[c].allowed_days
        )

    # Leclerc Rouffiac with camions frigos except on Tuesdays
    index_rouffiac = pb.get_index_by_name("LECLERC ROUFFIAC")
    model.add_bool_and(
        vars.visits[d, v, trip, index_rouffiac].negated()
        for d in range(pb.n_days)
        if d != 1
        for v in pb.allowed_vehicles
        if ProductType.F not in pb.vehicles[v].can_carry
        for trip in range(pb.params.max_trips)
        if (d, v, trip, index_rouffiac) in vars.visits
    )


def add_domsym_breaking(
    pb: Problem,
    week: DeliveryWeek,
    model,
    vars: VarContainer,
):
    # Dominance breaking
    for d, v, trip, c in vars.visits:
        if (d, v, trip, c) in pb.livraisons_de_ramasses:
            continue
        if c == 0 or c >= pb.n_centres:
            continue

        # Visited -> delivery is not 0
        model.add(
            sum(vars.delivers[d, v, trip, c, dem] for dem in range(len(pb.demands[c])))
            > 0
        ).only_enforce_if(vars.visits[d, v, trip, c])


def add_objective(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    allowed_vehicles = pb.allowed_vehicles

    vars.total_distance = sum(rd for rd in vars.route_dist.values())

    for v in allowed_vehicles:
        model.add_bool_and(
            vars.visits[d, v2, trip, c].negated()
            for (d, v2, trip, c) in vars.visits
            if v2 == v
        ).only_enforce_if(vars.used[v].negated())
        # model.add_bool_and(
        #     vars.visits[d, v, trip, 0].negated()
        #     for d in range(pb.n_days)
        #     for trip in range(pb.params.max_trips)
        #     if (d, v, trip, 0) in vars.visits
        # ).only_enforce_if(vars.used[v].negated())

    vars.variable_costs = (
        sum(
            vars.route_dist[d, v] * pb.vehicles[v].cost_per_km
            for d in range(pb.n_days)
            for v in allowed_vehicles
        )
        * 0.001
    )

    vars.fixed_costs = sum(
        vars.used[v] * pb.vehicles[v].fixed_cost
        for v in allowed_vehicles  # type: ignore
    )

    vars.total_costs = vars.variable_costs + vars.fixed_costs
    model.minimize(vars.total_costs)


def add_decision_strategies(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    pass


def add_manual_pruning(
    pb: Problem, week: DeliveryWeek, model: cp_model.CpModel, vars: VarContainer
):
    model.add(sum(vars.used.values()) <= 9)
    # model.add(sum(vars.used.values()) >= 7)
    # model.add(vars.total_distance <= 2_400_000)
    return


def solve_vrp(
    pb: Problem,
    week: DeliveryWeek,
    hint=None,
    time_limit=None,
    stop_at_first_solution=False,
    fix_hint=False,
    use_all_vehicles=False,
):
    """Creates a CP model and solves it with cp-sat.
    When done, re-optimizes the solution and writes it in JSON format to the specified file
    """

    print("[CP-SAT] Modeling... ", end="")

    model = cp_model.CpModel()
    vars = VarContainer()

    print("\r[CP-SAT] Creating variables...", end="")
    create_tour_variables(pb, week, model, vars, use_all_vehicles)
    create_load_variables(pb, week, model, vars)
    create_duration_variables(pb, week, model, vars)

    print("\r[CP-SAT] Adding constraints 0/10        ", end="")
    add_time_windows(pb, week, model, vars)
    add_cover_constraints(pb, week, model, vars, loose=False)
    add_duration_constraints(pb, week, model, vars, loose=False)

    print("\r[CP-SAT] Adding constraints 3/10        ", end="")
    add_load_constraints(pb, week, model, vars)
    add_capacity_constraints(pb, week, model, vars)
    add_demand_constraints(pb, week, model, vars)

    print("\r[CP-SAT] Adding constraints 6/10        ", end="")

    add_specific_requirements(pb, week, model, vars)
    add_domsym_breaking(pb, week, model, vars)
    add_decision_strategies(pb, week, model, vars)

    print("\r[CP-SAT] Adding objective        ", end="")
    add_objective(pb, week, model, vars)
    # add_manual_pruning(pb, week, model, vars)

    # Add input file as a hint
    if hint:
        add_hint(
            model,
            pb,
            vars,
            hint,
        )

    solver = cp_model.CpSolver()

    # set_current_tours(pb, model, vars)
    # solver.parameters.fix_variables_to_their_hinted_value = True
    if hint and fix_hint:
        solver.parameters.fix_variables_to_their_hinted_value = True

    solver.parameters.num_workers = 12
    solver.parameters.use_lns_only = True
    # solver.parameters.diversify_lns_params = True

    if time_limit:
        solver.parameters.max_time_in_seconds = time_limit
    if stop_at_first_solution:
        solver.parameters.stop_after_first_solution = True

    print("\r[CP-SAT] Solving...                    ")

    # solver.parameters.log_search_progress = True
    # status = solver.solve(model)
    status = solver.solve(model, LogPrinter(vars))

    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        print("No solution found")
        return status, None
    else:
        # sol = day_lns(pb, week, model, vars, solver)
        sol = make_solution(pb, week, vars, solver)
        return status, sol


def day_lns(
    pb: Problem,
    week: DeliveryWeek,
    model: cp_model.CpModel,
    vars: VarContainer,
    solver: cp_model.CpSolver,
):
    suboptimal_days = set(range(pb.n_days))
    obj = solver.objective_value
    for iteration in range(10):
        for d in suboptimal_days:
            print(
                "[LNS] iteration ",
                iteration,
                " obj = ",
                obj,
                "unfixing day ",
                d,
            )

            # Fix all days except this one
            new_model = model.clone()
            for d2 in range(pb.n_days):
                if d2 == d:
                    continue

                for v in pb.allowed_vehicles:
                    for trip in range(pb.params.max_trips):
                        for n1 in range(pb.n_nodes):
                            for n2 in range(pb.n_nodes):
                                if n1 == n2:
                                    continue

                                new_model.add(
                                    vars.arcs[d2, v, trip][n1, n2]
                                    == solver.value(vars.arcs[d2, v, trip][n1, n2])
                                )

                        # for c in pb.delivered_centres(week):
                        #     model.add_hint(
                        #         vars.visits[d2, v, trip, c],
                        #         solver.value(vars.visits[d2, v, trip, c]),
                        #     )

            solver = cp_model.CpSolver()

            if iteration % 2 == 0:
                solver.parameters.use_lns_only = True
            # solver.parameters.log_search_progress = True

            status = solver.solve(new_model, LogPrinter(vars))

            obj = solver.objective_value
            if status == cp_model.OPTIMAL:
                suboptimal_days.remove(d)
                print("Suboptimal days remaining : ", suboptimal_days)
                break


def solve_day_per_day(
    pb: Problem,
    week: DeliveryWeek,
):
    # Assign each centre to a single day
    for c in range(pb.n_centres):
        while len(pb.centres[c].allowed_days) > 1:
            pb.centres[c].allowed_days.pop()

    solutions: list[Solution] = []
    for d in range(pb.n_days):
        day_problem = deepcopy(pb)

        for c in range(pb.n_centres):
            if d not in day_problem.centres[c].allowed_days:
                match week:
                    case DeliveryWeek.ODD:
                        day_problem.centres[c].delivery_week = DeliveryWeek.EVEN
                    case DeliveryWeek.EVEN:
                        day_problem.centres[c].delivery_week = DeliveryWeek.ODD
        for p in range(pb.n_pdr):
            day_problem.pdrs[p].required_days = {d} & day_problem.pdrs[p].required_days

        status, sol = solve_vrp(day_problem, week, use_all_vehicles=True)
        if sol:
            solutions.append(sol)

    # Aggregate into a single solution
    solution = solutions[0]
    for d in range(1, len(solutions)):
        solution.variable_costs += solutions[d].variable_costs
        solution.total_costs += solutions[d].variable_costs
        solution.total_distance += solutions[d].total_distance
        solution.tour_durations.update(solutions[d].tour_durations)
        solution.tours.update(solutions[d].tours)

    print(solution.tour_durations)
    print(solution.tours)

    solution.adjust_durations(pb)

    return solution
