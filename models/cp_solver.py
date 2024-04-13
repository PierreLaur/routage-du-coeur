import json
from ortools.sat.python import cp_model
from math import inf, ceil
from utils.problem import Problem, Solution, DeliveryWeek, ProductType, Stop, StopType
from typing import Any
import numpy as np


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
        self.delivers: dict[
            tuple[int, int, int, int, ProductType], cp_model.IntVar
        ] = {}
        self.palettes: dict[tuple[int, int, int, int], cp_model.IntVar] = {}
        self.demi_palettes: dict[
            tuple[int, int, int, int, ProductType], cp_model.IntVar
        ] = {}
        self.norvegiennes: dict[tuple[int, int, int, int], cp_model.IntVar] = {}

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

                if (
                    b < pb.n_centres
                    and (d, v, trip, b) in vars.palettes
                    and pb.centres[b].delivery_week
                    in [DeliveryWeek.ANY, pb.params.week]
                ):
                    for product_type in ProductType:
                        if product_type in pb.vehicles[v].can_carry:
                            model.add_hint(
                                vars.delivers[d, v, trip, b, product_type],
                                stop.delivery[product_type.value],
                            )

                    if ProductType.A in pb.vehicles[v].can_carry:
                        model.add_hint(
                            vars.palettes[d, v, trip, b],
                            ceil(stop.palettes[0]),
                        )

                    if ProductType.F in pb.vehicles[v].can_carry:
                        model.add_hint(
                            vars.demi_palettes[d, v, trip, b, ProductType.F],
                            round(stop.palettes[1] * 2),
                        )

                    if ProductType.S in pb.vehicles[v].can_carry:
                        model.add_hint(
                            vars.norvegiennes[d, v, trip, b],
                            stop.norvegiennes,
                        )
                        if pb.vehicles[v].allows_isotherm_cover:
                            model.add_hint(
                                vars.demi_palettes[d, v, trip, b, ProductType.S],
                                round(stop.palettes[2] * 2),
                            )

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
    with open(
        f"data/current/tours_tournees_actuelles_w{pb.params.week.value}.json", "r"
    ) as f:
        tours = json.load(f)
        for d in range(pb.n_days):
            for v in range(pb.m):
                key = f"({d}, {v})"
                if key in tours:
                    for trip in range(pb.params.max_trips):
                        for c in range(1, pb.n_centres + pb.n_pdr):
                            if (d, v, c) not in vars.visits:
                                continue

                            if c in tours[key][1:-1]:
                                model.add_hint(vars.visits[d, v, trip, c], 1)
                            else:
                                model.add_hint(vars.visits[d, v, trip, c], 0)
                else:
                    for trip in range(pb.params.max_trips):
                        for c in range(pb.n_centres + pb.n_pdr):
                            if (d, v, c) in vars.visits:
                                model.add_hint(vars.visits[d, v, trip, c], 0)


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


def post_process(
    pb: Problem,
    tours,
    delivers,
    palettes,
    demi_palettes,
    norvegiennes,
    tour_duration,
):
    """Shaves any superfluous decisions (useless norvegiennes/palettes) from the solution
    Also adds the time required to come back to the depot to the tour duration"""

    for (d, v), tour in tours.items():
        # Remove any useless norvegiennes and demi palettes for s products
        for c in tour:
            if (d, v, c, ProductType.S) not in demi_palettes:
                continue

            new_n_norvegiennes = max(
                0,
                ceil(
                    (
                        delivers[d, v, c, ProductType.S]
                        - demi_palettes[d, v, c, ProductType.S]
                        * pb.params.demi_palette_capacity
                    )
                    / pb.params.norvegienne_capacity
                ),
            )
            if norvegiennes[d, v, c] != new_n_norvegiennes:
                norvegiennes[d, v, c] = new_n_norvegiennes

            new_n_demi_palettes_s = max(
                0,
                ceil(
                    (
                        delivers[d, v, c, ProductType.S]
                        - norvegiennes[d, v, c] * pb.params.norvegienne_capacity
                    )
                    / pb.params.demi_palette_capacity
                ),
            )

            if demi_palettes[d, v, c, ProductType.S] != new_n_demi_palettes_s:
                demi_palettes[d, v, c, ProductType.S] = new_n_demi_palettes_s

        # Remove any useless palettes for f products
        for c in tour:
            if (d, v, c, ProductType.F) not in demi_palettes:
                continue

            new_n_demi_palettes = max(
                0,
                ceil(
                    2
                    * delivers[d, v, c, ProductType.F]
                    / pb.params.demi_palette_capacity
                ),
            )

            if demi_palettes[d, v, c, ProductType.F] != new_n_demi_palettes:
                demi_palettes[d, v, c, ProductType.F] = new_n_demi_palettes

        # Remove any useless palettes for a products
        for c in tour:
            if (d, v, c) not in palettes:
                continue

            new_n_palettes = max(
                0,
                ceil(delivers[d, v, c, ProductType.A] / pb.params.max_palette_capacity),
            )

            if palettes[d, v, c] != new_n_palettes:
                palettes[d, v, c] = new_n_palettes


def make_solution(pb: Problem, vars: VarContainer, solver: cp_model.CpSolver):
    tours_flat = read_tours(pb, vars, solver)
    delivers = {k: solver.value(v) for k, v in vars.delivers.items()}
    palettes = {k: solver.value(v) for k, v in vars.palettes.items()}
    demi_palettes = {k: solver.value(v) for k, v in vars.demi_palettes.items()}
    norvegiennes = {k: solver.value(v) for k, v in vars.norvegiennes.items()}

    tour_duration = {
        k: solver.value(v) for k, v in vars.tour_duration.items() if k in tours_flat
    }

    post_process(
        pb,
        tours_flat,
        delivers,
        palettes,
        demi_palettes,
        norvegiennes,
        tour_duration,
    )

    tours = {}
    for (d, v), tour in tours_flat.items():
        if len(tour) <= 1:
            continue
        t = []
        trip = 0
        for node in tour[1:-1]:
            if node >= pb.n_centres:
                name = pb.pdrs[node - pb.n_centres].name
            else:
                name = pb.centres[node].name
            stop = Stop(node, name, StopType(node >= pb.n_centres))

            if node == 0:
                trip += 1
            elif node < pb.n_centres:
                stop.delivery = (
                    solver.value(delivers[d, v, trip, node, ProductType.A]),
                    solver.value(delivers[d, v, trip, node, ProductType.F]),
                    solver.value(delivers[d, v, trip, node, ProductType.S]),
                )
                stop.palettes = (
                    solver.value(palettes[d, v, trip, node]),
                    0.5 * solver.value(demi_palettes[d, v, trip, node, ProductType.F]),
                    0.5 * solver.value(demi_palettes[d, v, trip, node, ProductType.S]),
                )

                stop.norvegiennes = solver.value(norvegiennes[d, v, trip, node])
            t.append(stop)
        tours[d, v] = t

    sol = Solution(
        pb.params.week.value,
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


def create_tour_variables(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    """Creates arcs and visits variables and loads them into the VarContainer and the model"""
    n_nodes = pb.n_nodes
    used_nodes = pb.used_nodes
    allowed_vehicles = pb.allowed_vehicles

    # Optionally prune some long arcs to make solving easier
    pruned = prune_arcs(n_nodes, pb.distance_matrix, limit=1)
    if pruned:
        print(f"    Pruned arcs : {100*len(pruned)/(n_nodes)**2:.2f}%")

    vars.used = {v: model.new_bool_var("") for v in allowed_vehicles}

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
                    for node in pb.delivered_centres
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


def create_load_variables(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    """Creates the deliveries, palettes and norvegiennes variables and loads them into the VarContainer and the model"""

    ########## Load variables & constraints ###########
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            for trip in range(pb.params.max_trips):
                for c in pb.delivered_centres:
                    for t in ProductType:
                        max_demands = max(
                            scenario[t] for scenario in pb.centres[c].demands
                        )
                        vars.delivers[d, v, trip, c, t] = (
                            model.new_int_var(
                                0,
                                min(pb.vehicles[v].capacity, max_demands),
                                "",
                            )
                            if t in pb.vehicles[v].can_carry
                            else model.new_constant(0)
                        )

                    vars.palettes[d, v, trip, c] = model.new_int_var(
                        0, pb.vehicles[v].size, ""
                    )

                    vars.demi_palettes[d, v, trip, c, ProductType.F] = (
                        model.new_int_var(0, 2 * pb.vehicles[v].size, "")
                        if ProductType.F in pb.vehicles[v].can_carry
                        else model.new_constant(0)
                    )

                    vars.demi_palettes[d, v, trip, c, ProductType.S] = (
                        model.new_int_var(0, 2 * pb.vehicles[v].size, "")
                        if ProductType.S in pb.vehicles[v].can_carry
                        and pb.vehicles[v].allows_isotherm_cover
                        else model.new_constant(0)
                    )

                    vars.norvegiennes[d, v, trip, c] = (
                        model.new_int_var(0, pb.params.n_norvegiennes, "")
                        if ProductType.S in pb.vehicles[v].can_carry
                        else model.new_constant(0)
                    )


def create_duration_variables(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    delivered_centres = pb.delivered_centres
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


def add_load_constraints(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    """Adds the deliveries, palettes and norvegiennes constraints to the model"""
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            for trip in range(pb.params.max_trips):
                for c in pb.delivered_centres:
                    # Use an appropriate number of palettes for each type of product
                    model.add(
                        vars.palettes[d, v, trip, c] * pb.params.max_palette_capacity
                        >= vars.delivers[d, v, trip, c, ProductType.A]
                    )
                    # Use twice the required number of demi palettes for F - this is for fruits & vegetables
                    model.add(
                        vars.demi_palettes[d, v, trip, c, ProductType.F]
                        * pb.params.demi_palette_capacity
                        >= 2 * vars.delivers[d, v, trip, c, ProductType.F]  # type: ignore
                    )

                    model.add(
                        vars.demi_palettes[d, v, trip, c, ProductType.S]
                        * pb.params.demi_palette_capacity
                        + vars.norvegiennes[d, v, trip, c]
                        * pb.params.norvegienne_capacity
                        >= vars.delivers[d, v, trip, c, ProductType.S]  # type: ignore
                    )

                    # Not visited -> delivery is 0
                    model.add(vars.palettes[d, v, trip, c] == 0).only_enforce_if(
                        vars.visits[d, v, trip, c].negated()
                    )
                    for t in [ProductType.F, ProductType.S]:
                        model.add(
                            vars.demi_palettes[d, v, trip, c, t] == 0
                        ).only_enforce_if(vars.visits[d, v, trip, c].negated())
                    model.add(vars.norvegiennes[d, v, trip, c] == 0).only_enforce_if(
                        vars.visits[d, v, trip, c].negated()
                    )
                    for t in ProductType:
                        model.add(vars.delivers[d, v, trip, c, t] == 0).only_enforce_if(
                            vars.visits[d, v, trip, c].negated()
                        )
                        model.add(vars.delivers[d, v, trip, c, t] == 0).only_enforce_if(
                            vars.visits[d, v, trip, c].negated()
                        )
                        model.add(vars.delivers[d, v, trip, c, t] == 0).only_enforce_if(
                            vars.visits[d, v, trip, c].negated()
                        )


def add_time_windows(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
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
    for c in pb.delivered_centres:
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


def add_cover_constraints(pb: Problem, model, vars: VarContainer, loose: bool = False):
    fourgon_capacity = 1200
    for c in pb.delivered_centres:
        # Redundant Cover constraint
        max_sum_demands = max(
            sum(scenario.values()) for scenario in pb.centres[c].demands
        )

        if max_sum_demands > fourgon_capacity:
            min_visits = 2
        else:
            min_visits = 1

        # # Heuristic pruning rule (very effective)
        max_visits = ceil(max_sum_demands / pb.params.max_palette_capacity)

        # more than 3 visits is unreasonable in practice
        max_visits = min(max_visits, 3)

        if loose:
            # needed for the current tours to satisfy this constraint
            print("loose pruning")
            if c == 22:
                max_visits = 2
            if c == 23:
                max_visits = 3
            if c == 15:
                max_visits = 2
            if c == 6:
                max_visits = 2
            if c == 7:
                min_visits = 1
            if c == 16:
                max_visits = 2

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
    pb: Problem, model: cp_model.CpModel, vars: VarContainer, loose=False
):
    # Every customer must be served his demand
    for c in pb.delivered_centres:
        if loose:
            if c in [
                10,  # Fonsorbes AFS (+ Muret FS + Portet FS) trop pour 1 seul CF
                17,  # Pibrac : F/S trop pour 1 seul CF
                26,  # Les Arenes : F/S trop pour 1 seul CF
                25,  # Malepere FS + Escalquens AFS trop pour 2 CFs
            ]:  # These centres have too much demand for the current tours
                print("Relaxing ", pb.centres[c].name, pb.centres[c].demands)
                continue

        for type in ProductType:
            model.add(
                sum(
                    vars.delivers[d, v, trip, c, type]
                    for v in pb.allowed_vehicles
                    for trip in range(pb.params.max_trips)
                    for d in range(pb.n_days)
                )
                == pb.centres[c].demands[0][type]
            )


def add_robust_demand_constraints(
    pb: Problem, model: cp_model.CpModel, vars: VarContainer
):
    for c in pb.delivered_centres:
        for type in ProductType:
            model.add(
                sum(
                    vars.delivers[d, v, trip, c, type]
                    for v in pb.allowed_vehicles
                    for trip in range(pb.params.max_trips)
                    for d in range(pb.n_days)
                )
                == max(scenario[type] for scenario in pb.centres[c].demands)
            )


def add_stochastic_demand_constraints(
    pb: Problem, model: cp_model.CpModel, vars: VarContainer, percentage=1.0
):
    satisfied = [model.new_bool_var("") for _ in range(len(pb.centres[1].demands))]
    for c in pb.delivered_centres:
        for i, scenario in enumerate(pb.centres[c].demands):
            for type in ProductType:
                model.add(
                    sum(
                        vars.delivers[d, v, trip, c, type]
                        for v in pb.allowed_vehicles
                        for trip in range(pb.params.max_trips)
                        for d in range(pb.n_days)
                    )
                    >= scenario[type]
                ).only_enforce_if(satisfied[i])

    model.add(sum(satisfied) >= ceil(len(satisfied) * percentage))


def add_capacity_constraints(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    delivered_centres = pb.delivered_centres
    allowed_vehicles = pb.allowed_vehicles
    # Capacity constraints
    for d in range(pb.n_days):
        model.add(
            sum(
                vars.norvegiennes[d, v, trip, c]
                for v in allowed_vehicles
                for trip in range(pb.params.max_trips)
                for c in delivered_centres
            )
            <= pb.params.n_norvegiennes
        )

        for v in allowed_vehicles:
            for trip in range(pb.params.max_trips):
                # Delivery
                model.add(
                    sum(
                        vars.delivers[d, v, trip, c, t]
                        for c in delivered_centres
                        for t in ProductType
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
                        2 * vars.palettes[d, v, trip, c]
                        + vars.demi_palettes[d, v, trip, c, ProductType.F]
                        + vars.demi_palettes[d, v, trip, c, ProductType.S]  # type: ignore
                        for c in delivered_centres
                    )
                    <= 2 * pb.vehicles[v].size
                )

                # Two palettes per pickup
                model.add(
                    sum(2 * vars.visits[d, v, trip, p] for p in pb.pdr_nodes)
                    <= pb.vehicles[v].size
                )


def add_duration_constraints(
    pb: Problem, model: cp_model.CpModel, vars: VarContainer, loose=False
):
    delivered_centres = pb.delivered_centres
    pdr_nodes = pb.pdr_nodes

    # Duration constraints
    for d in range(pb.n_days):
        for v in pb.allowed_vehicles:
            if loose and (d, v) in [(1, 0), (1, 7), (3, 2)]:
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


def add_specific_requirements(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    """
    Adds a few "hard-coded" constraints
    These are not written in the input files but have been specified by the client
    """

    # The PL must visit Carrefour En Jacca on Wednesday
    # He must only do one pickup
    carrefour_en_jacca_index = pb.n_centres + 6
    model.add(vars.visits[2, 0, 0, carrefour_en_jacca_index] == 1)
    model.add_bool_and(
        vars.visits[2, 0, 0, p].negated()
        for p in pb.pdr_nodes
        if p != carrefour_en_jacca_index
    )

    # A camion Frigo must visit Carrefour la menude on Wednesday and Friday
    # To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    carrefour_menude_index = pb.n_centres + 5
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

    # # # Redundant
    # model.add_bool_and(
    #     vars.visits[d, v, carrefour_en_jacca_index].negated()
    #     for d in range(pb.n_days)
    #     for v in range(pb.m)
    #     if (d, v, carrefour_en_jacca_index) in vars.visits and (d, v) != (2, 0)
    # )

    # The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = 24
    model.add_bool_and(
        vars.visits[d, 0, trip, index_gde_bretagne].negated()
        for d in range(pb.n_days)
        for trip in range(pb.params.max_trips)
    )

    # Bessières cannot be delivered with camions frigos
    index_bessieres = 3
    model.add_bool_and(
        vars.visits[d, v, trip, index_bessieres].negated()
        for d in range(pb.n_days)
        for v in pb.allowed_vehicles
        for trip in range(pb.params.max_trips)
        if ProductType.F in pb.vehicles[v].can_carry
        and (d, v, trip, index_bessieres) in vars.visits
    )

    # St Orens must be picked up with a specific CF (the smallest one, for height limit reasons)
    index_st_orens = pb.n_centres + 1
    model.add_bool_and(
        vars.visits[d, v, 0, index_st_orens].negated()
        for d in range(pb.n_days)
        for v in pb.allowed_vehicles
        if v != 2
    )

    # "Livraisons de Ramasses" these deliver the food gathered in previous pickups. Only a single delivery can be done each time.
    for d, v, trip, c in pb.livraisons_de_ramasses:
        model.add_bool_and(
            vars.visits[d, v, trip, c2]
            if c2 == c
            else vars.visits[d, v, trip, c2].negated()
            for c2 in pb.delivered_centres
            if (d, v, trip, c2) in vars.visits
        )
        for t in ProductType:
            model.add(vars.delivers[d, v, trip, c, t] == 0)


def add_domsym_breaking(
    pb: Problem,
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
            sum(vars.delivers[d, v, trip, c, t] for t in ProductType) > 0
        ).only_enforce_if(vars.visits[d, v, trip, c])


def add_objective(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
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


def add_decision_strategies(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    return


def add_manual_pruning(pb: Problem, model: cp_model.CpModel, vars: VarContainer):
    model.add(sum(vars.used.values()) <= 8)
    model.add(vars.total_distance <= 2_500_000)
    return


def solve_vrp(
    pb: Problem,
    hint=None,
    time_limit=None,
    stop_at_first_solution=False,
    fix_hint=False,
):
    """Creates a CP model and solves it with cp-sat.
    When done, re-optimizes the solution and writes it in JSON format to the specified file
    """

    print("[CP-SAT] Modeling... ", end="")

    model = cp_model.CpModel()
    vars = VarContainer()

    loose = False

    print("\r[CP-SAT] Creating variables...", end="")
    create_tour_variables(pb, model, vars)
    create_load_variables(pb, model, vars)
    create_duration_variables(pb, model, vars)

    print("\r[CP-SAT] Adding constraints 0/10        ", end="")
    add_time_windows(pb, model, vars)
    add_cover_constraints(pb, model, vars, loose=loose)
    add_duration_constraints(pb, model, vars, loose=loose)

    print("\r[CP-SAT] Adding constraints 3/10        ", end="")
    add_load_constraints(pb, model, vars)
    add_capacity_constraints(pb, model, vars)
    # add_demand_constraints(pb, model, vars, loose=loose)
    add_robust_demand_constraints(pb, model, vars)
    # add_stochastic_demand_constraints(pb, model, vars, percentage=0.6)

    print("\r[CP-SAT] Adding constraints 6/10        ", end="")

    add_specific_requirements(pb, model, vars)
    add_domsym_breaking(pb, model, vars)
    add_decision_strategies(pb, model, vars)

    print("\r[CP-SAT] Adding objective        ", end="")
    add_objective(pb, model, vars)
    add_manual_pruning(pb, model, vars)

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
    if hint and fix_hint:
        solver.parameters.fix_variables_to_their_hinted_value = True

    solver.parameters.num_workers = 16
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
        sol = make_solution(pb, vars, solver)
        return status, sol
