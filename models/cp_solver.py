from ortools.sat.python import cp_model
from math import inf, ceil
import json
from utils.problem import Problem, Solution, DeliveryWeek, ProductType, Stop, StopType


class LogPrinter(cp_model.CpSolverSolutionCallback):
    """Print the current objective value and the optimality gap as the search progresses"""

    def __init__(self, fixed_costs, variable_costs, total_distance, used):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.best_solution = inf
        self.fixed_costs = fixed_costs
        self.variable_costs = variable_costs
        self.total_distance = total_distance
        self.used = used

    def on_solution_callback(self):
        obj = self.ObjectiveValue()

        if obj < self.best_solution:
            self.best_solution = obj
            if self.BestObjectiveBound() == 0:
                return
            gap = (
                100
                * (self.best_solution - self.BestObjectiveBound())
                / self.BestObjectiveBound()
            )
            print(
                f"{f'[{self.WallTime():.2f}s]':8}  cost={obj:>6.2f}€ ({self.Value(self.fixed_costs):.0f}F {self.Value(self.variable_costs):.0f}V) |"
                f" n_vehicles={self.Value(sum(self.used.values()))} | dist={self.Value(self.total_distance)/1000:.1f}km   gap = {gap:.2f}% [{self.BestObjectiveBound():.2f}, {obj:.2f}]"
            )


def add_domsym_breaking(
    pb: Problem,
    model,
    delivers,
    visits,
):
    # Dominance breaking
    for d, v, c in visits:

        if c == 0 or c >= pb.n_centres:
            continue

        # Visited -> delivery is not 0
        model.Add(sum(delivers[d, v, c, t] for t in ProductType) > 0).OnlyEnforceIf(
            visits[d, v, c]
        )

    # # Symmetry breaking
    for d in range(1, pb.n_days):
        for v in range(3, 9):
            if not pb.vehicles[v].allowed:
                continue
            model.AddImplication(visits[d, v, 0], visits[d, v - 1, 0])


def add_hint(
    model: cp_model.CpModel,
    problem: Problem,
    arcs,
    visits,
    delivers,
    palettes,
    demi_palettes,
    norvegiennes,
    init_sol: Solution,
) -> None:
    """
    Adds an initial solution as a hint to the model.

    Args:
        model (cp_model.CpModel): The constraint programming model.
        problem (Problem): The problem instance.
        arcs: The arcs in the problem instance.
        visits: The visits in the problem instance.
        delivers: The delivers in the problem instance.
        palettes: The palettes in the problem instance.
        demi_palettes: The demi_palettes in the problem instance.
        norvegiennes: The norvegiennes in the problem instance.
        init_sol: The file containing the hint data.
    """

    # Iterate over each day and vehicle
    for day in range(problem.n_days):
        for vehicle in range(problem.m):

            # Ignore unused vehicles
            if not problem.vehicles[vehicle].allowed:
                continue

            # If the hint for the current day and vehicle is missing or empty, hint the solver to not use the vehicle on this day
            if (day, vehicle) not in init_sol.tours:
                if problem.vehicles[vehicle].allowed:
                    model.AddHint(visits[day, vehicle, 0], 0)
                continue

            # Create a set of arcs to be assigned
            to_assign = set(
                (start, end)
                for start in range(problem.n_centres + problem.n_pdr)
                for end in range(problem.n_centres + problem.n_pdr)
                if start != end and (day, vehicle, start, end) in arcs
            )

            # hint the solver to use this vehicle on this day
            if (day, vehicle, 0) in visits:
                model.AddHint(visits[day, vehicle, 0], 1)
            current_location = 0

            # Iterate over each place in the hint tour
            tour = init_sol.tours[day, vehicle]
            for stop in tour:
                next_location = stop.index

                # Set the hint for the arc from the current location to the next location
                if (day, vehicle, current_location, next_location) in arcs:
                    model.AddHint(
                        arcs[day, vehicle, current_location, next_location][2], 1
                    )
                    to_assign.remove((current_location, next_location))
                current_location = next_location

                # Set the hints for the quantities to be delivered at the next location
                if (
                    next_location < problem.n_centres
                    and (day, vehicle, next_location) in palettes
                    and problem.centres[next_location].delivery_week
                    in [DeliveryWeek.ANY, problem.params.week]
                ):
                    for product_type in ProductType:
                        if product_type in problem.vehicles[vehicle].can_carry:
                            model.AddHint(
                                delivers[day, vehicle, next_location, product_type],
                                stop.delivery[product_type.value],
                            )

                    if ProductType.A in problem.vehicles[vehicle].can_carry:
                        model.AddHint(
                            palettes[day, vehicle, next_location],
                            ceil(stop.palettes[0]),
                        )

                    if ProductType.F in problem.vehicles[vehicle].can_carry:
                        model.AddHint(
                            demi_palettes[day, vehicle, next_location, ProductType.F],
                            round(stop.palettes[1] * 2),
                        )

                    if ProductType.S in problem.vehicles[vehicle].can_carry:
                        model.AddHint(
                            norvegiennes[day, vehicle, next_location], stop.norvegiennes
                        )
                        if problem.vehicles[vehicle].allows_isotherm_cover:
                            model.AddHint(
                                demi_palettes[
                                    day, vehicle, next_location, ProductType.S
                                ],
                                round(stop.palettes[2] * 2),
                            )

            # Set the hint for the arc from the current location back to the depot
            if (day, vehicle, current_location, 0) in arcs:
                model.AddHint(arcs[day, vehicle, current_location, 0][2], 1)
                to_assign.remove((current_location, 0))

            # Hint to not use any remaining arcs
            for start, end in to_assign:
                if (day, vehicle, start, end) in arcs:
                    model.AddHint(arcs[day, vehicle, start, end][2], 0)


def prune_arcs(n_nodes, matrix, limit=1.0):
    """Evaluates whether each arc is longer than a given ratio of the maximum distance between nodes in the matrix, and returns the list of arcs that should be pruned"""
    max_dist = matrix.max().max()
    pruned = []
    for node1 in range(n_nodes):
        for node2 in range(n_nodes):
            if matrix.iloc[node1, node2] > limit * max_dist:
                pruned.append((node1, node2))
    return pruned


def read_tours(pb: Problem, arcs, solver):
    """Reads tours from the solver
    Returns a dictionary of tours with keys (day, vehicle)"""
    tours = {}
    for d in range(pb.n_days):
        for v in range(pb.m):
            tour = [0]

            while True:
                a = tour[-1]
                goes_to = set()

                for b in range(1, pb.n_centres + pb.n_pdr):
                    if a == b:
                        continue
                    if (d, v, a, b) not in arcs:
                        continue
                    arc = arcs[d, v, a, b]

                    if solver.Value(arc[2]):
                        goes_to.add(b)

                if len(goes_to) > 1:
                    print("Warning : ", d, v, a, "goes to ", goes_to)
                if len(goes_to) == 0:
                    break

                tour.append(goes_to.pop())
                if tour[-1] == 0:
                    break

            tour.append(0)

            if len(tour) == 2:
                continue
            tours[d, v] = tour
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
                    delivers[d, v, c, ProductType.F] / pb.params.demi_palette_capacity
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

    for d, v in tour_duration:
        tour_duration[d, v] += pb.duration_matrix.iloc[tours[d, v][-2], 0]


def add_specific_requirements(pb: Problem, model: cp_model.CpModel, visits, arcs):
    """
    Adds a few "hard-coded" constraints
    These are not written in the input files but have been specified by the client
    """

    # The PL must visit Carrefour Centrale on Wednesday
    # He must only do one pickup
    carrefour_centrale_index = pb.n_centres + 5
    carrefour_centrale_a_index = pb.n_centres + 6
    model.Add(visits[2, 0, carrefour_centrale_a_index] == 1)
    model.AddBoolAnd(
        visits[2, 0, p].Not()
        for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
        if p != carrefour_centrale_a_index
    )

    # The first camion Frigo must visit Carrefour Centrale on Wednesday and Friday
    # To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    model.Add(visits[2, 2, carrefour_centrale_index] == 1)
    model.AddBoolAnd(
        visits[2, 2, p].Not()
        for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
        if p != carrefour_centrale_index
    )

    model.Add(visits[4, 2, carrefour_centrale_index] == 1)
    model.AddBoolAnd(
        visits[4, 2, p].Not()
        for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
        if p != carrefour_centrale_index
    )

    # # Redundant
    model.AddBoolAnd(
        visits[d, v, carrefour_centrale_a_index].Not()
        for d in range(pb.n_days)
        for v in range(pb.m)
        if (d, v, carrefour_centrale_a_index) in visits and (d, v) != (2, 0)
    )
    # for vd in [0 + 2 * pb.m, 2 + 2 * pb.m, 2 + 4 * pb.m]:
    #     if not arcs:
    #         break
    #     model.AddBoolAnd(
    #         (
    #             arcs[vd, p1, p2][2].Not()
    #             if p1 != p2 or p1 == carrefour_centrale_index
    #             else arcs[vd, p1, p2][2]
    #         )
    #         for p1 in range(pb.n_centres, pb.n_centres + pb.n_pdr)
    #         for p2 in range(pb.n_centres, pb.n_centres + pb.n_pdr)
    #     )
    #     model.Add(arcs[vd, carrefour_centrale_index, 0][2] == 1)

    # The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = 24
    model.AddBoolAnd(visits[d, 0, index_gde_bretagne].Not() for d in range(pb.n_days))
    # model.AddBoolAnd(
    #     arcs[0 + d * pb.m, p1, p2][2].Not()
    #     for d in range(pb.n_days)
    #     for p1 in range(pb.n_centres, pb.n_centres + pb.n_pdr)
    #     for p2 in range(pb.n_centres, pb.n_centres + pb.n_pdr)
    #     if p1 == index_gde_bretagne or p2 == index_gde_bretagne
    # )

    # A camion frigo must deliver Revel on Tuesdays and do nothing else
    index_revel = 20
    model.AddBoolAnd(
        visits[1, 2, node] if node == index_revel else visits[1, 2, node].Not()
        for node in range(1, pb.n_centres + pb.n_pdr)
        if (1, 2, node) in visits
    )
    # vd = 2 + 1 * pb.m
    # model.AddBoolAnd(
    #     (
    #         arcs[vd, a, b][2]
    #         if (a, b) in [(0, index_revel), (index_revel, 0)]
    #         else arcs[vd, a, b][2].Not()
    #     )
    #     for (vd2, a, b) in arcs
    #     if vd2 == vd and a != b
    # )
    model.AddBoolAnd(
        visits[d, v, index_revel].Not()
        for d in range(pb.n_days)
        for v in range(pb.m)
        if pb.vehicles[v].allowed and v != 2
    )

    # A Camion Frigo must deliver Les Arènes on Friday then pickup Leclerc Blagnac
    index_lc_blagnac = pb.n_centres + 3
    index_arenes = 26
    model.Add(visits[4, 3, index_arenes] == 1)
    model.Add(visits[4, 3, index_lc_blagnac] == 1)
    # model.Add(arcs[3 + 4 * pb.m, index_arenes, index_lc_blagnac][2] == 1)


def solve_vrp(
    pb: Problem,
    hint=None,
    outfile=None,
    time_limit=None,
    violation_cost=None,
    predefined_visits=None,
):
    """Creates a CP model and solves it with cp-sat.
    When done, re-optimizes the solution and writes it in JSON format to the specified file
    """

    print("Modeling...")

    model = cp_model.CpModel()

    n_nodes = pb.n_centres + pb.n_pdr
    allowed_vehicles = [v for v in range(pb.m) if pb.vehicles[v].allowed]
    delivered_centres = [
        c
        for c in range(1, pb.n_centres)
        if pb.centres[c].delivery_week in [DeliveryWeek.ANY, pb.params.week]
    ]
    pdr_nodes = range(pb.n_centres, pb.n_centres + pb.n_pdr)
    used_nodes = [0] + delivered_centres + list(pdr_nodes)

    # Optionally prune some long arcs to make solving easier
    pruned = prune_arcs(n_nodes, pb.distance_matrix, limit=1)
    if pruned:
        print(f"    Pruned arcs : {100*len(pruned)/(n_nodes)**2:.2f}%")

    ########## Arcs & Visit variables ###########
    visits = {}
    arcs = {}

    for d in range(pb.n_days):
        for v in allowed_vehicles:
            for node in used_nodes:
                if predefined_visits and node != 0:
                    visits[d, v, node] = model.NewConstant(
                        predefined_visits[d, v, node]
                    )
                else:
                    visits[d, v, node] = model.NewBoolVar(
                        "visits_%i_%i_%i" % (d, v, node)
                    )

                # Add an arc from the node to itself (necessary for Circuit constraints)
                # If used, the node is not visited
                arcs[d, v, node, node] = (node, node, visits[d, v, node].Not())

                for node2 in used_nodes:
                    if node == node2:
                        continue
                    if (node, node2) in pruned:
                        continue
                    # No arc from pdr to centre (except depot) - we do deliveries first then pickups
                    elif node >= pb.n_centres and 0 < node2 and node2 < pb.n_centres:
                        continue
                    else:
                        lit = model.NewBoolVar("")
                        arcs[d, v, node, node2] = (node, node2, lit)

            # If you don't visit the depot, visit nothing
            model.AddBoolAnd(
                visits[d, v, node].Not()
                for node in delivered_centres
                + list(range(pb.n_centres, pb.n_centres + pb.n_pdr))
            ).OnlyEnforceIf(visits[d, v, 0].Not())

            model.AddCircuit(
                [
                    arcs[d, v, node1, node2]
                    for node1 in used_nodes
                    for node2 in used_nodes
                    if (d, v, node1, node2) in arcs
                ]
            )

    ########## Load variables & constraints ###########
    delivers = {}
    palettes = {}
    demi_palettes = {}
    norvegiennes = {}
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            for c in delivered_centres:

                for t in ProductType:
                    delivers[d, v, c, t] = (
                        model.NewIntVar(
                            0,
                            min(pb.vehicles[v].capacity, pb.centres[c].demands[t]),
                            f"",
                        )
                        if t in pb.vehicles[v].can_carry
                        else model.NewConstant(0)
                    )

                palettes[d, v, c] = model.NewIntVar(0, pb.vehicles[v].size, "")

                demi_palettes[d, v, c, ProductType.F] = (
                    model.NewIntVar(0, 2 * pb.vehicles[v].size, "")
                    if ProductType.F in pb.vehicles[v].can_carry
                    else model.NewConstant(0)
                )

                demi_palettes[d, v, c, ProductType.S] = (
                    model.NewIntVar(0, 2 * pb.vehicles[v].size, "")
                    if ProductType.S in pb.vehicles[v].can_carry
                    and pb.vehicles[v].allows_isotherm_cover
                    else model.NewConstant(0)
                )

                norvegiennes[d, v, c] = (
                    model.NewIntVar(0, pb.params.n_norvegiennes, "")
                    if ProductType.S in pb.vehicles[v].can_carry
                    else model.NewConstant(0)
                )

                # Use an appropriate number of palettes for each type of product
                model.Add(
                    palettes[d, v, c] * pb.params.max_palette_capacity
                    >= delivers[d, v, c, ProductType.A]
                )
                model.Add(
                    demi_palettes[d, v, c, ProductType.F]
                    * pb.params.demi_palette_capacity
                    >= delivers[d, v, c, ProductType.F]
                )
                model.Add(
                    demi_palettes[d, v, c, ProductType.S]
                    * pb.params.demi_palette_capacity
                    + norvegiennes[d, v, c] * pb.params.norvegienne_capacity
                    >= delivers[d, v, c, ProductType.S]
                )

                # Not visited -> delivery is 0
                model.Add(palettes[d, v, c] == 0).OnlyEnforceIf(visits[d, v, c].Not())
                for t in [ProductType.F, ProductType.S]:
                    model.Add(demi_palettes[d, v, c, t] == 0).OnlyEnforceIf(
                        visits[d, v, c].Not()
                    )
                model.Add(norvegiennes[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                for t in ProductType:
                    model.Add(delivers[d, v, c, t] == 0).OnlyEnforceIf(
                        visits[d, v, c].Not()
                    )
                    model.Add(delivers[d, v, c, t] == 0).OnlyEnforceIf(
                        visits[d, v, c].Not()
                    )
                    model.Add(delivers[d, v, c, t] == 0).OnlyEnforceIf(
                        visits[d, v, c].Not()
                    )

    ########## Time  windows ###########
    violations = (
        {c: model.NewBoolVar("") for c in delivered_centres} if violation_cost else {}
    )

    for d in range(pb.n_days):
        for pdr in range(pb.n_pdr):
            # Visit each exactly once on each required days with an appropriate vehicle
            if d in pb.pdrs[pdr].required_days:
                model.AddExactlyOne(
                    visits[d, v, pdr + pb.n_centres]
                    for v in allowed_vehicles
                    if pb.pdrs[pdr].product_type in pb.vehicles[v].can_carry
                )

                # Disallow any ramasses with vehicles that cannot carry that type of product
                model.AddBoolAnd(
                    visits[d, v, pdr + pb.n_centres].Not()
                    for v in allowed_vehicles
                    if pb.pdrs[pdr].product_type not in pb.vehicles[v].can_carry
                )
            else:
                # Do not visit them if not required (this improves pruning)
                model.AddBoolAnd(
                    [visits[d, v, pdr + pb.n_centres].Not() for v in allowed_vehicles]
                )

    # Only deliver centres on allowed days
    for c in delivered_centres:
        for d in range(pb.n_days):
            if d not in pb.centres[c].allowed_days:
                if violation_cost:
                    model.AddBoolAnd(
                        visits[d, v, c].Not() for v in allowed_vehicles
                    ).OnlyEnforceIf(violations[c].Not())
                else:
                    model.AddBoolAnd([visits[d, v, c].Not() for v in allowed_vehicles])

    for c in delivered_centres:

        # Redundant Cover constraint
        min_visits = max(
            ceil(sum(pb.centres[c].demands.values()) / pb.vehicles[2].capacity),
            ceil(
                (
                    2
                    * ceil(
                        pb.centres[c].demands[ProductType.A]
                        / pb.params.max_palette_capacity
                    )
                    + ceil(
                        pb.centres[c].demands[ProductType.F]
                        / pb.params.demi_palette_capacity
                    )
                )
                / pb.vehicles[2].size
            ),
            1,
        )

        # Heuristic pruning rule (very effective)
        sum_demands = sum(pb.centres[c].demands.values())
        max_visits = ceil(sum_demands / pb.params.max_palette_capacity)

        model.AddLinearConstraint(
            sum(
                visits[d, v, c]
                for d in range(pb.n_days)
                for v in allowed_vehicles
                if d in pb.centres[c].allowed_days
            ),
            min_visits,
            max_visits,
        )

    add_specific_requirements(pb, model, visits, arcs)
    add_domsym_breaking(pb, model, delivers, visits)

    # Add input file as a hint
    if hint:
        add_hint(
            model,
            pb,
            arcs,
            visits,
            delivers,
            palettes,
            demi_palettes,
            norvegiennes,
            hint,
        )

    # Every customer must be served his demand
    for c in delivered_centres:
        for type in ProductType:
            model.Add(
                sum(
                    delivers[d, v, c, type]
                    for v in allowed_vehicles
                    for d in range(pb.n_days)
                )
                == pb.centres[c].demands[type]
            )

    # Capacity constraints
    for d in range(pb.n_days):

        model.Add(
            sum(
                norvegiennes[d, v, c]
                for v in allowed_vehicles
                for c in delivered_centres
            )
            <= pb.params.n_norvegiennes
        )

        for v in allowed_vehicles:
            # Delivery
            model.Add(
                sum(
                    delivers[d, v, c, t] for c in delivered_centres for t in ProductType
                )
                <= pb.vehicles[v].capacity
            )

            # Pickup
            model.Add(
                sum(
                    pb.pdrs[p].weight * visits[d, v, p + pb.n_centres]
                    for p in range(pb.n_pdr)
                )
                <= pb.vehicles[v].capacity
            )

            # Size limits
            model.Add(
                sum(
                    2 * palettes[d, v, c]
                    + demi_palettes[d, v, c, ProductType.F]
                    + demi_palettes[d, v, c, ProductType.S]
                    for c in delivered_centres
                )
                <= 2 * pb.vehicles[v].size
            )

            # Two palettes per pickup
            model.Add(
                sum(2 * visits[d, v, p] for p in pdr_nodes) <= pb.vehicles[v].size
            )

    # Duration constraints
    route_dist = {}
    tour_duration = {}
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            vd_arcs = [
                arcs[d2, v2, a, b]
                for (d2, v2, a, b) in arcs.keys()
                if d2 == d and v2 == v
            ]

            route_dist[d, v] = sum(
                arc[2] * pb.distance_matrix.iloc[arc[0], arc[1]] for arc in vd_arcs
            )

            tour_duration[d, v] = sum(
                arc[2] * pb.duration_matrix.iloc[arc[0], arc[1]]
                for arc in vd_arcs
                if arc[1] != 0  # The trip back to the depot is not counted
            ) + 60 * (
                pb.params.wait_at_centres
                * sum(visits[d, v, c] for c in delivered_centres)
                + pb.params.wait_at_pdrs * sum(visits[d, v, p] for p in pdr_nodes)
            )

            model.Add(tour_duration[d, v] <= pb.params.max_tour_duration * 60)
            model.Add(
                sum(visits[d, v, node] for node in delivered_centres + list(pdr_nodes))
                <= pb.params.max_stops
            )

            # time_to_first_pickup = sum(
            #     arc[2] * pb.duration_matrix.iloc[arc[0], arc[1]]
            #     for arc in vd_arcs
            #     if (arc[0] < pb.n)
            # ) + pb.params.wait_at_centres * 60 * sum(
            #     visits[d, v, c] for c in delivered_centres
            # )
            # ram = model.NewBoolVar("")
            # model.Add(
            #     time_to_first_pickup <= pb.max_first_pickup_time * 60
            # ).OnlyEnforceIf(ram)

            # model.AddBoolAnd(
            #     visits[d, v, p].Not() for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
            # ).OnlyEnforceIf(ram.Not())
            # model.AddBoolOr(
            #     visits[d, v, p] for p in range(pb.n_centres, pb.n_centres + pb.n_pdr)
            # ).OnlyEnforceIf(ram)

    # Used vehicles
    used = {v: model.NewBoolVar("") for v in allowed_vehicles}
    for v in allowed_vehicles:
        model.AddBoolAnd(
            visits[d, v2, c].Not() for (d, v2, c) in visits if v2 == v
        ).OnlyEnforceIf(used[v].Not())

    # Objectives
    total_distance = sum(
        arc[2] * pb.distance_matrix.iloc[arc[0], arc[1]]
        for arc in arcs.values()
        if arc[0] != arc[1]
    )

    variable_costs = (
        sum(
            route_dist[d, v] * pb.vehicles[v].cost_per_km
            for d in range(pb.n_days)
            for v in allowed_vehicles
        )
        * 0.001
    )

    fixed_costs = sum(used[v] * pb.vehicles[v].fixed_cost for v in allowed_vehicles)

    total_costs = variable_costs + fixed_costs
    if violation_cost:
        total_costs += violation_cost * sum(violations.values())
    model.Minimize(total_costs)

    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 24

    solver.parameters.use_lns_only = True
    # solver.parameters.diversify_lns_params = True
    # solver.parameters.min_num_lns_workers = 24

    # solver.parameters.repair_hint = True
    # solver.parameters.hint_conflict_limit = 9999999
    # solver.parameters.fix_variables_to_their_hinted_value = True

    # solver.parameters.ignore_subsolvers.extend(
    #     [
    #         "lb_tree_search",
    #         "max_lp",
    #         "pseudo_costs",
    #         "reduced_costs",
    #         "objective_lb_search",
    #         "objective_lb_search_max_lp",
    #         "objective_lb_search_no_lp",
    #         "objective_shaving_search_max_lp",
    #         "objective_shaving_search_no_lp",
    #         "probing",
    #         "probing_max_lp",
    #         "quick_restart",
    #         "default_lp",
    #         "core",
    #     ]
    # )

    """
    All solvers : [
        core,
        default_lp,
        max_lp,
        no_lp,
        objective_lb_search,
        objective_lb_search_max_lp,
        objective_lb_search_no_lp,
        objective_shaving_search_max_lp,
        objective_shaving_search_no_lp,
        probing,
        probing_max_lp,
        pseudo_costs,
        quick_restart,
        quick_restart_no_lp,
        reduced_costs,
    ]
    """

    if time_limit:
        solver.parameters.max_time_in_seconds = time_limit

    print("Solving...")
    solver.parameters.log_search_progress = True
    status = solver.Solve(model)
    # status = solver.Solve(
    #     model, LogPrinter(fixed_costs, variable_costs, total_distance, used)
    # )

    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        print("No solution found")
        return [], 0

    if violation_cost:
        for c in violations:
            if solver.Value(violations[c]):
                print("Violation on centre ", c)

    if not outfile:
        return

    tours_flat = read_tours(pb, arcs, solver)
    delivers = {k: solver.Value(v) for k, v in delivers.items()}
    palettes = {k: solver.Value(v) for k, v in palettes.items()}
    demi_palettes = {k: solver.Value(v) for k, v in demi_palettes.items()}
    norvegiennes = {k: solver.Value(v) for k, v in norvegiennes.items()}

    tour_duration = {
        k: solver.Value(v) for k, v in tour_duration.items() if k in tours_flat
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
        t = []
        for node in tour[1:-1]:
            name = pb.distance_matrix.index[node]
            stop = Stop(node, name, StopType(node >= pb.n_centres))  # type:ignore
            if node < pb.n_centres:
                stop.delivery = (
                    solver.Value(delivers[d, v, node, ProductType.A]),
                    solver.Value(delivers[d, v, node, ProductType.F]),
                    solver.Value(delivers[d, v, node, ProductType.S]),
                )
                stop.palettes = (
                    solver.Value(palettes[d, v, node]),
                    0.5 * solver.Value(demi_palettes[d, v, node, ProductType.F]),
                    0.5 * solver.Value(demi_palettes[d, v, node, ProductType.S]),
                )
                stop.norvegiennes = solver.Value(norvegiennes[d, v, node])
            t.append(stop)
        tours[d, v] = t

    sol = Solution(
        pb.params.week.value,
        solver.Value(total_costs),  # type: ignore
        solver.Value(fixed_costs),  # type: ignore
        solver.Value(variable_costs),  # type: ignore
        solver.Value(total_distance),
        [
            bool(solver.Value(used[v])) if v in allowed_vehicles else False
            for v in range(pb.m)
        ],
        tours,
        {k: solver.Value(v) for k, v in tour_duration.items() if k in tours_flat},
    )
    sol.adjust_durations(pb.duration_matrix, pb.no_traffic_duration_matrix, pb.params)

    sol.write_as_json(outfile)
