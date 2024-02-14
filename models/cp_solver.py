from ortools.sat.python import cp_model
from math import inf, ceil
import json


class LogPrinter(cp_model.CpSolverSolutionCallback):
    """Print the current objective value and the optimality gap as the search progresses"""

    def __init__(self, fuel, total_distance, used):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.best_solution = inf
        self.fuel = fuel
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
                f"[{self.WallTime():.2f}s]  cost={obj:<7.2f}E n_vehicles={self.Value(sum(self.used.values()))} fuel={self.Value(self.fuel):.1f}L dist={self.Value(self.total_distance)/1000:.1f}km   gap = {gap:.2f}% [{self.BestObjectiveBound():.2f}, {obj:.2f}]"
            )


def add_domsym_breaking(
    pb,
    model,
    delivers,
    visits,
):
    # Dominance breaking
    for d, v, c in delivers["a"]:
        # Visited -> delivery is not 0
        model.Add(
            delivers["a"][d, v, c] + delivers["f"][d, v, c] + delivers["s"][d, v, c] > 0
        ).OnlyEnforceIf(visits[d, v, c])

    # # Symmetry breaking
    for d in range(1, pb.n_days):
        for v in range(3, 8):
            if not pb.vehicle_allowed[v]:
                continue
            model.AddImplication(visits[d, v, 0], visits[d, v - 1, 0])
    return


def add_hints(
    model: cp_model.CpModel,
    problem,
    arcs,
    visits,
    delivers,
    palettes,
    demi_palettes,
    demi_palettes_s,
    norvegiennes,
    hint_file,
    force=False,
):
    """Adds an initial solution as a hint to the model. If force is True, add constraints to force the solution to be the same."""

    print("Adding hints...", end="")

    hint = json.load(open(hint_file))["tours"]

    for d in range(problem.n_days):
        for v in range(problem.m):
            if not problem.vehicle_allowed[v]:
                continue
            vd = v + d * problem.m
            key = str(d) + ", " + str(v)

            if force:
                if key in hint.keys() and hint[key]:
                    to_assign = set(
                        (a, b)
                        for a in range(problem.n + problem.n_pdr)
                        for b in range(problem.n + problem.n_pdr)
                        if a != b and (vd, a, b) in arcs
                    )
                    model.Add(visits[d, v, 0] == 1)
                    a = 0
                    for place in hint[key]:
                        b = place["index"]
                        model.Add(arcs[vd, a, b][2] == 1)
                        to_assign.remove((a, b))
                        a = b

                        if b < problem.n:
                            model.Add(palettes[d, v, b] == place["palettes"][0])
                            model.Add(norvegiennes[d, v, b] == place["norvegiennes"])
                            # model.Add(delivers["a"][d, v, b] == place["delivery"][0])
                            if problem.frais[v]:
                                # model.Add(
                                #     delivers["f"][d, v, b] == place["delivery"][1]
                                # )
                                model.Add(
                                    demi_palettes[d, v, b] == place["palettes"][1]
                                )
                                model.Add(
                                    demi_palettes_s[d, v, b] == place["palettes"][2]
                                )
                            # model.Add(delivers["s"][d, v, b] == place["delivery"][2])

                    model.Add(arcs[vd, a, 0][2] == 1)
                    to_assign.remove((a, 0))

                    for a, b in to_assign:
                        model.Add(arcs[vd, a, b][2] == 0)
                else:
                    model.Add(visits[d, v, 0] == 0)
            else:
                if key in hint.keys() and hint[key]:
                    to_assign = set(
                        (a, b)
                        for a in range(problem.n + problem.n_pdr)
                        for b in range(problem.n + problem.n_pdr)
                        if a != b and (vd, a, b) in arcs
                    )
                    if (d, v, 0) in visits:
                        model.AddHint(visits[d, v, 0], 1)
                    a = 0
                    for place in hint[key]:
                        b = place["index"]
                        if (vd, a, b) in arcs:
                            model.AddHint(arcs[vd, a, b][2], 1)
                            to_assign.remove((a, b))
                        a = b

                        if (
                            b < problem.n
                            and (d, v, b) in palettes
                            and problem.demands["a"][b] > 0
                        ):
                            model.AddHint(palettes[d, v, b], place["palettes"][0])
                            if v != 0 or not problem.disallow_norvegiennes_in_PL:
                                model.AddHint(
                                    norvegiennes[d, v, b], place["norvegiennes"]
                                )
                            model.AddHint(delivers["a"][d, v, b], place["delivery"][0])
                            if problem.frais[v]:
                                model.AddHint(
                                    delivers["f"][d, v, b], place["delivery"][1]
                                )
                                model.AddHint(
                                    demi_palettes[d, v, b], place["palettes"][1]
                                )
                                model.AddHint(
                                    demi_palettes_s[d, v, b], place["palettes"][2]
                                )
                            model.AddHint(delivers["s"][d, v, b], place["delivery"][2])

                    if (vd, a, 0) in arcs:
                        model.AddHint(arcs[vd, a, 0][2], 1)
                        to_assign.remove((a, 0))

                    for a, b in to_assign:
                        if (vd, a, b) in arcs:
                            model.AddHint(arcs[vd, a, b][2], 0)
                else:
                    if problem.vehicle_allowed[v]:
                        model.AddHint(visits[d, v, 0], 0)

    print("\r")


def prune_arcs(n_nodes, matrix, limit=1):  #:limit=0.28):
    """Evaluates whether each arc is longer than a given ratio of the maximum distance between nodes in the matrix, and returns the list of arcs that should be pruned"""
    max_dist = matrix.max().max()
    pruned = []
    for node1 in range(n_nodes):
        for node2 in range(n_nodes):
            if matrix.iloc[node1, node2] > limit * max_dist:
                pruned.append((node1, node2))
    return pruned


def read_tours(pb, arcs, solver):
    """Reads tours from the solver
    Returns a dictionary of tours with keys (day, vehicle)"""
    tours = {}
    for d in range(pb.n_days):
        for v in range(pb.m):
            tour = [0]
            vd = v + d * pb.m

            while True:
                a = tour[-1]
                goes_to = set()

                for b in range(1, pb.n + pb.n_pdr):
                    if a == b:
                        continue
                    if (vd, a, b) not in arcs:
                        continue
                    arc = arcs[vd, a, b]

                    if solver.Value(arc[2]):
                        goes_to.add(b)

                if len(goes_to) > 1:
                    print("Warning : ", vd, a, "goes to ", goes_to)
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


def write_sol(
    pb,
    solver,
    total_distance,
    fuel_consumption,
    total_cost,
    tour_duration,
    vehicles_used,
    tours,
    delivers,
    palettes,
    demi_palettes,
    demi_palettes_s,
    norvegiennes,
    outfile,
):
    """Writes the solution in JSON format to the specified file"""
    sol = {
        "total_distance": total_distance,
        "fuel_consumption": round(fuel_consumption, 5),
        "total_cost": total_cost,
        "vehicles_used": vehicles_used,
    }
    sol["tours"] = {}
    sol["tour_durations"] = {
        str(d) + ", " + str(v): int(td) for (d, v), td in tour_duration.items()
    }
    for (d, v), tour in tours.items():

        t = []
        for node in tour[1:-1]:
            place = {
                "index": node,
                "name": pb.distance_matrix.index[node],
            }
            if node < pb.n:
                place["type"] = "livraison"
                place["delivery"] = [
                    solver.Value(delivers[t][d, v, node]) for t in ["a", "f", "s"]
                ]
                place["palettes"] = [
                    solver.Value(palettes[d, v, node]),
                    solver.Value(demi_palettes[d, v, node]),
                    solver.Value(demi_palettes_s[d, v, node]),
                ]
                place["norvegiennes"] = solver.Value(norvegiennes[d, v, node])
            else:
                place["type"] = "ramasse"
            t.append(place)
        sol["tours"][str(d) + ", " + str(v)] = t

    json.dump(sol, open(outfile, "w"))


def post_process(
    pb,
    tours,
    delivers,
    palettes,
    demi_palettes,
    demi_palettes_s,
    norvegiennes,
    tour_duration,
):
    """Shaves any superfluous decisions (useless norvegiennes/palettes) from the solution
    Also adds the time required to come back to the depot to the tour duration"""

    for (d, v), tour in tours.items():

        # Remove any useless norvegiennes and demi palettes for s products
        for c in tour:
            if (d, v, c) not in demi_palettes_s:
                continue

            new_n_norvegiennes = max(
                0,
                ceil(
                    (
                        delivers["s"][d, v, c]
                        - demi_palettes_s[d, v, c] * pb.demi_palette_capacity
                    )
                    / pb.norvegienne_capacity
                ),
            )
            if norvegiennes[d, v, c] != new_n_norvegiennes:
                norvegiennes[d, v, c] = new_n_norvegiennes

            new_n_demi_palettes_s = max(
                0,
                ceil(
                    (
                        delivers["s"][d, v, c]
                        - norvegiennes[d, v, c] * pb.norvegienne_capacity
                    )
                    / pb.demi_palette_capacity
                ),
            )

            if demi_palettes_s[d, v, c] != new_n_demi_palettes_s:
                demi_palettes_s[d, v, c] = new_n_demi_palettes_s

        # Remove any useless palettes for f products
        for c in tour:
            if (d, v, c) not in demi_palettes:
                continue

            new_n_demi_palettes = max(
                0,
                ceil(delivers["f"][d, v, c] / pb.demi_palette_capacity),
            )

            if demi_palettes[d, v, c] != new_n_demi_palettes:
                demi_palettes[d, v, c] = new_n_demi_palettes

        # Remove any useless palettes for a products
        for c in tour:
            if (d, v, c) not in palettes:
                continue

            new_n_palettes = max(
                0,
                ceil(delivers["a"][d, v, c] / pb.max_palette_capacity),
            )

            if palettes[d, v, c] != new_n_palettes:
                palettes[d, v, c] = new_n_palettes

    for (d, v), duration in tour_duration.items():
        tour_duration[d, v] += pb.duration_matrix.iloc[tours[d, v][-2], 0]


def add_specific_requirements(pb, model: cp_model.CpModel, visits, arcs):
    """
    Adds a few "hard-coded" constraints
    These are not written in the input files but have been specified by the client
    """

    # The PL must visit Carrefour Centrale on Wednesday
    # He must only do one pickup
    carrefour_centrale_index = pb.n + 5
    model.Add(visits[2, 0, carrefour_centrale_index] == 1)
    model.AddBoolAnd(
        visits[2, 0, p].Not()
        for p in range(pb.n, pb.n + pb.n_pdr)
        if p != carrefour_centrale_index
    )
    # vd = 0 + 2 * pb.m
    # model.Add(arcs[vd, carrefour_centrale_index, 0][2] == 1)

    # The first camion Frigo must visit Carrefour Centrale on Wednesday and Friday
    model.Add(visits[2, 2, carrefour_centrale_index] == 1)
    model.Add(visits[4, 2, carrefour_centrale_index] == 1)

    # To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    model.AddBoolAnd(
        visits[2, 2, p].Not()
        for p in range(pb.n, pb.n + pb.n_pdr)
        if p != carrefour_centrale_index
    )
    model.AddBoolAnd(
        visits[4, 2, p].Not()
        for p in range(pb.n, pb.n + pb.n_pdr)
        if p != carrefour_centrale_index
    )

    # The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = 24
    model.AddBoolAnd(visits[d, 0, index_gde_bretagne].Not() for d in range(pb.n_days))

    # A camion frigo must deliver Revel on Tuesdays and do nothing else
    index_revel = 20
    model.AddBoolAnd(
        visits[1, 2, node] if node == index_revel else visits[1, 2, node].Not()
        for node in range(1, pb.n + pb.n_pdr)
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
        if pb.vehicle_allowed[v] and v != 2
    )

    # A Camion Frigo must deliver Les Arènes on Friday then pickup Leclerc Blagnac
    index_lc_blagnac = pb.n + 3
    index_arenes = 26
    model.Add(visits[4, 3, index_arenes] == 1)
    model.Add(visits[4, 3, index_lc_blagnac] == 1)


def solve_vrp(
    pb,
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

    pruned = prune_arcs(pb.n + pb.n_pdr, pb.distance_matrix, limit=1)
    if pruned:
        print(f"    Pruned arcs : {100*len(pruned)/(pb.n+pb.n_pdr)**2:.2f}%")

    ########## Arcs & Visit variables ###########
    visits = {}
    arcs = {}
    allowed_vehicles = [v for v in range(pb.m) if pb.vehicle_allowed[v]]
    delivered_centres = [
        c
        for c in range(pb.n)
        if pb.demands["a"][c] + pb.demands["f"][c] + pb.demands["s"][c] > 0
    ]
    used_nodes = [0] + delivered_centres + list(range(pb.n, pb.n + pb.n_pdr))

    for d in range(pb.n_days):
        for v in allowed_vehicles:

            # Identifies uniquely the vehicle-day pair (easier for arcs)
            vd = v + d * pb.m
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
                arcs[vd, node, node] = (node, node, visits[d, v, node].Not())

                for node2 in used_nodes:
                    if node == node2:
                        continue
                    if (node, node2) in pruned:
                        continue
                    # No arc from pdr to centre (except depot) - we do deliveries first then pickups
                    elif node >= pb.n and 0 < node2 and node2 < pb.n:
                        continue
                    else:
                        lit = model.NewBoolVar("arc_%i_%i_%i" % (vd, node, node2))
                        arcs[vd, node, node2] = (node, node2, lit)

            # # # If you don't visit the depot, visit nothing
            model.AddBoolAnd(
                visits[d, v, node].Not()
                for node in delivered_centres + list(range(pb.n, pb.n + pb.n_pdr))
            ).OnlyEnforceIf(visits[d, v, 0].Not())

            model.AddCircuit(
                [
                    arcs[vd, node1, node2]
                    for node1 in used_nodes
                    for node2 in used_nodes
                    if (vd, node1, node2) in arcs
                ]
            )

    ########## Load variables & constraints ###########
    delivers = {
        "a": {},
        "f": {},
        "s": {},
    }
    palettes = {}
    demi_palettes = {}
    demi_palettes_s = {}
    norvegiennes = {}
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            for c in delivered_centres:

                delivers["a"][d, v, c] = model.NewIntVar(
                    0,
                    min(pb.capacities[v], pb.demands["a"][c]),
                    f"delivers_a_{v}_{d}_{c}",
                )

                palettes[d, v, c] = model.NewIntVar(
                    0, pb.sizes[v], "palettes_%i_%i_%i" % (d, v, c)
                )

                # Frais must be delivered with appropriate vehicles
                if pb.frais[v]:
                    delivers["f"][d, v, c] = model.NewIntVar(
                        0,
                        min(pb.capacities[v], pb.demands["f"][c]),
                        "delivers_f_%i_%i_%i" % (d, v, c),
                    )
                    demi_palettes[d, v, c] = model.NewIntVar(0, 2 * pb.sizes[v], "")
                    demi_palettes_s[d, v, c] = model.NewIntVar(0, 2 * pb.sizes[v], "")
                else:
                    delivers["f"][d, v, c] = 0
                    demi_palettes[d, v, c] = 0
                    demi_palettes_s[d, v, c] = 0

                delivers["s"][d, v, c] = model.NewIntVar(
                    0,
                    min(pb.capacities[v], pb.demands["s"][c]),
                    "delivers_s_%i_%i_%i" % (d, v, c),
                )

                if pb.disallow_norvegiennes_in_PL and v == 0:
                    norvegiennes[d, v, c] = model.NewConstant(0)
                else:
                    norvegiennes[d, v, c] = model.NewIntVar(0, pb.n_norvegiennes, "")

                # Use an appropriate number of palettes
                model.Add(
                    palettes[d, v, c] * pb.max_palette_capacity
                    >= delivers["a"][d, v, c]
                )
                model.Add(
                    demi_palettes[d, v, c] * pb.demi_palette_capacity
                    >= delivers["f"][d, v, c]
                )
                model.Add(
                    demi_palettes_s[d, v, c] * pb.demi_palette_capacity
                    + norvegiennes[d, v, c] * pb.norvegienne_capacity
                    >= delivers["s"][d, v, c]
                )

                # Not visited -> delivery is 0
                model.Add(palettes[d, v, c] == 0).OnlyEnforceIf(visits[d, v, c].Not())
                model.Add(demi_palettes[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(demi_palettes_s[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(norvegiennes[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["a"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["f"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["s"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )

    ########## Time  windows ###########
    if violation_cost:
        violations = {c: model.NewBoolVar("") for c in delivered_centres}
    for d in range(pb.n_days):
        # Visit pdrs at required days
        for pdr in range(pb.n_pdr):
            if d in pb.j_de_ramasse[pdr]:
                model.AddExactlyOne(
                    visits[d, v, pdr + pb.n] for v in allowed_vehicles if pb.frais[v]
                )
            else:
                model.AddBoolAnd(
                    [
                        visits[d, v, pdr + pb.n].Not()
                        for v in allowed_vehicles
                        if (d, v) != (2, 0)
                    ]
                )

            # No ramasses with camions not frais
            model.AddBoolAnd(
                visits[d, v, pdr + pb.n].Not()
                for v in allowed_vehicles
                if not pb.frais[v] and (d, v) != (2, 0)
            )

        # Only deliver centres on allowed days
        for c in delivered_centres:
            if d not in pb.j_de_livraison_possibles[c]:
                if violation_cost:
                    model.AddBoolAnd(
                        [visits[d, v, c].Not() for v in allowed_vehicles]
                    ).OnlyEnforceIf(violations[c].Not())
                else:
                    model.AddBoolAnd([visits[d, v, c].Not() for v in allowed_vehicles])

    add_specific_requirements(pb, model, visits, arcs)
    add_domsym_breaking(pb, model, delivers, visits)

    # # Add hints
    if hint:
        add_hints(
            model,
            pb,
            arcs,
            visits,
            delivers,
            palettes,
            demi_palettes,
            demi_palettes_s,
            norvegiennes,
            hint,
            force=False,
        )

    # Every customer must be served his demand
    for c in delivered_centres:
        model.Add(
            sum(
                delivers["a"][d, v, c]
                for v in allowed_vehicles
                for d in range(pb.n_days)
            )
            >= pb.demands["a"][c]
        )
        model.Add(
            sum(
                delivers["f"][d, v, c]
                for v in allowed_vehicles
                for d in range(pb.n_days)
            )
            >= pb.demands["f"][c]
        )
        model.Add(
            sum(
                delivers["s"][d, v, c]
                for v in allowed_vehicles
                for d in range(pb.n_days)
            )
            >= pb.demands["s"][c]
        )

    # Capacity constraints
    for d in range(pb.n_days):

        model.Add(
            sum(
                norvegiennes[d, v, c]
                for v in allowed_vehicles
                for c in delivered_centres
            )
            <= pb.n_norvegiennes
        )

        for v in allowed_vehicles:
            # Delivery
            model.Add(
                sum(
                    [
                        delivers["a"][d, v, c]
                        + delivers["f"][d, v, c]
                        + delivers["s"][d, v, c]
                        for c in delivered_centres
                    ]
                )
                <= pb.capacities[v]
            )

            # Pickup
            model.Add(
                sum([pb.weights[p] * visits[d, v, p + pb.n] for p in range(pb.n_pdr)])
                <= pb.capacities[v]
            )

            # Size limits
            model.Add(
                sum(
                    [
                        2 * palettes[d, v, c]
                        + demi_palettes[d, v, c]
                        + demi_palettes_s[d, v, c]
                        for c in delivered_centres
                    ]
                )
                <= 2 * pb.sizes[v]
            )

            # Two palettes per pickup
            model.Add(
                sum(2 * visits[d, v, p] for p in range(pb.n, pb.n + pb.n_pdr))
                <= pb.sizes[v]
            )

    # Duration constraints
    route_dist = {}
    tour_duration = {}
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            vd = v + d * pb.m

            vd_arcs = [arcs[i, a, b] for (i, a, b) in arcs.keys() if i == vd]

            route_dist[d, v] = sum(
                arc[2] * pb.distance_matrix.iloc[arc[0], arc[1]] for arc in vd_arcs
            )

            tour_duration[d, v] = (
                sum(
                    arc[2] * pb.duration_matrix.iloc[arc[0], arc[1]]
                    for arc in vd_arcs
                    if arc[1] != 0
                )
                + pb.wait_at_centres
                * 60
                * sum(visits[d, v, c] for c in delivered_centres)
                + pb.wait_at_pdrs
                * 60
                * sum(visits[d, v, p] for p in range(pb.n, pb.n + pb.n_pdr))
            )
            model.Add(tour_duration[d, v] <= pb.max_tour_duration * 60)
            model.Add(
                sum(
                    visits[d, v, node]
                    for node in delivered_centres + list(range(pb.n, pb.n + pb.n_pdr))
                )
                <= pb.max_stops
            )

            # time_to_first_pickup = sum(
            #     arc[2] * pb.duration_matrix.iloc[arc[0], arc[1]]
            #     for arc in vd_arcs
            #     if (arc[0] < pb.n)
            # ) + pb.wait_at_centres * 60 * sum(
            #     visits[d, v, c] for c in delivered_centres
            # )
            # ram = model.NewBoolVar("")
            # model.Add(
            #     time_to_first_pickup <= pb.max_first_pickup_time * 60
            # ).OnlyEnforceIf(ram)

            # model.AddBoolAnd(
            #     visits[d, v, p].Not() for p in range(pb.n, pb.n + pb.n_pdr)
            # ).OnlyEnforceIf(ram.Not())
            # model.AddBoolOr(
            #     visits[d, v, p] for p in range(pb.n, pb.n + pb.n_pdr)
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

    fuel_consumption = (
        sum(
            route_dist[d, v] * pb.consumptions[v]
            for d in range(pb.n_days)
            for v in allowed_vehicles
        )
        * 0.00001
    )

    total_cost = fuel_consumption * pb.fuel_cost + pb.weekly_fixed_cost * sum(
        used.values()
    )
    if violation_cost:
        total_cost += violation_cost * sum(violations.values())
    model.Minimize(total_cost)

    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 16
    # solver.parameters.use_lns_only = True
    # solver.parameters.min_num_lns_workers = 12

    if time_limit:
        solver.parameters.max_time_in_seconds = time_limit

    print("Solving...")
    solver.parameters.log_search_progress = True
    status = solver.Solve(model)
    # status = solver.Solve(model, LogPrinter(fuel_consumption, total_distance, used))

    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        print("No solution found")
        return [], 0

    if violation_cost:
        for c in violations:
            if solver.Value(violations[c]):
                print("Violation on centre ", c)

    tours = read_tours(pb, arcs, solver)
    total_distance = solver.Value(total_distance)
    fuel_consumption = solver.Value(fuel_consumption)
    total_cost = solver.Value(total_cost)
    total_distance = solver.Value(total_distance)
    vehicles_used = [
        solver.Value(used[v]) if v in allowed_vehicles else 0 for v in range(pb.m)
    ]
    delivers = {
        "a": {k: solver.Value(v) for k, v in delivers["a"].items()},
        "f": {k: solver.Value(v) for k, v in delivers["f"].items()},
        "s": {k: solver.Value(v) for k, v in delivers["s"].items()},
    }
    palettes = {k: solver.Value(v) for k, v in palettes.items()}
    demi_palettes = {k: solver.Value(v) for k, v in demi_palettes.items()}
    demi_palettes_s = {k: solver.Value(v) for k, v in demi_palettes_s.items()}
    norvegiennes = {k: solver.Value(v) for k, v in norvegiennes.items()}
    tour_duration = {k: solver.Value(v) for k, v in tour_duration.items() if k in tours}

    post_process(
        pb,
        tours,
        delivers,
        palettes,
        demi_palettes,
        demi_palettes_s,
        norvegiennes,
        tour_duration,
    )

    if outfile:
        write_sol(
            pb,
            solver,
            total_distance,
            fuel_consumption,
            total_cost,
            tour_duration,
            vehicles_used,
            tours,
            delivers,
            palettes,
            demi_palettes,
            demi_palettes_s,
            norvegiennes,
            outfile,
        )

    return (
        tours,
        fuel_consumption,
    )


def cluster_nodes(pb):

    allowed_vehicles = [v for v in range(pb.m) if pb.vehicle_allowed[v]]
    delivered_centres = [
        c
        for c in range(pb.n)
        if pb.demands["a"][c] + pb.demands["f"][c] + pb.demands["s"][c] > 0
    ]
    used_nodes = [0] + delivered_centres + list(range(pb.n, pb.n + pb.n_pdr))

    model = cp_model.CpModel()

    visits = {
        (d, v, node): model.NewBoolVar("")
        for d in range(pb.n_days)
        for v in allowed_vehicles
        for node in used_nodes
    }

    same = {}
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            for node1 in used_nodes:
                for node2 in used_nodes:
                    if node1 == node2:
                        continue
                    same[d, v, node1, node2] = model.NewBoolVar("")

                    model.AddBoolOr(
                        visits[d, v, node1].Not(), visits[d, v, node2].Not()
                    ).OnlyEnforceIf(same[d, v, node1, node2].Not())

                    model.Add(same[d, v, node1, node2] == 1).OnlyEnforceIf(
                        visits[d, v, node1], visits[d, v, node2]
                    )

    ########## Time  windows ###########
    for d in range(pb.n_days):
        # Visit pdrs at required days
        for pdr in range(pb.n_pdr):
            if d in pb.j_de_ramasse[pdr]:
                model.AddExactlyOne(
                    visits[d, v, pdr + pb.n] for v in allowed_vehicles if pb.frais[v]
                )
            else:
                model.AddBoolAnd(
                    [
                        visits[d, v, pdr + pb.n].Not()
                        for v in allowed_vehicles
                        if (d, v) != (2, 0)
                    ]
                )

            # No ramasses with camions not frais
            model.AddBoolAnd(
                visits[d, v, pdr + pb.n].Not()
                for v in allowed_vehicles
                if not pb.frais[v] and (d, v) != (2, 0)
            )

        # Only deliver centres on allowed days
        for c in delivered_centres:
            if d not in pb.j_de_livraison_possibles[c]:
                model.AddBoolAnd([visits[d, v, c].Not() for v in allowed_vehicles])

    ########## Load variables & constraints ###########
    delivers = {
        "a": {},
        "f": {},
        "s": {},
    }
    palettes = {}
    demi_palettes = {}
    demi_palettes_s = {}
    norvegiennes = {}
    for d in range(pb.n_days):
        for v in allowed_vehicles:
            for c in delivered_centres:

                delivers["a"][d, v, c] = model.NewIntVar(
                    0,
                    min(pb.capacities[v], pb.demands["a"][c]),
                    f"delivers_a_{v}_{d}_{c}",
                )

                palettes[d, v, c] = model.NewIntVar(
                    0, pb.sizes[v], "palettes_%i_%i_%i" % (d, v, c)
                )

                # Frais must be delivered with appropriate vehicles
                if pb.frais[v]:
                    delivers["f"][d, v, c] = model.NewIntVar(
                        0,
                        min(pb.capacities[v], pb.demands["f"][c]),
                        "delivers_f_%i_%i_%i" % (d, v, c),
                    )
                    demi_palettes[d, v, c] = model.NewIntVar(0, 2 * pb.sizes[v], "")
                    demi_palettes_s[d, v, c] = model.NewIntVar(0, 2 * pb.sizes[v], "")
                else:
                    delivers["f"][d, v, c] = 0
                    demi_palettes[d, v, c] = 0
                    demi_palettes_s[d, v, c] = 0

                delivers["s"][d, v, c] = model.NewIntVar(
                    0,
                    min(pb.capacities[v], pb.demands["s"][c]),
                    "delivers_s_%i_%i_%i" % (d, v, c),
                )

                if pb.disallow_norvegiennes_in_PL and v == 0:
                    norvegiennes[d, v, c] = model.NewConstant(0)
                else:
                    norvegiennes[d, v, c] = model.NewIntVar(0, pb.n_norvegiennes, "")

                # Use an appropriate number of palettes
                model.Add(
                    palettes[d, v, c] * pb.max_palette_capacity
                    >= delivers["a"][d, v, c]
                )
                model.Add(
                    demi_palettes[d, v, c] * pb.demi_palette_capacity
                    >= delivers["f"][d, v, c]
                )
                model.Add(
                    demi_palettes_s[d, v, c] * pb.demi_palette_capacity
                    + norvegiennes[d, v, c] * pb.norvegienne_capacity
                    >= delivers["s"][d, v, c]
                )

                # Not visited -> delivery is 0
                model.Add(palettes[d, v, c] == 0).OnlyEnforceIf(visits[d, v, c].Not())
                model.Add(demi_palettes[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(demi_palettes_s[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(norvegiennes[d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["a"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["f"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["s"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )

    add_specific_requirements(pb, model, visits, None)

    # Every customer must be served his demand
    for c in delivered_centres:
        model.Add(
            sum(
                delivers["a"][d, v, c]
                for v in allowed_vehicles
                for d in range(pb.n_days)
            )
            == pb.demands["a"][c]
        )
        model.Add(
            sum(
                delivers["f"][d, v, c]
                for v in allowed_vehicles
                for d in range(pb.n_days)
            )
            == pb.demands["f"][c]
        )
        model.Add(
            sum(
                delivers["s"][d, v, c]
                for v in allowed_vehicles
                for d in range(pb.n_days)
            )
            == pb.demands["s"][c]
        )

    # Capacity constraints
    for d in range(pb.n_days):

        model.Add(
            sum(
                norvegiennes[d, v, c]
                for v in allowed_vehicles
                for c in delivered_centres
            )
            <= pb.n_norvegiennes
        )

        for v in allowed_vehicles:
            # Delivery
            model.Add(
                sum(
                    [
                        delivers["a"][d, v, c]
                        + delivers["f"][d, v, c]
                        + delivers["s"][d, v, c]
                        for c in delivered_centres
                    ]
                )
                <= pb.capacities[v]
            )

            # Pickup
            model.Add(
                sum([pb.weights[p] * visits[d, v, p + pb.n] for p in range(pb.n_pdr)])
                <= pb.capacities[v]
            )

            # Size limits
            model.Add(
                sum(
                    [
                        2 * palettes[d, v, c]
                        + demi_palettes[d, v, c]
                        + demi_palettes_s[d, v, c]
                        for c in delivered_centres
                    ]
                )
                <= 2 * pb.sizes[v]
            )

            # Two palettes per pickup
            model.Add(
                sum(2 * visits[d, v, p] for p in range(pb.n, pb.n + pb.n_pdr))
                <= pb.sizes[v]
            )

    total_dist = sum(
        pb.distance_matrix.iloc[node1, node2] * same[d, v, node1, node2]
        for (d, v, node1, node2) in same
    )
    model.Minimize(total_dist)

    solver = cp_model.CpSolver()
    solver.parameters.log_search_progress = True
    solver.Solve(model)

    tours = {(d, v): set() for d in range(pb.n_days) for v in allowed_vehicles}
    for d, v, node in visits:
        if solver.Value(visits[d, v, node]):
            tours[d, v].add(node)

    for d in range(pb.n_days):
        print("Day ", d)
        for v in allowed_vehicles:
            if len(tours[d, v]) > 1:
                print("\tVehicle ", v, " : ", tours[d, v])

    return {k: solver.Value(v) for k, v in visits.items()}
