from ortools.sat.python import cp_model
from math import inf


class LogPrinter(cp_model.CpSolverSolutionCallback):
    """Print the current objective value and the optimality gap as the search progresses"""

    def __init__(self):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.best_solution = inf

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
                f"[{self.WallTime():.2f}s]  obj = {obj/1000:<7.2f}km     gap = {gap:.2f}%"
            )


def add_domsym_breaking(
    model,
    n,
    m,
    n_days,
    delivers,
    visits,
):
    # Dominance breaking
    for v in range(m):
        for d in range(n_days):
            for c in range(1, n):
                # Visited -> delivery is not 0
                model.Add(
                    delivers["a"][v, d, c]
                    + delivers["f"][v, d, c]
                    + delivers["s"][v, d, c]
                    > 0
                ).OnlyEnforceIf(visits[v, d, c])

    # # Symmetry breaking
    # for d in range(1, n_days):
    #     for v in range(m):
    #         model.AddImplication(visits[v, d, 0], visits[v, d - 1, 0])
    return


def add_hints(model: cp_model.CpModel, visits, arcs, hint, n, n_pdr, m, n_days):
    print("Adding hints...", end="")
    hint_tours, hint_arcs = hint
    for v in range(m):
        for d in range(n_days):
            if (v, d) in hint_tours.keys():
                model.AddHint(visits[v, d, 0], 0)
                # model.Add(visits[v, d, 0] == 1)
                for node in range(1, n + n_pdr):
                    # # Add provided tours as constraints
                    # model.Add(visits[v, d, node] == (node in hint_tours[v, d]))
                    # Add them as hints
                    model.AddHint(visits[v, d, node], (node in hint_tours[v, d]))

            vd = v + d * m
            for node1 in range(n + n_days):
                for node2 in range(n + n_days):
                    if node1 == node2:
                        continue
                    if (vd, node1, node2) in hint_arcs.keys() and (
                        vd,
                        node1,
                        node2,
                    ) in arcs.keys():
                        model.AddHint(
                            arcs[vd, node1, node2][2], hint_arcs[vd, node1, node2]
                        )
                        # model.Add(
                        #     arcs[vd, node1, node2][2] == hint_arcs[vd, node1, node2]
                        # )
    print("\r")


def prune_arcs(n_nodes, matrix, limit=0.28):
    max_dist = matrix.max().max()
    pruned = []
    for node1 in range(n_nodes):
        for node2 in range(n_nodes):
            if matrix.iloc[node1, node2] >= limit * max_dist:
                pruned.append((node1, node2))
    return pruned


def solve_vrp(
    matrix,
    n,
    n_pdr,
    m,
    n_days,
    demands,
    freqs_pdr,
    capacities,
    sizes,
    frais,
    max_palette_capacity,
    hint=None,
):
    print("Modeling...")

    # TODO : dispatch semi-hebdomadaire to the appropriate week

    model = cp_model.CpModel()

    pruned = prune_arcs(n + n_pdr, matrix)
    print(f"    Pruned arcs : {100*len(pruned)/(n+n_pdr)**2:.2}%")

    ########## Arcs & Visit variables ###########
    visits = {}
    arcs = {}
    for v in range(m):
        for d in range(n_days):
            # Identifies uniquely the vehicle-day pair (easier for arcs)
            vd = v + d * m
            for node in range(n + n_pdr):
                visits[v, d, node] = model.NewBoolVar("visits_%i_%i_%i" % (v, d, node))

                # Add an arc from the node to itself (necessary for Circuit constraints)
                # If used, the node is not visited
                arcs[vd, node, node] = (node, node, visits[v, d, node].Not())

                for node2 in range(n + n_pdr):
                    if node == node2:
                        continue
                    if (node, node2) in pruned:
                        continue
                    # No arc from pdr to centre (except depot)
                    elif node >= n and 0 < node2 and node2 < n:
                        continue
                    else:
                        lit = model.NewBoolVar("arc_%i_%i_%i" % (vd, node, node2))
                        arcs[vd, node, node2] = (node, node2, lit)

            # If you don't visit the depot, visit nothing
            model.AddBoolAnd(visits[v, d, c].Not() for c in range(n)).OnlyEnforceIf(
                visits[v, d, 0].Not()
            )
            # for c in range(1, n):
            #     model.AddImplication(visits[v, 0].Not(), visits[v, c].Not())

            model.AddCircuit(
                [
                    arcs[vd, node1, node2]
                    for node1 in range(n + n_pdr)
                    for node2 in range(n + n_pdr)
                    if (node2 == 0 or node1 < n or node2 >= n)
                    and (node1, node2) not in pruned
                ]
            )

    ########## Load variables & constraints ###########
    delivers = {
        "a": {},
        "f": {},
        "s": {},
    }
    palettes = {}
    for v in range(m):
        for d in range(n_days):
            for c in range(n):
                delivers["a"][v, d, c] = model.NewIntVar(
                    0,
                    min(max(capacities), max(demands["a"])),
                    f"delivers_a_{v}_{d}_{c}",
                )

                # Frais must be delivered with appropriate vehicles
                if frais[v]:
                    delivers["f"][v, d, c] = model.NewIntVar(
                        0,
                        min(max(capacities), max(demands["f"])),
                        "delivers_f_%i_%i_%i" % (v, d, c),
                    )
                else:
                    delivers["f"][v, d, c] = 0

                delivers["s"][v, d, c] = model.NewIntVar(
                    0,
                    min(max(capacities), max(demands["s"])),
                    "delivers_s_%i_%i_%i" % (v, d, c),
                )

                palettes[v, d, c] = model.NewIntVar(
                    0, sizes[v], "palettes_%i_%i_%i" % (v, d, c)
                )

                # Use an appropriate number of palettes
                model.Add(
                    palettes[v, d, c] * max_palette_capacity
                    >= delivers["a"][v, d, c] + delivers["f"][v, d, c]
                )

                # Not visited -> delivery is 0
                model.Add(delivers["a"][v, d, c] == 0).OnlyEnforceIf(
                    visits[v, d, c].Not()
                )
                model.Add(delivers["f"][v, d, c] == 0).OnlyEnforceIf(
                    visits[v, d, c].Not()
                )
                model.Add(delivers["s"][v, d, c] == 0).OnlyEnforceIf(
                    visits[v, d, c].Not()
                )

    for pdr in range(n, n + n_pdr):
        model.Add(
            sum(visits[v, d, pdr] for v in range(m) for d in range(n_days))
            >= freqs_pdr[pdr - n]
        )

    # Add hints
    if hint:
        add_hints(model, visits, arcs, hint, n, n_pdr, m, n_days)

    # Every customer must be served exactly his demand
    for c in range(1, n):
        model.Add(
            sum(delivers["a"][v, d, c] for v in range(m) for d in range(n_days))
            == demands["a"][c]
        )
        model.Add(
            sum(delivers["f"][v, d, c] for v in range(m) for d in range(n_days))
            == demands["f"][c]
        )
        model.Add(
            sum(delivers["s"][v, d, c] for v in range(m) for d in range(n_days))
            == demands["s"][c]
        )

    # add_domsym_breaking(model, n, m, n_days, delivers, visits)

    # Capacity constraints
    for v in range(m):
        for d in range(n_days):
            model.Add(
                sum([delivers["a"][v, d, c] + delivers["f"][v, d, c] for c in range(n)])
                <= capacities[v]
            )

            model.Add(sum([palettes[v, d, c] for c in range(n)]) <= sizes[v])

    # Minimize total distance
    total_distance = sum(
        arc[2] * matrix.iloc[arc[0], arc[1]]
        for arc in arcs.values()
        if arc[0] != arc[1]
    )
    model.Minimize(total_distance)

    solver = cp_model.CpSolver()

    solver.parameters.num_workers = 12
    # solver.parameters.max_time_in_seconds = 180
    # solver.parameters.min_num_lns_workers = 8

    solver.parameters.log_search_progress = True
    # status = solver.Solve(model, LogPrinter())
    status = solver.Solve(model)

    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        print("No solution found")
        exit()

    tours = {}
    obj = solver.Value(total_distance) / 1000
    for d in range(n_days):
        for v in range(m):
            tour = [0]
            vd = v + d * m

            while True:
                a = tour[-1]
                found = False
                for b in range(n + n_pdr):
                    if a == b:
                        continue
                    if a >= n and b <= n and b > 0:
                        continue
                    if (vd, a, b) not in arcs:
                        continue
                    arc = arcs[vd, a, b]
                    # print(v, "\t", arc, a, b, solver.Value(arc[2]))
                    if solver.Value(arc[2]):
                        tour.append(b)
                        found = True
                        break

                if not found:
                    break
                if tour[-1] == 0:
                    break

            tours[v, d] = tour

    deliveries = (
        {k: solver.Value(v) for k, v in delivers["a"].items()},
        {k: solver.Value(v) for k, v in delivers["f"].items()},
        {k: solver.Value(v) for k, v in delivers["s"].items()},
    )

    return (
        tours,
        obj,
        deliveries,
        {k: solver.Value(v) for k, v in visits.items()},
        {k: solver.Value(v[2]) for k, v in arcs.items()},
    )
