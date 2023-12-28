from ortools.sat.python import cp_model
from math import inf
import json


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
    for d in range(n_days):
        for v in range(m):
            for c in range(1, n):
                # Visited -> delivery is not 0
                model.Add(
                    delivers["a"][d, v, c]
                    + delivers["f"][d, v, c]
                    + delivers["s"][d, v, c]
                    > 0
                ).OnlyEnforceIf(visits[d, v, c])

    # # Symmetry breaking
    # for d in range(1, n_days):
    #     for v in range(m):
    #         model.AddImplication(visits[d, v, 0], visits[d, v - 1, 0])
    return


def add_hints(model: cp_model.CpModel, arcs, visits, problem, hint):
    # TODO : fix this
    print("Adding hints...", end="")
    hint_tours, hint_arcs = hint
    for d in range(problem.n_days):
        for v in range(problem.m):
            if (d, v) in hint_tours.keys():
                # model.AddHint(visits[d, v, 0], 0)
                model.Add(visits[d, v, 0] == 1)
                for node in range(1, problem.n + problem.n_pdr):
                    # Add provided tours as constraints
                    model.Add(visits[d, v, node] == (node in hint_tours[d, v]))
                    # Add them as hints
                    # model.AddHint(visits[d, v, node], (node in hint_tours[d, v]))
                    pass

            vd = v + d * problem.m
            for node1 in range(problem.n + problem.n_pdr):
                for node2 in range(problem.n + problem.n_pdr):
                    if node1 == node2:
                        continue
                    if (vd, node1, node2) in hint_arcs.keys() and (
                        vd,
                        node1,
                        node2,
                    ) in arcs.keys():
                        # model.AddHint(
                        #     arcs[vd, node1, node2][2], hint_arcs[vd, node1, node2]
                        # )

                        pass
                        # if hint_arcs[vd, node1, node2]:
                        #     model.Add(arcs[vd, node1, node2][2] == 1)
                        # else:
                        #     model.Add(arcs[vd, node1, node2][2] == 0)
    print("\r")


def prune_arcs(n_nodes, matrix, limit=0.28):  #:limit=0.28):
    max_dist = matrix.max().max()
    pruned = []
    for node1 in range(n_nodes):
        for node2 in range(n_nodes):
            if matrix.iloc[node1, node2] >= limit * max_dist:
                pruned.append((node1, node2))
    return pruned


def read_tours(pb, arcs, solver):
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

            tours[d, v] = tour
    return tours


def write_sol(pb, solver, obj, tours, delivers, palettes):
    sol = {"total_distance": obj}
    sol["tours"] = {}
    for (d, v), tour in tours.items():
        t = []

        for node in tour[1:-1]:
            place = {
                "index": node,
                "name": pb.matrix.index[node],
            }
            if node < pb.n:
                place["type"] = "livraison"
                place["delivery"] = [
                    solver.Value(delivers[t][d, v, node]) for t in ["a", "f", "s"]
                ]
                place["palettes"] = solver.Value(palettes[d, v, node])
            else:
                place["type"] = "ramasse"
            t.append(place)
        sol["tours"][str(d) + ", " + str(v)] = t

    json.dump(sol, open("solutions/test2.json", "w"))


def reoptimize_tours(pb, solver, tours, delivers, obj):
    # Re-optimize tours
    for d in range(pb.n_days):
        for v in range(pb.m):
            if len(tours[d, v]) == 1:
                continue

            i = 1
            while i < len(tours[d, v]) - 1:
                t = tours[d, v][i]
                if t == 0 or t >= pb.n:
                    i += 1
                    continue
                # If the city is visited for nothing, skip it and adjust the objective value
                elif (
                    sum(solver.Value(delivers[i][d, v, t]) for i in ["a", "f", "s"])
                    == 0
                ):
                    obj = (
                        obj
                        - pb.matrix.iloc[tours[d, v][i - 1], t]
                        - pb.matrix.iloc[t, tours[d, v][i + 1]]
                        + pb.matrix.iloc[tours[d, v][i - 1], tours[d, v][i + 1]]
                    )
                    tours[d, v].pop(i)
                    print(f"\rRe-optimized tours to {obj/1000:.2f}km", end="")
                else:
                    i += 1
    print()
    return tours, int(obj)


def solve_vrp(
    pb,
    hint=None,
):
    print("Modeling...")

    model = cp_model.CpModel()

    pruned = prune_arcs(pb.n + pb.n_pdr, pb.matrix)
    print(f"    Pruned arcs : {100*len(pruned)/(pb.n+pb.n_pdr)**2:.2f}%")

    ########## Arcs & Visit variables ###########
    visits = {}
    arcs = {}
    for d in range(pb.n_days):
        for v in range(pb.m):
            # Identifies uniquely the vehicle-day pair (easier for arcs)
            vd = v + d * pb.m
            for node in range(pb.n + pb.n_pdr):
                visits[d, v, node] = model.NewBoolVar("visits_%i_%i_%i" % (d, v, node))

                # Add an arc from the node to itself (necessary for Circuit constraints)
                # If used, the node is not visited
                arcs[vd, node, node] = (node, node, visits[d, v, node].Not())

                for node2 in range(pb.n + pb.n_pdr):
                    if node == node2:
                        continue
                    if (node, node2) in pruned:
                        continue
                    # No arc from pdr to centre (except depot)
                    elif node >= pb.n and 0 < node2 and node2 < pb.n:
                        continue
                    else:
                        lit = model.NewBoolVar("arc_%i_%i_%i" % (vd, node, node2))
                        arcs[vd, node, node2] = (node, node2, lit)

            # If you don't visit the depot, visit nothing
            model.AddBoolAnd(
                visits[d, v, node].Not() for node in range(pb.n + pb.n_pdr)
            ).OnlyEnforceIf(visits[d, v, 0].Not())
            # for c in range(1, n):
            #     model.AddImplication(visits[v, 0].Not(), visits[v, c].Not())

            model.AddCircuit(
                [
                    arcs[vd, node1, node2]
                    for node1 in range(pb.n + pb.n_pdr)
                    for node2 in range(pb.n + pb.n_pdr)
                    if (node2 == 0 or node1 < pb.n or node2 >= pb.n)
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
    for d in range(pb.n_days):
        for v in range(pb.m):
            for c in range(pb.n):
                delivers["a"][d, v, c] = model.NewIntVar(
                    0,
                    min(pb.capacities[v], pb.demands["a"][c]),
                    f"delivers_a_{v}_{d}_{c}",
                )

                # Frais must be delivered with appropriate vehicles
                # if True:
                if pb.frais[v]:
                    delivers["f"][d, v, c] = model.NewIntVar(
                        0,
                        min(pb.capacities[v], pb.demands["f"][c]),
                        "delivers_f_%i_%i_%i" % (d, v, c),
                    )
                else:
                    delivers["f"][d, v, c] = 0

                delivers["s"][d, v, c] = model.NewIntVar(
                    0,
                    min(pb.capacities[v], pb.demands["s"][c]),
                    "delivers_s_%i_%i_%i" % (d, v, c),
                )

                palettes[d, v, c] = model.NewIntVar(
                    0, pb.sizes[v], "palettes_%i_%i_%i" % (d, v, c)
                )

                # Use an appropriate number of palettes
                model.Add(
                    palettes[d, v, c] * pb.max_palette_capacity
                    >= delivers["a"][d, v, c] + delivers["f"][d, v, c]
                )

                # Not visited -> delivery is 0
                model.Add(delivers["a"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["f"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )
                model.Add(delivers["s"][d, v, c] == 0).OnlyEnforceIf(
                    visits[d, v, c].Not()
                )

    for pdr in range(pb.n, pb.n + pb.n_pdr):
        # Visit pdrs enough times
        model.Add(
            sum(visits[d, v, pdr] for v in range(pb.m) for d in range(pb.n_days))
            >= pb.freqs_pdr[pdr - pb.n]
        )

        # But not twice a day
        for d in range(pb.n_days):
            model.AddAtMostOne(visits[d, v, pdr] for v in range(pb.m))

    # # Add hints
    # if hint:
    #     add_hints(model, arcs, visits, pb, hint)

    # Every customer must be served exactly his demand
    for c in range(1, pb.n):
        model.Add(
            sum(delivers["a"][d, v, c] for v in range(pb.m) for d in range(pb.n_days))
            == pb.demands["a"][c]
        )
        model.Add(
            sum(delivers["f"][d, v, c] for v in range(pb.m) for d in range(pb.n_days))
            == pb.demands["f"][c]
        )
        model.Add(
            sum(delivers["s"][d, v, c] for v in range(pb.m) for d in range(pb.n_days))
            == pb.demands["s"][c]
        )
        # model.Add(
        #     sum(palettes[d, v, c] for v in range(m) for d in range(pb.n_days))
        #     * max_palette_capacity
        #     >= pb.demands["a"][c] + pb.demands["f"][c]
        # )

    # add_domsym_breaking(model, n, m, n_days, delivers, visits)

    # Capacity constraints
    for d in range(pb.n_days):
        for v in range(pb.m):
            model.Add(
                sum(
                    [
                        delivers["a"][d, v, c]
                        + delivers["f"][d, v, c]
                        + delivers["s"][d, v, c]
                        for c in range(pb.n)
                    ]
                )
                <= pb.capacities[v]
            )

            model.Add(sum([palettes[d, v, c] for c in range(pb.n)]) <= pb.sizes[v])

    # Minimize total distance
    total_distance = sum(
        arc[2] * pb.matrix.iloc[arc[0], arc[1]]
        for arc in arcs.values()
        if arc[0] != arc[1]
    )
    model.Minimize(total_distance)

    solver = cp_model.CpSolver()

    solver.parameters.num_workers = 12
    # solver.parameters.max_time_in_seconds = 180
    # solver.parameters.min_num_lns_workers = 8

    print("Solving...")
    # solver.parameters.log_search_progress = True
    # status = solver.Solve(model)
    status = solver.Solve(model, LogPrinter())

    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        print("No solution found")
        exit()

    obj = solver.Value(total_distance)

    tours = read_tours(pb, arcs, solver)

    tours, obj = reoptimize_tours(pb, solver, tours, delivers, obj)

    write_sol(pb, solver, obj, tours, delivers, palettes)

    return (
        tours,
        obj,
    )
