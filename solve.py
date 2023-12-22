from ortools.sat.python import cp_model
import pandas as pd
from math import inf
from plots import plot_tours


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


def check_solution(
    centres, vehicles, tours, obj, deliveries, vehicle_types, visits, arcs
):
    n = len(centres.index)
    m = vehicles["Nombre"].sum() * 10
    n_vehicle_types = len(vehicles.index)

    demands_a = centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist()
    demands_f = centres["Tonnage Frais (kg)"].fillna(0).astype(int).tolist()
    demands_s = centres["Tonnage Surgelé (kg)"].fillna(0).astype(int).tolist()

    capacities = vehicles["Capacité (kg)"].astype(int).tolist()
    sizes = vehicles["Taille (Palettes)"].astype(int).tolist()
    print(sizes)
    exit()

    da, df, ds = deliveries

    for v in range(m):
        for c in range(n):
            if visits[v, c]:
                assert sum([arcs[v, c2, c] for c2 in range(n) if c2 != c]) == 1
                assert arcs[v, c, c] == 0
                if da[v, c] + df[v, c] + ds[v, c] == 0:
                    f"Warning : vehicle {v} visits centre {c} with no deliveries"
            else:
                assert sum([arcs[v, c2, c] for c2 in range(n) if c2 != c]) == 0
                assert arcs[v, c, c] == 1
                assert da[v, c] + df[v, c] + ds[v, c] == 0

        for i, cap in enumerate(capacities):
            if vehicle_types[v] == i:
                assert sum([da[v, c] + df[v, c] for c in range(n)]) <= cap

        if vehicle_types[v] != 2:
            assert sum([df[v, c] for c in range(n)]) == 0

        # TODO : check vehicle types

    for v in range(m):
        for c1 in range(n):
            for c2 in range(n):
                if c1 == c2:
                    continue

    for c in range(1, n):
        assert sum([da[v, c] for v in range(m)]) == demands_a[c]
        assert sum([df[v, c] for v in range(m)]) == demands_f[c]
        assert sum([ds[v, c] for v in range(m)]) == demands_s[c]

    print("Assertions ok")


def add_domsym_breaking(
    model,
    n,
    m,
    delivers_a,
    delivers_f,
    delivers_s,
    visits,
    n_vehicle_types,
    vehicle_is_of_type,
):
    # # Dominance breaking
    # for v in range(m):
    #     for c in range(1, n):
    #         # Visited -> delivery is not 0
    #         model.Add(
    #             delivers_a[v, c] + delivers_f[v, c] + delivers_s[v, c] > 0
    #         ).OnlyEnforceIf(visits[v, c])

    # # Symmetry breaking
    # for v in range(1, m):
    #     model.AddImplication(visits[v, 0], visits[v - 1, 0])

    for v in range(1, m):
        for i in range(n_vehicle_types - 1):
            for j in range(i + 1, n_vehicle_types):
                model.AddImplication(
                    vehicle_is_of_type[v, i], vehicle_is_of_type[v - 1, j].Not()
                )


def solve_vrp(centres, vehicles, matrix):
    n = len(centres.index)
    m = vehicles["Nombre"].sum() * 10
    n_vehicle_types = len(vehicles.index)

    capacities = vehicles["Capacité (kg)"].astype(int).tolist()

    demands_a = centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist()
    demands_f = centres["Tonnage Frais (kg)"].fillna(0).astype(int).tolist()
    demands_s = centres["Tonnage Surgelé (kg)"].fillna(0).astype(int).tolist()

    model = cp_model.CpModel()

    visits = {}
    delivers_a = {}
    delivers_f = {}
    delivers_s = {}
    arcs = {}
    for v in range(m):
        for i in range(n):
            visits[v, i] = model.NewBoolVar("visits_%i_%i" % (v, i))

            delivers_a[v, i] = model.NewIntVar(
                0, min(max(capacities), max(demands_a)), "delivers_a_%i_%i" % (v, i)
            )
            delivers_f[v, i] = model.NewIntVar(
                0, min(max(capacities), max(demands_f)), "delivers_f_%i_%i" % (v, i)
            )
            delivers_s[v, i] = model.NewIntVar(
                0, min(max(capacities), max(demands_s)), "delivers_s_%i_%i" % (v, i)
            )

            # Add an arc from the node to itself (necessary for Circuit constraints)
            # If used, the node is not visited
            arcs[v, i, i] = (i, i, visits[v, i].Not())

            for j in range(n):
                if i == j:
                    continue
                else:
                    lit = model.NewBoolVar("arc_%i_%i_%i" % (v, i, j))
                    arcs[v, i, j] = (i, j, lit)

            # Not visited -> delivery is 0
            model.Add(delivers_a[v, i] == 0).OnlyEnforceIf(visits[v, i].Not())
            model.Add(delivers_f[v, i] == 0).OnlyEnforceIf(visits[v, i].Not())
            model.Add(delivers_s[v, i] == 0).OnlyEnforceIf(visits[v, i].Not())

        # If you don't visit the depot, visit nothing
        model.AddBoolAnd(visits[v, c].Not() for c in range(n)).OnlyEnforceIf(
            visits[v, 0].Not()
        )
        # for c in range(1, n):
        #     model.AddImplication(visits[v, 0].Not(), visits[v, c].Not())

        model.AddCircuit([arcs[v, i, j] for i in range(n) for j in range(n)])

    vehicle_is_of_type = {
        (v, i): model.NewBoolVar("vehicle_type_%i_%i" % (v, i))
        for v in range(m)
        for i in range(n_vehicle_types)
    }

    # add_domsym_breaking(
    #     model,
    #     n,
    #     m,
    #     delivers_a,
    #     delivers_f,
    #     delivers_s,
    #     visits,
    #     n_vehicle_types,
    #     vehicle_is_of_type,
    # )

    for v in range(m):
        model.AddExactlyOne(vehicle_is_of_type[v, i] for i in range(n_vehicle_types))

    # Every customer must be served exactly his demand
    for c in range(1, n):
        model.Add(sum(delivers_a[v, c] for v in range(m)) == demands_a[c])
        model.Add(sum(delivers_f[v, c] for v in range(m)) == demands_f[c])
        model.Add(sum(delivers_s[v, c] for v in range(m)) == demands_s[c])

    # Capacity constraints
    for v in range(m):
        for i, cap in enumerate(capacities):
            model.Add(
                sum([delivers_a[v, c] + delivers_f[v, c] for c in range(n)]) <= cap
            ).OnlyEnforceIf(vehicle_is_of_type[v, i])
        # model.Add(
        #     sum([delivers_a[v, c] + delivers_f[v, c] for c in range(n)])
        #     <= sum(cap * vehicle_is_of_type[v, i] for i, cap in enumerate(capacities))
        # )

    # Frais must be served with type 2 vehicle (camion frigo)
    for v in range(m):
        model.Add(sum([delivers_f[v, c] for c in range(n)]) == 0).OnlyEnforceIf(
            vehicle_is_of_type[v, 2].Not()
        )

        # for c in range(n):
        #     model.Add(delivers_f[v, c] == 0).OnlyEnforceIf(
        #         vehicle_is_of_type[v, 2].Not()
        #     )

    # Minimize total distance
    total_distance = sum(
        arc[2] * matrix.iloc[arc[0], arc[1]]
        for arc in arcs.values()
        if arc[0] != arc[1]
    )
    model.Minimize(total_distance)

    solver = cp_model.CpSolver()

    solver.parameters.num_workers = 16
    # solver.parameters.max_time_in_seconds = 180
    # solver.parameters.min_num_lns_workers = 12

    solver.parameters.log_search_progress = True
    # status = solver.Solve(model, LogPrinter())
    status = solver.Solve(model)

    tours = []
    vehicle_types = []
    obj = solver.Value(total_distance) / 1000
    for v in range(m):
        tour = [0]

        while True:
            a = tour[-1]
            found = False
            for b in range(n):
                if a == b:
                    continue

                arc = arcs[v, a, b]
                # print(v, "\t", arc, a, b, solver.Value(arc[2]))
                if solver.Value(arc[2]):
                    tour.append(b)
                    found = True
                    break

            if not found:
                break
            if tour[-1] == 0:
                break

        if len(tour) > 1:
            tours.append(tour)

        for i in range(n_vehicle_types):
            if solver.Value(vehicle_is_of_type[v, i]):
                vehicle_types.append(i)

    deliveries = (
        {k: solver.Value(v) for k, v in delivers_a.items()},
        {k: solver.Value(v) for k, v in delivers_f.items()},
        {k: solver.Value(v) for k, v in delivers_s.items()},
    )

    return (
        tours,
        obj,
        deliveries,
        vehicle_types,
        {k: solver.Value(v) for k, v in visits.items()},
        {k: solver.Value(v[2]) for k, v in arcs.items()},
    )


if __name__ == "__main__":
    centres = pd.read_excel("data/centres.xlsx", index_col=0)
    vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
    matrix = pd.read_excel("data/matrix.xlsx", index_col=0)

    tours, obj, deliveries, vehicle_types, visits, arcs = solve_vrp(
        centres, vehicles, matrix
    )

    check_solution(
        centres, vehicles, tours, obj, deliveries, vehicle_types, visits, arcs
    )

    for v in range(len(tours)):
        print(f"Vehicule {v} ({vehicles.index[vehicle_types[v]]}) tour : ")
        print("\t", tours[v])

    plot_tours(tours)
