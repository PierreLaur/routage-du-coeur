from ortools.sat.python import cp_model
import pandas as pd
from math import inf


class LogPrinter(cp_model.CpSolverSolutionCallback):
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


def solve_vrp(centres, vehicles, matrix):
    n = len(centres.index)

    m = vehicles["Nombre"].sum()
    cap = vehicles["Capacité (kg)"].max()
    demands = centres["Tonnage Ambiant (kg)"].fillna(0).astype(int).tolist()

    model = cp_model.CpModel()

    visits = {}
    delivers = {}
    arcs = {}
    for v in range(m):
        arcs[v] = []

        for i in range(n):
            visits[v, i] = model.NewBoolVar("visits_%i_%i" % (v, i))
            delivers[v, i] = model.NewIntVar(0, cap, "delivers_%i_%i" % (v, i))

            # Add an arc from the node to itself (necessary for Circuit constraints)
            arcs[v].append((i, i, visits[v, i].Not()))

            for j in range(n):
                if i == j:
                    continue
                lit = model.NewBoolVar("arc_%i_%i_%i" % (v, i, j))
                arcs[v].append((i, j, lit))

            # Not visited -> delivery is 0
            model.Add(delivers[v, i] == 0).OnlyEnforceIf(visits[v, i].Not())

        model.AddCircuit(arcs[v])

    # Every customer must be served exactly his demand
    for c in range(1, n):
        model.Add(sum(delivers[v, c] for v in range(m)) == demands[c])

    # Capacity constraints
    for v in range(m):
        model.Add(sum([delivers[v, c] for c in range(n)]) <= cap)

    # Minimize total distance
    total_distance = sum(
        arc[2] * matrix.iloc[arc[0], arc[1]] for arcsv in arcs.values() for arc in arcsv
    )
    model.Minimize(total_distance)

    solver = cp_model.CpSolver()

    solver.parameters.num_workers = 16
    # solver.parameters.log_search_progress = True
    solver.parameters.max_time_in_seconds = 120

    status = solver.Solve(model, LogPrinter())

    tours = []
    deliveries = []
    obj = solver.Value(total_distance) / 1000
    for v in range(m):
        tour = [0]
        delivery = [0]
        for arc in arcs[v]:
            if arc[0] == arc[1]:
                continue

            if solver.Value(arc[2]):
                tour.append(arc[1])
                delivery.append(solver.Value(delivers[v, arc[1]]))
        tours.append(tour)
        deliveries.append(delivery)

    return tours, obj, deliveries


if __name__ == "__main__":
    centres = pd.read_excel("data/centres.xlsx", index_col=0)
    vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
    matrix = pd.read_excel("data/matrix.xlsx", index_col=0)

    tours, obj, deliveries = solve_vrp(centres, vehicles, matrix)

    for v in range(len(tours)):
        print(f"Vehicule {v} tour : ")
        print("\t", tours[v])
        print("\t", deliveries[v])
