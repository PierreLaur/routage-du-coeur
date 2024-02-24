from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
from utils.problem import Problem
from time import time
from math import inf, ceil
from functools import cache


def make_routing_monitor(pb, routing_model: pywrapcp.RoutingModel, scaling_factor):
    class RoutingMonitor:
        def __init__(self, pb: Problem, model: pywrapcp.RoutingModel, scaling_factor):
            self.pb = pb
            self.model = model
            self._best_objective = inf
            self.start_time = time()
            self.scaling_factor = scaling_factor
            self.fixed_cost = sum(pb.vehicle_allowed) * pb.weekly_fixed_cost

        def __call__(self):
            curr = self.model.CostVar().Max()
            if curr < self._best_objective:
                self._best_objective = min(self._best_objective, curr)
                print(
                    f"[{time() - self.start_time :.2f}s] Current best : {self.fixed_cost + self._best_objective / (self.scaling_factor * 100000):.2f}E"
                )

    return RoutingMonitor(pb, routing_model, scaling_factor)


def solve(pb: Problem, unit_delivery_size=100):

    ### Preprocessing
    scaling_factor = 10

    # Duplicate customers as many times as needed for split delivery
    real_index = {0: 0}
    node_index = {c: [] for c in range(pb.n + pb.n_pdr)}
    node_type = {}
    i = 1
    for c in range(pb.n):
        for type in ["a", "f", "s"]:
            n_unit_deliveries = ceil(pb.demands[type][c] / unit_delivery_size)
            for unit in range(n_unit_deliveries):
                real_index[i] = c
                node_index[c] += [i]
                node_type[i] = type
                i += 1
    for p in range(pb.n, pb.n + pb.n_pdr):
        for d in pb.j_de_ramasse[p - pb.n]:
            real_index[i] = p
            node_index[p] += [i]
            node_type[i] = "r"
            i += 1

    n_nodes = i
    dist_matrix = [
        [pb.distance_matrix.iloc[real_index[i], real_index[j]] for j in range(n_nodes)]
        for i in range(n_nodes)
    ]

    allowed_vehicles = [v for v in range(pb.m) if pb.vehicle_allowed[v]]
    n_vehicles = len(allowed_vehicles)

    vehicle_days = range(n_vehicles * pb.n_days)
    vehicle_day = lambda d, v: v + d * n_vehicles

    # Allowed vehicle-days for each node
    allowed_vehicle_days_for_node = {i: set() for i in range(n_nodes)}

    for c in range(1, pb.n):
        for j in pb.j_de_livraison_possibles[c]:
            for node in node_index[c]:
                for v_index, v in enumerate(allowed_vehicles):
                    if v == 0 and node_type[node] != "a":
                        continue
                    allowed_vehicle_days_for_node[node].add(vehicle_day(j, v_index))
    for p in range(pb.n_pdr):
        days = pb.j_de_ramasse[p]
        for i, node in enumerate(node_index[p + pb.n]):
            for v in range(n_vehicles):
                allowed_vehicle_days_for_node[node].add(vehicle_day(days[i], v))

    def print_solution(pb, manager, routing, solution):
        """Prints solution on console."""
        print(f"Objective: {solution.ObjectiveValue() * 0.00001 / scaling_factor:.2f}E")
        for d in range(pb.n_days):
            print(f"- - - - DAY {d} - - - ")
            for v in range(n_vehicles):
                vehicle_id = vehicle_day(d, v)
                index = routing.Start(vehicle_id)
                plan_output = f"\tRoute for vehicle {v}:\n\t"
                previous_index = None
                while not routing.IsEnd(index):
                    if not previous_index or (
                        real_index[manager.IndexToNode(index)]
                        != real_index[manager.IndexToNode(previous_index)]
                    ):
                        plan_output += f" {real_index[manager.IndexToNode(index)]} -> "
                    previous_index = index
                    index = solution.Value(routing.NextVar(index))
                plan_output += f"{real_index[manager.IndexToNode(index)]}"
                print(plan_output)

    ##########################################################################################################

    manager = pywrapcp.RoutingIndexManager(n_nodes, n_vehicles * pb.n_days, 0)
    routing = pywrapcp.RoutingModel(manager)

    # Disallow going from pdrs to centres
    for i in range(n_nodes):
        for j in range(n_nodes):
            if real_index[i] >= pb.n and real_index[j] != 0 and real_index[j] < pb.n:
                i_index = manager.NodeToIndex(i)
                j_index = manager.NodeToIndex(j)
                routing.NextVar(i_index).RemoveValue(j_index)
                # dist_matrix[i][j] = 999999999

    ### Objectives ###
    arc_costs = [
        cache(
            lambda i, j: dist_matrix[manager.IndexToNode(i)][manager.IndexToNode(j)]
            * pb.consumptions[v]
            * round(scaling_factor * pb.fuel_cost)
        )
        for d in range(pb.n_days)
        for v in allowed_vehicles
    ]
    for vd in vehicle_days:
        transit_callback_index = routing.RegisterTransitCallback(arc_costs[vd])
        routing.SetArcCostEvaluatorOfVehicle(transit_callback_index, vd)
    # routing.SetFixedCostOfAllVehicles(pb.weekly_fixed_cost)

    ### Capacity Constraints ###
    @cache
    def demand_callback(i):
        node_i = manager.IndexToNode(i)
        return (
            0
            if node_i == 0
            else (
                pb.weights[real_index[node_i] - pb.n]
                if node_type[node_i] == "r"
                else unit_delivery_size
            )
        )

    routing.AddDimensionWithVehicleCapacity(
        routing.RegisterUnaryTransitCallback(demand_callback),
        0,  # null capacity slack
        [
            pb.capacities[v] for d in range(pb.n_days) for v in allowed_vehicles
        ],  # vehicle maximum capacities
        True,  # start cumul to zero
        "Capacity",
    )

    for node in range(1, n_nodes):
        routing.SetAllowedVehiclesForIndex(
            allowed_vehicle_days_for_node[node], manager.NodeToIndex(node)
        )

    ### Parameters ###
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PARALLEL_CHEAPEST_INSERTION
    )
    search_parameters.local_search_metaheuristic = (
        routing_enums_pb2.LocalSearchMetaheuristic.SIMULATED_ANNEALING
    )
    search_parameters.time_limit.seconds = 10
    search_parameters.use_cp_sat = True

    ##########################################################################################################

    routing.AddAtSolutionCallback(make_routing_monitor(pb, routing, scaling_factor))
    solution = routing.SolveWithParameters(search_parameters)

    if solution:
        pass
        print_solution(pb, manager, routing, solution)
    else:
        print("No solution found !")

    return [], 0
