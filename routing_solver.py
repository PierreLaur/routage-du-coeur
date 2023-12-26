from ortools.constraint_solver import pywrapcp, routing_enums_pb2
from math import inf
from time import time


def route_vrp(
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
    matrix = [
        [int(matrix.iloc[r, c]) for c in range(len(matrix.index))]
        for r in range(len(matrix.index))
    ]
    manager = pywrapcp.RoutingIndexManager(n + n_pdr, m * n_days, 0)
    model = pywrapcp.RoutingModel(manager)

    #################################################################

    def distance_callback(i, j):
        node_i = manager.IndexToNode(i)
        node_j = manager.IndexToNode(j)
        return matrix[node_i][node_j]

    transit_callback_index = model.RegisterTransitCallback(distance_callback)
    model.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    def demand_callback(i):
        node = manager.IndexToNode(i)
        if node >= n:
            return 0
        return demands["a"][node] + demands["f"][node]

    demand_callback_index = model.RegisterUnaryTransitCallback(demand_callback)
    model.AddDimensionWithVehicleCapacity(
        demand_callback_index,
        0,  # null capacity slack
        capacities * n_days,  # vehicle maximum capacities
        True,  # start cumul to zero
        "Capacity",
    )

    #################################################################

    def make_routing_monitor(routing_model: pywrapcp.RoutingModel) -> callable:
        class RoutingMonitor:
            def __init__(self, model: pywrapcp.RoutingModel):
                self.model = model
                self._best_objective = inf
                self.start_time = time()

            def __call__(self):
                curr = self.model.CostVar().Max()
                if curr < self._best_objective:
                    self._best_objective = min(self._best_objective, curr)
                    print(
                        f"[{time() - self.start_time :.2f}s] Current best : ",
                        self._best_objective / 1000,
                    )

        return RoutingMonitor(routing_model)

    def read_solution(manager, routing, solution):
        """Prints solution on console."""
        obj = solution.ObjectiveValue()
        tours = {}
        for vehicle_id in range(m * n_days):
            index = routing.Start(vehicle_id)
            route_distance = 0
            tour = []
            while not routing.IsEnd(index):
                tour.append(manager.IndexToNode(index))
                previous_index = index
                index = solution.Value(routing.NextVar(index))
                route_distance += routing.GetArcCostForVehicle(
                    previous_index, index, vehicle_id
                )
            tour.append(manager.IndexToNode(index))
            # plan_output += f"Distance of the route: {route_distance}m\n"

            d, v = divmod(vehicle_id, m)
            tours[v, d] = tour
        return obj, tours

    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PARALLEL_CHEAPEST_INSERTION
    )
    search_parameters.local_search_metaheuristic = (
        routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
    )
    search_parameters.time_limit.seconds = 5
    # search_parameters.log_search = True

    model.AddAtSolutionCallback(make_routing_monitor(model))

    # Solve the problem.
    # model.CloseModelWithParameters(search_parameters)
    solution = model.SolveWithParameters(search_parameters)

    # Print solution on console.
    if solution:
        obj, tours = read_solution(manager, model, solution)
        pass
    else:
        print("No solution found !")

    return obj, tours
