import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("route", help="route")
    parser.add_argument("-t", help="with traffic", action="store_true")

    args = parser.parse_args()
    path = (
        "data/matrices/traffic_duration_matrix.xlsx"
        if args.t
        else "data/matrices/no_traffic_duration_matrix.xlsx"
    )

    matrix = pd.read_excel(path, index_col=0)

    route = args.route.split(",")
    if not route[0].isnumeric():
        route = (
            [0]
            + [matrix.index.tolist().index(route[i]) for i in range(len(route))]
            + [0]
        )
    else:
        route = [0] + [int(x) for x in args.route.split()] + [0]

    matrix = matrix.values

    time_seconds = 0
    last_ramasse = 0
    for i in range(len(route) - 1):
        print(route[i], route[i + 1], matrix[route[i], route[i + 1]])
        time_seconds += matrix[route[i], route[i + 1]]

        if route[i + 1] == 0:
            break

        if route[i + 1] < 29:
            time_seconds += 30 * 60
        else:
            last_ramasse = time_seconds
            time_seconds += 30 * 60

    time = pd.Timedelta(seconds=time_seconds)
    print(f"Time : {time}")
    print("Last ramasse : ", pd.Timedelta(seconds=last_ramasse))
