import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import sys
import json


def read_solution(file):
    with open(file) as f:
        sol = json.load(f)

        obj = sol["total_distance"]

        tours = sol["tours"]


def pretty_print_solution(file):
    with open(file) as f:
        sol = json.load(f)

        obj = sol["total_distance"]
        tours = sol["tours"]

        print("")

    # TODO
    # for d in range(n_days):
    #     print(f"-- JOUR {d} --")
    #     for v in range(m):
    #         if len(tours[d, v]) == 1:
    #             continue
    #         print(f"Vehicule {v} ({vehicles.index[v]}) tour : \n", end="")

    #         for t in tours[d, v]:
    #             delivery = None
    #             pals = 0
    #             name = matrix.index[t]
    #             if t == 0 or t >= n:
    #                 continue
    #             else:
    #                 delivery = tuple(deliveries[i][d, v, t] for i in range(3))
    #                 pals = palettes[d, v, t]
    #             print(
    #                 f"\t{name:20} {str(delivery) if delivery else ''} {' - ' + str(pals)+' palettes'}"
    #             )
    pass


def plot_solution(file):
    # TODO
    coords_centres = pd.read_excel("data/centres.xlsx")
    coords_pdr = pd.read_excel("data/points_de_ramasse.xlsx")

    lats = coords_centres["Latitude"].tolist() + coords_pdr["Latitude"].tolist()
    longs = coords_centres["Longitude"].tolist() + coords_pdr["Longitude"].tolist()
    colors = [
        "red",
        "blue",
        "green",
        "purple",
        "orange",
        "pink",
        "cyan",
        "brown",
        "gray",
        "olive",
    ]

    fig = px.scatter(
        y=lats,
        x=longs,
        text=coords_centres["Nom"].tolist() + coords_pdr["Nom"].tolist(),
        title="Restos du Cœur",
    )

    for i, tour in enumerate(tours):
        color = colors[i % len(colors)]

        fig.add_traces(
            go.Scatter(
                y=[lats[t] for t in tour],
                x=[longs[t] for t in tour],
                mode="lines+markers",
                line=dict(color=color),
                name=f"Tour {i}",
            )
        )
    fig.show()


if __name__ == "__main__":
    read_solution(sys.argv[1])
