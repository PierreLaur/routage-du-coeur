import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import sys
import json
from problem import Problem, read_problem


def pretty_print_solution(file, week):
    """Generates a pretty printed version of the solution and returns it as a string"""
    with open(file) as f:
        sol = json.load(f)

        centres = pd.read_excel("data/centres.xlsx")
        vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)

        pb = read_problem(
            "data/centres.xlsx",
            "data/points_de_ramasse.xlsx",
            "data/vehicules.xlsx",
            "data/euclidean_matrix.xlsx",
            week,
        )

        obj = sol["total_distance"]
        tours_strkey = sol["tours"]
        tours = {}

        for key, tour in tours_strkey.items():
            d, v = map(int, key.split(", "))
            tours[d, v] = tour

        jours_map = {0: "Lundi", 1: "Mardi", 2: "Mercredi", 3: "Jeudi", 4: "Vendredi"}

        output = ""
        output += f"- - - - - - TOURNEES SEMAINE {week} - - - - - -\n"
        for d in range(pb.n_days):
            if not any((d, v) in tours for v in range(pb.m)):
                continue

            output += f"\n- - - - {jours_map[d].upper()} - - - -\n\n"
            for v in range(pb.m):
                if not (d, v) in tours:
                    continue

                output += f"\tVéhicule {v} ({vehicles.index[v]}) \n"
                tour = tours[d, v]

                for place in tour:
                    product_types = ""
                    palettes = ""
                    if place["type"] == "livraison":
                        for i, product_type in enumerate("AFS"):
                            if place["delivery"][i] > 0:
                                product_types += product_type

                        palettes = f"{place['palettes']} palette{'s' if place['palettes'] > 1 else ''}"

                    else:
                        product_types = "Ramasse"

                    output += f"\t\t{place['name']:30}\t{product_types}\t{palettes}\n"

        output += f"\nDistance totale : {round(obj/1000):d}km"

        return output


def plot_solution(file):
    # TODO
    coords_centres = pd.read_excel("data/centres.xlsx")
    coords_pdr = pd.read_excel("data/points_de_ramasse.xlsx")

    # lats = coords_centres["Latitude"].tolist() + coords_pdr["Latitude"].tolist()
    # longs = coords_centres["Longitude"].tolist() + coords_pdr["Longitude"].tolist()
    # colors = [
    #     "red",
    #     "blue",
    #     "green",
    #     "purple",
    #     "orange",
    #     "pink",
    #     "cyan",
    #     "brown",
    #     "gray",
    #     "olive",
    # ]

    # fig = px.scatter(
    #     y=lats,
    #     x=longs,
    #     text=coords_centres["Nom"].tolist() + coords_pdr["Nom"].tolist(),
    #     title="Restos du Cœur",
    # )

    # for i, tour in enumerate(tours):
    #     color = colors[i % len(colors)]

    #     fig.add_traces(
    #         go.Scatter(
    #             y=[lats[t] for t in tour],
    #             x=[longs[t] for t in tour],
    #             mode="lines+markers",
    #             line=dict(color=color),
    #             name=f"Tour {i}",
    #         )
    #     )
    # fig.show()


if __name__ == "__main__":
    output = pretty_print_solution(sys.argv[1], week=int(sys.argv[2]))

    with open(f"solutions/week_{sys.argv[2]}.txt", "w") as f:
        f.write(output)
    print(output)
