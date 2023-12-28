import plotly.express as px
import pandas as pd
import plotly.graph_objects as go


def plot_tours(tours):
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
