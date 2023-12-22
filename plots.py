import plotly.express as px
import pandas as pd
import plotly.graph_objects as go


def plot_tours(tours):
    coords = pd.read_excel("data/centres.xlsx")
    lats = coords["Latitude"]
    longs = coords["Longitude"]
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
        text=coords["Centre"],
        title="Restos du Cœur",
    )

    for i, tour in enumerate(tours):
        color = colors[i % len(colors)]

        fig.add_traces(
            go.Scatter(
                y=lats[tour],
                x=longs[tour],
                mode="lines+markers",
                line=dict(color=color),
                name=f"Tour {i}",
            )
        )
    fig.show()
