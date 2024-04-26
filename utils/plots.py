import pandas as pd
from utils.problem import Solution, Stop, StopType
import folium
import argparse
from dash import Dash, html
import dash_bootstrap_components as dbc
from datetime import date


def make_dashboard(sol: Solution):
    vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
    vehicle_names = vehicles.index

    def create_tour_card(key, tour: list[Stop]):
        card_items = []
        vehicle_name = vehicle_names[key[1]]

        for stop in tour:
            if stop.stop_type == StopType.Livraison:
                card_items.append(
                    html.Tr(
                        [
                            html.Td(
                                html.B(f"Livraison {stop.name}"),
                                style={
                                    "width": "50%",
                                    "background-color": "#f0f0f0",
                                },
                            ),
                            html.Td(
                                dbc.Table(
                                    [
                                        html.Tr(
                                            [
                                                html.Td(
                                                    "Palettes",
                                                    style={
                                                        "width": "31%",
                                                        "padding": "0 0px 0 5px",
                                                    },
                                                ),
                                                *[
                                                    html.Td(p, style={"width": "23%"})
                                                    for p in stop.palettes  # type: ignore
                                                ],
                                            ]
                                        ),
                                        html.Tr(
                                            [
                                                html.Td(
                                                    "Quantités",
                                                    style={
                                                        "width": "31%",
                                                        "padding": "0 0px 0 5px",
                                                    },
                                                ),
                                                *[
                                                    html.Td(q, style={"width": "23%"})
                                                    for q in stop.delivery  # type: ignore
                                                ],
                                            ]
                                        ),
                                    ],
                                    size="sm",
                                    style={"font-family": "monospace", "margin": "0"},
                                ),
                                style={
                                    "padding": "0",
                                    "background-color": "#f0f0f0",
                                },
                            ),
                        ],
                        style={
                            "margin": "0",
                            "padding": "0",
                        },
                    )
                )
            else:
                card_items.append(
                    html.Tr(
                        [
                            html.Td(
                                html.B(f"Ramasse {stop.name}"),
                                style={"width": "50%", "height": "100%"},
                            ),
                            html.Td("", style={"width": "50%", "height": "100%"}),
                        ],
                        style={"height": "100%"},
                    )
                )

        return dbc.Card(
            dbc.Row(
                [
                    dbc.Col(
                        html.B(
                            f"{vehicle_name}",
                            style={"margin": "0 0 0 3px", "padding": "0"},
                        ),
                        style={
                            "display": "flex",
                            "alignItems": "center",
                            "text-align": "center",
                            "max-width": "10%",
                            "min-height": "100%",
                        },
                    ),
                    dbc.Col(
                        dbc.Table(
                            html.Tbody(card_items, style={"min-height": "100%"}),
                            bordered=True,
                            size="sm",
                            style={
                                "margin": "0",
                                "overflow": "hidden",
                                "min-width": "90%",
                                "min-height": "100%",
                            },
                        ),
                    ),
                ],
            ),
            style={
                "margin": "0px 0% 3px 0%",
                "width": "100%",
                "border": "1px solid black",
                "overflow": "hidden",
            },
        )

    app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

    day_names = ["Lundi", "Mardi", "Mercredi", "Jeudi", "Vendredi"]
    tour_cards = {}
    for (d, v), tour in sol.tours.items():
        tour_cards[d, v] = create_tour_card((d, v), tour)

    days = [
        html.Div(
            [
                html.H5(
                    children=f"{day}",
                    style={"text-align": "center", "font-weight": "bold"},
                ),
                dbc.Row(
                    [
                        dbc.Col(),
                        dbc.Col(children=["A"], width=1, style={"margin-right": "7px"}),
                        dbc.Col(children=["F"], width=1, style={"margin-right": "5px"}),
                        dbc.Col(
                            children=["S"], width=2, style={"margin-right": "-7px"}
                        ),
                    ]
                ),
                *[tour_cards[d, v] for (d2, v) in tour_cards if d2 == d],
            ]
        )
        for d, day in enumerate(day_names)
    ]

    infos = dbc.Card(
        [
            html.H4(f"Date : {date.today()}"),
            html.H4(""),
            html.H4(f"Semaine : {sol.week}"),
            html.H4(
                f"Coûts totaux : {sol.total_costs:.0f}€     (Fixes {sol.fixed_costs:.0f} | Variables {sol.variable_costs:.0f})"
            ),
            html.H4(f"Distance totale : {sol.total_distance/1000:.0f}km"),
            html.H4(f"Véhicules : {sum(sol.vehicles_used)}"),
        ],
        outline=True,
        style={"border": "none", "margin-top": "20%"},
    )

    app.layout = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [days[0], days[3]],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding-left": "15px",
                            "height": "100%",
                        },
                    ),
                    dbc.Col(
                        [days[1], days[4]],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px",
                            "height": "100%",
                        },
                    ),
                    dbc.Col(
                        [days[2], infos],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding-right": "15px",
                            "height": "100%",
                        },
                    ),
                ],
                style={
                    "height": "100%",
                    "width": "100%",
                    "padding": "0",
                    "margin": "0",
                },
            )
        ],
        style={
            "height": "1050px",
            "width": "1485px",
            "margin": "0",
            "padding": "0",
            "fontSize": 11,
            "border": "1px solid black",
        },
    )
    app.run(debug=False)


def print_to_txt(sol: Solution, output_file_path, display=False):
    """Generates a pretty printed txt version of the solution and prints it to a specified file"""

    vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
    vehicle_names = vehicles.index

    jours_map = {0: "Lundi", 1: "Mardi", 2: "Mercredi", 3: "Jeudi", 4: "Vendredi"}

    output = ""
    output += f"- - - - - - TOURNEES SEMAINE {sol.week.value} - - - - - -\n"

    n_days = len(set(k[0] for k in sol.tours))

    for d in range(n_days):
        if not any(d2 == d for (d2, v) in sol.tours):
            continue

        output += f"\n- - - - {jours_map[d].upper()} - - - -\n\n"
        day_tours = {v: tour for (d2, v), tour in sol.tours.items() if d2 == d}
        for v, tour in day_tours.items():
            est_time = sol.tour_durations[d, v] / 60
            est_time = f"{est_time//60:.0f}h{est_time%60:.0f}m"
            output += f"\tVéhicule {v} ({vehicle_names[v]}) - est. {est_time}\n"

            for stop in tour:
                product_types = ""
                palettes = ""
                name = stop.name
                if stop.index == 0:
                    name = "    [Retour Dépôt]"
                elif stop.stop_type == StopType.Ramasse:
                    product_types = "Ramasse"
                else:
                    if stop.stop_type == StopType.Liv_Ramasse:
                        product_types = "Liv. ramasse"

                    for i, product_type in enumerate("AFS"):
                        deliv = stop.delivery[i]
                        if deliv > 0:
                            product_types += f"{f'{deliv}'+product_type:5}"
                        else:
                            product_types += f"{'':5}"

                    pals = stop.palettes
                    palettes = f"{pals[0]}PA" if pals[0] else "   "
                    palettes += f" {pals[1]:.1f}PF" if pals[1] else "      "
                    palettes += f" {pals[2]:.1f}PS" if pals[2] else "      "
                    palettes += (
                        f" +{stop.norvegiennes} norvégienne{'s' if stop.norvegiennes > 1 else ''}"
                        if stop.norvegiennes > 0
                        else ""
                    )

                output += f"\t\t{name:40}\t{product_types}\t{palettes}\n"

    output += f"\nCoûts totaux : {sol.total_costs:.0f}€     (Fixes {sol.fixed_costs:.0f}€ | Variables {sol.variable_costs:.0f}€)"
    output += f"\nDistance totale : {round(sol.total_distance/1000):d}km"

    vehicles_used = vehicles.index.where(sol.vehicles_used)
    output += f"\nVéhicules utilisés : {' - '.join(f'{v} {k}' for k, v in vehicles_used.value_counts().items())}"

    output += f"\n\nDate : {date.today()}"

    with open(output_file_path, "w") as txt_file:
        txt_file.write(output)

    if display:
        print(output)
    print(f"Wrote solution to {output_file_path}")


def plot_solution(sol: Solution, output_file_path=None, day=None, silent=False):
    """Makes a folium plot of the solution
    If output_file_path, writes it to the specified file
    If day, only plots the tours for the specified day"""

    centres = pd.read_excel("data/centres.xlsx")
    vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
    pdr = pd.read_excel("data/points_de_ramasse.xlsx")

    n_centres = len(centres.index)
    n_pdr = len(pdr.index)
    m = len(vehicles.index)
    n_days = 5

    jours_map = {0: "Lundi", 1: "Mardi", 2: "Mercredi", 3: "Jeudi", 4: "Vendredi"}
    colors = [
        "red",
        "orange",
        "grey",
        "blue",
        "darkblue",
        "black",
        "purple",
        "green",
        "darkgreen",
        "white",
        "white",
        "white",
    ]

    lats = centres["Latitude"].tolist() + pdr["Latitude"].tolist()
    longs = centres["Longitude"].tolist() + pdr["Longitude"].tolist()
    coords_depot = (lats[0], longs[0])

    folmap = folium.Map(
        tiles="CartoDB Positron",
        location=(lats[0], longs[0]),
        max_bounds=True,
        min_lat=min(lats) - 0.5,
        max_lat=max(lats) + 0.5,
        min_lon=min(longs) - 1,
        max_lon=max(longs) + 1,
    )

    # Add the centres and points de ramasse as markers
    for i in range(n_centres):
        folium.CircleMarker(
            location=(lats[i], longs[i]),
            tooltip=centres["Nom"][i],
            color="green" if i == 0 else "red",
        ).add_to(folmap)

    for i in range(n_pdr):
        folium.CircleMarker(
            location=(lats[n_centres + i], longs[n_centres + i]),
            tooltip=pdr["Nom"][i],
            color="blue",
        ).add_to(folmap)

    # Add the tours
    days_to_plot = [day] if day is not None else range(n_days)

    group = {d: folium.FeatureGroup(jours_map[d]).add_to(folmap) for d in days_to_plot}

    for (d, v), tour in sol.tours.items():
        if d not in days_to_plot:
            continue

        tour_coords = []
        for stop in tour:
            tour_coords.append((lats[stop.index], longs[stop.index]))
        tour_coords = [coords_depot] + tour_coords + [coords_depot]

        folium.PolyLine(
            tour_coords,
            tooltip=f"Véhicule {v} ({vehicles.index[v]})",
            color=colors[v],
        ).add_to(group[d])

    folium.LayerControl().add_to(folmap)

    if output_file_path:
        folmap.save(output_file_path)
        if not silent:
            print("Map saved to ", output_file_path)

    return m


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="json solution file")
    parser.add_argument(
        "--map",
        "-m",
        action="store_true",
        help="make a html map",
    )
    parser.add_argument(
        "--txt",
        "-t",
        action="store_true",
        help="make a txt file",
    )
    parser.add_argument(
        "--dashboard",
        "-d",
        action="store_true",
        help="make a dashboard",
    )

    args = parser.parse_args()

    file_name = args.infile.split(".")[0]
    sol = Solution.from_json(args.infile)

    if args.txt:
        print_to_txt(sol, file_name + ".txt")

    if args.map:
        plot_solution(sol, file_name + ".html")

    if args.dashboard:
        make_dashboard(sol)
