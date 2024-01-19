import pandas as pd
import sys
import json
from problem import Problem, read_problem
import folium
import yaml
import argparse
from dash import Dash, html, dcc
import dash_bootstrap_components as dbc


def make_dashboard(file, week, output_file):
    def create_tour_card(name, tour):
        card_items = []

        for item in tour:
            if isinstance(item, dict):
                tour_key, tour_info = list(item.items())[0]
                palettes, quantities = tour_info["Palettes"], tour_info["Quantites"]
                card_items.append(
                    html.Tr(
                        [
                            html.Td(
                                html.B(f"{tour_key}"),
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
                                                    f"Palettes",
                                                    style={
                                                        "width": "31%",
                                                        "padding": "0 0px 0 5px",
                                                    },
                                                ),
                                                *[
                                                    html.Td(p, style={"width": "23%"})
                                                    for p in palettes.split()
                                                ],
                                            ]
                                        ),
                                        html.Tr(
                                            [
                                                html.Td(
                                                    f"Quantités",
                                                    style={
                                                        "width": "31%",
                                                        "padding": "0 0px 0 5px",
                                                    },
                                                ),
                                                *[
                                                    html.Td(q, style={"width": "23%"})
                                                    for q in quantities.split()
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
                                html.B(item), style={"width": "50%", "height": "100%"}
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
                            f"{name}", style={"margin": "0 0 0 3px", "padding": "0"}
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
                "margin": "0px 1% 5px 1%",
                "width": "98%",
                "border": "1px solid black",
                "overflow": "hidden",
            },
        )

    def create_map(day):
        return (
            html.Iframe(
                srcDoc=open(f"solutions/test_{day}.html", "r").read(),
                style={
                    "height": "100%",
                    "width": "100%",
                    "border": "1px solid black",
                },
            ),
        )

    app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

    data = yaml.safe_load(open(file))[f"Tournees Semaine {week}"]

    days = []
    day_names = ["Lundi", "Mardi", "Mercredi", "Jeudi", "Vendredi"]
    for day in day_names:
        vehicules = []
        for k, v in data["Tournees"][day].items():
            card = create_tour_card(k, v)
            vehicules.append(card)
        days.append(
            [
                html.H5(
                    children=f"{day}",
                    style={"text-align": "center", "font-weight": "bold"},
                ),
                *vehicules,
            ]
        )

    plot_solution("solutions/week_1.json", 1, "solutions/test_2.html", 1)

    app.layout = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        days[0],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px"
                        },
                    ),
                    dbc.Col(
                        days[1],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px"
                        },
                    ),
                    dbc.Col(
                        days[2],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px"
                        },
                    ),
                ],
                className="grid-container",
                style={"height": "30%", "margin": "0px"},
            ),
            dbc.Row(
                [
                    dbc.Col(
                        create_map(1),
                        width=4,
                        style={"margin": "0px", "padding": "0px"},
                    ),
                    dbc.Col(
                        create_map(2),
                        width=4,
                        style={"margin": "0px", "padding": "0px"},
                    ),
                    dbc.Col(
                        create_map(1),
                        width=4,
                        style={"margin": "0px", "padding": "0px"},
                    ),
                ],
                style={"height": "20%", "margin": "0px"},
            ),
            dbc.Row(
                [
                    dbc.Col(
                        days[3],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px",
                        },
                    ),
                    dbc.Col(
                        days[4],
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px",
                        },
                    ),
                    dbc.Col(
                        "infos",
                        width=4,
                        style={
                            # "border": "1px solid #ddd",
                            "padding": "0px",
                        },
                    ),
                ],
                className="grid-container",
                style={"height": "30%", "margin": "0px"},
            ),
            dbc.Row(
                [
                    dbc.Col(
                        create_map(1),
                        width=4,
                        style={"margin": "0px", "padding": "0px"},
                    ),
                    dbc.Col(
                        create_map(1),
                        width=4,
                        style={"margin": "0px", "padding": "0px"},
                    ),
                    dbc.Col(
                        "",
                        width=4,
                        style={"margin": "0px", "padding": "0px"},
                    ),
                ],
                style={"height": "20%", "margin": "0px"},
            ),
        ],
        style={
            "height": "2100px",
            "width": "1485px",
            "margin": "10px",
            "padding": "0",
            "fontSize": 11,
            "border": "1px solid black",
        },
    )
    app.run(debug=True)


def print_to_yaml(file, week, output_file):
    """Prints the solution to a yaml file"""
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

        total_distance = sol["total_distance"]
        fuel_consumption = sol["fuel_consumption"]
        tours_strkey = sol["tours"]
        tours = {}

        for key, tour in tours_strkey.items():
            d, v = map(int, key.split(", "))
            tours[d, v] = tour

        jours_map = {0: "Lundi", 1: "Mardi", 2: "Mercredi", 3: "Jeudi", 4: "Vendredi"}
        names_short = [
            "Toulouse/Seminaire",
            "Auterive",
            "Bagneres de Luchon",
            "Bessieres",
            "Blagnac",
            "Carbonne",
            "Cazeres",
            "Cugnaux",
            "Escalquens",
            "Fenouillet",
            "Fonsorbes",
            "Fronton",
            "L Isle en Dodon",
            "Leguevin",
            "Levignac",
            "Montrejeau",
            "Muret",
            "Pibrac",
            "Plaisance du Touch",
            "Portet sur Garonne",
            "Revel",
            "Rieumes",
            "Saint-Gaudens",
            "Salies du Salat/Mane",
            "Toulouse/Casselardit",
            "Toulouse/Malepere",
            "Toulouse/Negogousses",
            "Tournefeuille",
            "Villef. de Lauragais",
            "Auchan Gramont",
            "Leclerc St Orens",
            "Super U Flourens",
            "Leclerc Blagnac",
            "Leclerc Rouffiac",
            "Carrefour Centrale",
        ]

        vehicles_used = {}
        days = {}
        for d in range(pb.n_days):
            vehicles_used[d] = []
            if not any((d, v) in tours for v in range(pb.m)):
                continue

            days[jours_map[d]] = {}
            index_vehicule = 1
            for v in range(pb.m):
                if not (d, v) in tours:
                    continue

                vehicles_used[d].append(vehicles.index[v])
                tour = tours[d, v]

                tour_list = []
                for place in tour:
                    product_types = ""
                    palettes = ""
                    quantities = ""
                    if place["type"] == "livraison":
                        for i, product_type in enumerate("AFS"):
                            if place["delivery"][i] > 0:
                                quantities += (
                                    f" {str(place['delivery'][i])+product_type:>5}"
                                )
                                product_types += product_type
                            else:
                                quantities += "     "

                        pals = place["palettes"]
                        palettes = f"  {str(pals[0])+'A':>5}" if pals[0] else "      "
                        pals[1] = pals[1] / 2 if pals[1] % 2 != 0 else pals[1] // 2
                        palettes += f"{str(pals[1])+'F':>6}" if pals[1] else "     "
                        pals[2] = pals[2] / 2 if pals[2] % 2 != 0 else pals[2] // 2
                        palettes += f"{str(pals[2])+'S':>6}" if pals[2] else "     "
                        palettes += (
                            f"  +{place['norvegiennes']} norvegienne{'s' if place['norvegiennes'] > 1 else ''}"
                            if place["norvegiennes"]
                            else ""
                        )
                        tour_list.append(
                            {
                                f"Livraison - {names_short[place['index']]} {product_types}": {
                                    "Quantites": quantities,
                                    "Palettes": palettes,
                                }
                            }
                        )
                    else:
                        tour_list.append(
                            f"Ramasse - {names_short[place['index']]}",
                        )

                days[jours_map[d]][
                    f"{vehicles.index[v]} {vehicles_used[d].count(vehicles.index[v])}"
                ] = tour_list

        vehicles_used = {
            v: max(vehicles_used[d].count(v) for d in range(pb.n_days))
            for v in vehicles.index
        }

        output = {
            f"Tournees Semaine {week}": {
                "Consommation totale": f"{fuel_consumption:.2f}L",
                "Distance totale": f"{round(total_distance/1000):d}km",
                "Vehicules": " - ".join(f"{v} {k}" for k, v in vehicles_used.items()),
                "Tournees": days,
            }
        }

        with open(output_file, "w") as yaml_data:
            yaml.dump(output, yaml_data, default_flow_style=False)
            yaml.dump(output, sys.stdout, default_flow_style=False)


def print_to_txt(file, week, output_file):
    """Generates a pretty printed txt version of the solution and prints it to a specified file"""
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

        total_distance = sol["total_distance"]
        fuel_consumption = sol["fuel_consumption"]
        tours_strkey = sol["tours"]
        tours = {}

        for key, tour in tours_strkey.items():
            d, v = map(int, key.split(", "))
            tours[d, v] = tour

        jours_map = {0: "Lundi", 1: "Mardi", 2: "Mercredi", 3: "Jeudi", 4: "Vendredi"}

        output = ""
        output += f"- - - - - - TOURNEES SEMAINE {week} - - - - - -\n"
        vehicles_used = {}
        for d in range(pb.n_days):
            vehicles_used[d] = []
            if not any((d, v) in tours for v in range(pb.m)):
                continue

            output += f"\n- - - - {jours_map[d].upper()} - - - -\n\n"
            for v in range(pb.m):
                if not (d, v) in tours:
                    continue

                output += f"\tVéhicule {v} ({vehicles.index[v]}) \n"
                vehicles_used[d].append(vehicles.index[v])
                tour = tours[d, v]

                for place in tour:
                    product_types = ""
                    palettes = ""
                    if place["type"] == "livraison":
                        for i, product_type in enumerate("AFS"):
                            if place["delivery"][i] > 0:
                                product_types += (
                                    f"{str(place['delivery'][i])+product_type:5}"
                                )
                            else:
                                product_types += f"{'':5}"

                        pals = place["palettes"]
                        palettes = f"{pals[0]}PA" if pals[0] else "   "
                        palettes += f" {pals[1]/2:.1f}PF" if pals[1] else "      "
                        palettes += f" {pals[2]/2:.1f}PS" if pals[2] else "      "
                        palettes += (
                            f" +{place['norvegiennes']} norvégienne{'s' if place['norvegiennes'] > 1 else ''}"
                            if place["norvegiennes"]
                            else ""
                        )
                    else:
                        product_types = "Ramasse"

                    output += f"\t\t{place['name']:40}\t{product_types}\t{palettes}\n"

        output += f"\nConsommation totale : {fuel_consumption:.2f}L"
        output += f"\nDistance totale : {round(total_distance/1000):d}km"

        vehicles_used = {
            v: max(vehicles_used[d].count(v) for d in range(pb.n_days))
            for v in vehicles.index
        }
        output += f"\nVéhicules utilisés : {' - '.join(f'{v} {k}' for k, v in vehicles_used.items())}"

        with open(output_file, "w") as txt_file:
            txt_file.write(output)
        print(output)


def plot_solution(file, week, output_file, specific_day=None):
    """Makes a folium plot of the solution and saves it to a specified file"""
    with open(file) as f:
        sol = json.load(f)

        centres = pd.read_excel("data/centres.xlsx")
        vehicles = pd.read_excel("data/vehicules.xlsx", index_col=0)
        pdr = pd.read_excel("data/points_de_ramasse.xlsx")

        pb = read_problem(
            "data/centres.xlsx",
            "data/points_de_ramasse.xlsx",
            "data/vehicules.xlsx",
            "data/euclidean_matrix.xlsx",
            week,
        )

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

        obj = sol["total_distance"]
        tours_strkey = sol["tours"]
        tours = {}

        for key, tour in tours_strkey.items():
            d, v = map(int, key.split(", "))
            tours[d, v] = tour

        m = folium.Map(
            tiles="CartoDB Positron",
            location=(lats[0], longs[0]),
            max_bounds=True,
            min_lat=min(lats) - 0.5,
            max_lat=max(lats) + 0.5,
            min_lon=min(longs) - 1,
            max_lon=max(longs) + 1,
        )

        week = int(week)

        # Add the centres and points de ramasse as markers
        for i in range(len(centres)):
            folium.CircleMarker(
                location=(lats[i], longs[i]),
                tooltip=centres["Nom"][i],
                color="green" if i == 0 else "red",
            ).add_to(m)
        for i in range(len(pdr)):
            folium.CircleMarker(
                location=(lats[pb.n + i], longs[pb.n + i]),
                tooltip=pdr["Nom"][i],
                color="blue",
            ).add_to(m)

        vehicles_used = {}
        days_to_plot = [specific_day] if specific_day is not None else range(pb.n_days)
        for d in days_to_plot:
            vehicles_used[d] = []
            if not any((d, v) in tours for v in range(pb.m)):
                continue

            group = folium.FeatureGroup(jours_map[d]).add_to(m)

            for v in range(pb.m):
                if not (d, v) in tours:
                    continue

                vehicles_used[d].append(vehicles.index[v])
                tour = tours[d, v]

                tour_coords = []
                for place in tour:
                    index, name = place["index"], place["name"]
                    tour_coords.append((lats[index], longs[index]))
                tour_coords = [coords_depot] + tour_coords + [coords_depot]

                folium.PolyLine(
                    tour_coords,
                    tooltip=f"Véhicule {v} ({vehicles.index[v]})",
                    color=colors[v],
                ).add_to(group)

        folium.LayerControl().add_to(m)

        m.save(output_file)
        print("Map saved to ", output_file)

        return m


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="json solution file")
    parser.add_argument("week", type=int, help="week number (1 or 2)")
    parser.add_argument("outfile", type=str, help="desired output file path")

    args = parser.parse_args()

    file_name = args.outfile.split(".")[0]

    plot_solution(args.infile, args.week, file_name + ".html")
    print_to_txt(args.infile, args.week, file_name + ".txt")
    print_to_yaml(args.infile, args.week, file_name + ".yaml")
    # make_dashboard("solutions/week_1.yaml", 1, "dashboard.html")
