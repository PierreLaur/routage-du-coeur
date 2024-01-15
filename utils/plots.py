import pandas as pd
import plotly.graph_objects as go
import sys
import json
from problem import Problem, read_problem
import folium
import yaml
import argparse


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

        obj = sol["total_distance"]
        tours_strkey = sol["tours"]
        tours = {}

        for key, tour in tours_strkey.items():
            d, v = map(int, key.split(", "))
            tours[d, v] = tour

        jours_map = {0: "Lundi", 1: "Mardi", 2: "Mercredi", 3: "Jeudi", 4: "Vendredi"}

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
                                f"Livraison - {place['name']} {product_types}": {
                                    "Quantites": quantities,
                                    "Palettes": palettes,
                                }
                            }
                        )
                    else:
                        tour_list.append(
                            f"Ramasse - {place['name']}",
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
                "Distance totale": f"{round(obj/1000):d}km",
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

        obj = sol["total_distance"]
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
                        palettes += f" {pals[1]}/2PF" if pals[1] else "      "
                        palettes += f" {pals[2]}/2PS" if pals[2] else "      "
                        palettes += (
                            f" +{place['norvegiennes']} norvégienne{'s' if place['norvegiennes'] > 1 else ''}"
                            if place["norvegiennes"]
                            else ""
                        )
                    else:
                        product_types = "Ramasse"

                    output += f"\t\t{place['name']:40}\t{product_types}\t{palettes}\n"

        output += f"\nDistance totale : {round(obj/1000):d}km"

        vehicles_used = {
            v: max(vehicles_used[d].count(v) for d in range(pb.n_days))
            for v in vehicles.index
        }
        output += f"\nVéhicules utilisés : {' - '.join(f'{v} {k}' for k, v in vehicles_used.items())}"

        with open(output_file, "w") as txt_file:
            txt_file.write(output)
        print(output)


def plot_solution(file, week, output_file):
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
            "beige",
            "beige",
            "beige",
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
            folium.Marker(
                location=(lats[i], longs[i]),
                tooltip=centres["Nom"][i],
                icon=folium.Icon(color="green" if i == 0 else "red"),
            ).add_to(m)
        for i in range(len(pdr)):
            folium.Marker(
                location=(lats[pb.n + i], longs[pb.n + i]),
                tooltip=pdr["Nom"][i],
                icon=folium.Icon(color="blue"),
            ).add_to(m)

        vehicles_used = {}
        for d in range(pb.n_days):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="json solution file")
    parser.add_argument("week", type=int, help="week number (1 or 2)")
    parser.add_argument("outfile", type=str, help="desired output file path")

    args = parser.parse_args()

    extension = args.outfile.split(".")[-1]

    if extension == "html":
        plot_solution(args.infile, args.week, args.outfile)
    elif extension == "txt":
        print_to_txt(args.infile, args.week, args.outfile)
    elif extension == "yaml":
        print_to_yaml(args.infile, args.week, args.outfile)
    else:
        print("File extension must be html, txt or yaml")
