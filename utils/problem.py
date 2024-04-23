from dataclasses import dataclass
import json
import pandas as pd
from math import ceil
from utils.prepare_demands import prepare_demands
from utils.datatypes import (
    DeliveryWeek,
    ProductType,
    StopType,
    Stop,
    Vehicle,
    Centre,
    PDR,
    Demand,
    Params,
)


@dataclass
class Problem:
    distance_matrix: pd.DataFrame
    duration_matrix: pd.DataFrame
    no_traffic_duration_matrix: pd.DataFrame
    n_centres: int
    n_pdr: int
    m: int
    n_days: int
    vehicles: list[Vehicle]
    centres: list[Centre]
    pdrs: list[PDR]
    demands: dict[int, list[Demand]]
    livraisons_de_ramasses: set[tuple[int, int, int, int]]
    params: Params

    def __eq__(self, other):
        if not isinstance(other, Problem):
            return False
        return (
            self.distance_matrix.equals(other.distance_matrix)
            and self.duration_matrix.equals(other.duration_matrix)
            and self.no_traffic_duration_matrix.equals(other.no_traffic_duration_matrix)
            and self.n_centres == other.n_centres
            and self.n_pdr == other.n_pdr
            and self.m == other.m
            and self.n_days == other.n_days
            and self.vehicles == other.vehicles
            and self.centres == other.centres
            and self.pdrs == other.pdrs
            and self.demands == other.demands
            and self.livraisons_de_ramasses == other.livraisons_de_ramasses
            and self.params == other.params
        )

    def to_json(self, path):
        """Writes the problem to a json file"""
        dict = {
            "distance_matrix": self.distance_matrix.values.tolist(),
            "duration_matrix": self.duration_matrix.values.tolist(),
            "no_traffic_duration_matrix": self.no_traffic_duration_matrix.values.tolist(),
            "n_centres": self.n_centres,
            "n_pdr": self.n_pdr,
            "m": self.m,
            "n_days": self.n_days,
            "vehicles": [v.to_dict() for v in self.vehicles],
            "centres": [c.to_dict() for c in self.centres],
            "pdrs": [p.to_dict() for p in self.pdrs],
            "demands": {
                c: [d.to_dict() for d in demands] for c, demands in self.demands.items()
            },
            "livraisons_de_ramasses": list(self.livraisons_de_ramasses),
            "params": self.params.to_dict(),
        }
        json.dump(dict, open(path, "w"))

    @classmethod
    def from_json(cls, path):
        dict = json.load(open(path))
        return cls(
            pd.DataFrame(dict["distance_matrix"]),
            pd.DataFrame(dict["duration_matrix"]),
            pd.DataFrame(dict["no_traffic_duration_matrix"]),
            dict["n_centres"],
            dict["n_pdr"],
            dict["m"],
            dict["n_days"],
            [Vehicle.from_dict(v) for v in dict["vehicles"]],
            [Centre.from_dict(c) for c in dict["centres"]],
            [PDR.from_dict(p) for p in dict["pdrs"]],
            {
                int(c): [Demand.from_dict(d) for d in demands]
                for c, demands in dict["demands"].items()
            },
            set([(d, v, trip, c) for d, v, trip, c in dict["livraisons_de_ramasses"]]),
            Params.from_dict(dict["params"]),
        )

    @classmethod
    def from_files(
        cls,
        demands_file,
        allowed_days_file,
        week_assignments_file,
        points_de_ramasse_file,
        vehicles_file,
        distance_matrix_file,
        duration_matrix_file,
        no_traffic_duration_matrix_file,
        params_file,
    ):
        """Reads the problem from the usual input files and returns a Problem object"""

        demands_df = pd.read_excel(demands_file, index_col=0)
        allowed_days_df = pd.read_excel(allowed_days_file)
        week_assignments_df = pd.read_excel(week_assignments_file)
        points_de_ramasse_df = pd.read_excel(points_de_ramasse_file, index_col=0)
        vehicles_df = pd.read_excel(vehicles_file, index_col=0)
        distance_matrix_df = pd.read_excel(distance_matrix_file, index_col=0)
        duration_matrix_df = pd.read_excel(duration_matrix_file, index_col=0)
        no_traffic_duration_matrix_df = pd.read_excel(
            no_traffic_duration_matrix_file, index_col=0
        )

        n_centres = len(demands_df.index)
        n_pdr = len(points_de_ramasse_df.index)
        n_days = 5  # Single week scheduling
        m = len(vehicles_df.index)

        params_json = json.load(open(params_file, "r"))
        params = Params.from_dict(params_json)

        ## Process the vehicles info

        capacities = vehicles_df["Capacité (kg)"].astype(int).to_list()
        sizes = vehicles_df["Taille(Palettes)"].astype(int).to_list()
        frais = vehicles_df["Réfrigérant"].to_list()
        consumptions = vehicles_df["Consommation (L/100km)"].astype(int).to_list()
        frais_entretien = vehicles_df["Frais Entretien (€/km)"].astype(float).to_list()
        frais_ct = (
            vehicles_df["Frais Contrôle Technique (€/an)"].astype(float).to_list()
        )
        frais_assurance = vehicles_df["Frais Assurance (€/an)"].astype(float).to_list()
        costs_per_km = [
            params_json["fuel_cost"] * consumptions[v] / 100 + frais_entretien[v]
            for v in range(m)
        ]
        fixed_costs = [frais_ct[v] / 52 + frais_assurance[v] / 52 for v in range(m)]

        vehicles = [
            Vehicle(
                i,
                vehicles_df.index[i],
                bool(params_json["vehicle_allowed"][i]),
                capacities[i],
                consumptions[i],
                sizes[i],
                can_carry={ProductType.A, ProductType.F, ProductType.S},
                allows_isotherm_cover=frais[i] == "Oui",
                cost_per_km=costs_per_km[i],
                fixed_cost=fixed_costs[i],
            )
            for i in range(m)
        ]
        for v in range(m):
            if not frais[v] == "Oui":
                vehicles[v].can_carry.remove(ProductType.F)
        vehicles[0].can_carry.remove(
            ProductType.S
        )  # The PL (index 0) can't carry S products

        ## Process the centres info

        jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4}

        jours_de_livraison = allowed_days_df["AllowedDays"].tolist()
        jours_de_livraison = [[]] + [
            j.replace(" ", "").split(",") for j in jours_de_livraison[1:]
        ]
        j_de_livraison_possibles = [
            set(map(lambda x: jours_map[x], j)) for j in jours_de_livraison
        ]

        demands = prepare_demands(demands_file, params)

        centres = [
            Centre(
                i,
                week_assignments_df["Centre"][i],
                j_de_livraison_possibles[i],
                DeliveryWeek(week_assignments_df["Week"][i]),
            )
            for i in range(n_centres)
        ]

        weights = (
            points_de_ramasse_df["Poids par ramasse(kg)"].fillna(0).astype(int).tolist()
        )
        for p in range(n_pdr):
            weights[p] = ceil(weights[p])

        # Mandatory pickup days
        jours_de_ramasse = points_de_ramasse_df["Jours de Ramasse"].tolist()
        j_de_ramasse = []
        for pdr in range(n_pdr):
            jours = jours_de_ramasse[pdr].split(", ")
            for i, j in enumerate(jours):
                jours[i] = jours_map[jours[i]]
            j_de_ramasse.append(jours)

        product_types = points_de_ramasse_df["Type de produits"].to_list()
        pdrs = [
            PDR(
                i,
                points_de_ramasse_df.index[i],
                j_de_ramasse[i],
                weights[i],
                ProductType[product_types[i]],
            )
            for i in range(n_pdr)
        ]

        index_arenes = 25
        index_malepere = 26
        livraisons_de_ramasses = {
            (0, 2, 0, index_malepere),
            (1, 2, 0, index_malepere),
            (2, 2, 0, index_malepere),
            (4, 2, 0, index_malepere),
            (4, 3, 1, index_arenes),
        }

        problem = cls(
            distance_matrix_df,
            duration_matrix_df,
            no_traffic_duration_matrix_df,
            n_centres,
            n_pdr,
            m,
            n_days,
            vehicles,
            centres,
            pdrs,
            demands,
            livraisons_de_ramasses,
            params,
        )

        return problem

    @property
    def allowed_vehicles(self):
        return [v for v in range(self.m) if self.vehicles[v].allowed]

    @property
    def n_nodes(self):
        return self.n_centres + self.n_pdr

    def delivered_centres(self, week):
        return [
            c
            for c in range(1, self.n_centres)
            if self.centres[c].delivery_week in [DeliveryWeek.ANY, week]
        ]

    @property
    def pdr_nodes(self):
        return range(self.n_centres, self.n_centres + self.n_pdr)

    def used_nodes(self, week):
        return [0] + self.delivered_centres(week) + list(self.pdr_nodes)


@dataclass
class Solution:
    week: DeliveryWeek
    total_costs: float
    fixed_costs: float
    variable_costs: float
    total_distance: int
    vehicles_used: list[bool]
    tours: dict[tuple[int, int], list[Stop]]
    tour_durations: dict[tuple[int, int], int]

    def to_json(self, outfile):
        sol = {
            "week": self.week.value,
            "total_costs": self.total_costs,
            "fixed_costs": self.fixed_costs,
            "variable_costs": self.variable_costs,
            "total_distance": self.total_distance,
            "vehicles_used": self.vehicles_used,
        }
        sol["tours"] = {
            str(d) + ", " + str(v): [stop.to_dict() for stop in tour]
            for (d, v), tour in self.tours.items()
        }
        sol["tour_durations"] = {
            str(d) + ", " + str(v): int(round(td))
            for (d, v), td in self.tour_durations.items()
        }

        json.dump(sol, open(outfile, "w"))
        print(f"Wrote solution to {outfile}")

    @classmethod
    def from_json(cls, file_path):
        data = json.load(open(file_path))

        def key_to_ints(key_str: str):
            return (
                int(key_str.split(", ")[0]),
                int(key_str.split(", ")[1]),
            )

        tour_durations = {key_to_ints(k): v for k, v in data["tour_durations"].items()}

        tours = {
            key_to_ints(k): [Stop.from_dict(stop) for stop in tour]
            for k, tour in data["tours"].items()
        }

        sol = cls(
            DeliveryWeek(data["week"]),
            data["total_costs"],
            data["fixed_costs"],
            data["variable_costs"],
            data["total_distance"],
            data["vehicles_used"],
            tours,
            tour_durations,
        )
        return sol

    def adjust_durations(self, pb: Problem):
        def duration(a, b, tour_duration):
            if tour_duration <= 90 * 60:
                arc_duration = pb.duration_matrix.iloc[a, b]
            else:
                arc_duration = pb.no_traffic_duration_matrix.iloc[a, b]
            return arc_duration

        for (d, v), tour in self.tours.items():
            tour_duration = 0
            trip = 0
            a = 0
            for stop in tour:
                if stop.index == 0:
                    trip += 1
                    # Trip back to depot
                    arc_duration = duration(a, 0, tour_duration)
                    tour_duration += arc_duration  # type: ignore
                    a = 0
                    trip = 0
                    continue

                if trip > 0 and a == 0:
                    tour_duration += pb.params.wait_between_trips * 60

                b = stop.index
                arc_duration = duration(a, b, tour_duration)
                tour_duration += arc_duration  # type: ignore
                if stop.type == StopType.Livraison:
                    tour_duration += pb.params.wait_at_centres * 60
                elif stop.type == StopType.Ramasse:
                    tour_duration += pb.params.wait_at_pdrs * 60
                a = b

            # Trip back to depot
            arc_duration = duration(a, 0, tour_duration)

            tour_duration += arc_duration  # type: ignore

            self.tour_durations[d, v] = round(tour_duration)
