from dataclasses import dataclass, asdict
import pandas as pd
from math import ceil
import json
from enum import Enum


class ProductType(Enum):
    A = 0
    F = 1
    S = 2


class DeliveryWeek(Enum):
    ANY = 0
    ODD = 1
    EVEN = 2


@dataclass
class Centre:
    index: int
    name: str
    allowed_days: set[int]
    delivery_week: DeliveryWeek
    demands: list[dict[ProductType, int]]

    def to_dict(self):
        d = asdict(self)
        d["allowed_days"] = list(self.allowed_days)
        d["delivery_week"] = self.delivery_week.name
        d["demands"] = [
            {t.name: int(demand) for t, demand in scenario.items()}
            for scenario in self.demands
        ]
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["index"],
            d["name"],
            set(d["allowed_days"]),
            DeliveryWeek[d["delivery_week"]],
            [
                {ProductType[k]: v for k, v in scenario.items()}
                for scenario in d["demands"]
            ],
        )

    def __str__(self) -> str:
        return f"Centre {self.index} {self.name} [allowed_days={self.allowed_days} delivery_week={self.delivery_week.name} demands={[[str(demand)+t.name for t, demand in scenario.items()] for scenario in self.demands]}"

    def __repr__(self) -> str:
        return f"Centre {self.index} {self.name} [allowed_days={self.allowed_days} delivery_week={self.delivery_week.name} demands={[[str(demand)+t.name for t, demand in scenario.items()] for scenario in self.demands]}"


@dataclass
class PDR:
    index: int
    name: str
    required_days: set[int]
    weight: int
    product_type: ProductType

    def to_dict(self):
        d = asdict(self)
        d["required_days"] = list(self.required_days)
        d["product_type"] = self.product_type.name
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["index"],
            d["name"],
            set(d["required_days"]),
            d["weight"],
            ProductType[d["product_type"]],
        )


@dataclass
class Vehicle:
    index: int
    name: str
    allowed: bool
    capacity: int
    consumption: int
    size: int
    can_carry: set[ProductType]
    allows_isotherm_cover: bool
    cost_per_km: float
    fixed_cost: float

    def to_dict(self):
        d = asdict(self)
        d["can_carry"] = [p.name for p in self.can_carry]
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["index"],
            d["name"],
            d["allowed"],
            d["capacity"],
            d["consumption"],
            d["size"],
            set([ProductType[p] for p in d["can_carry"]]),
            d["allows_isotherm_cover"],
            d["cost_per_km"],
            d["fixed_cost"],
        )


@dataclass
class Params:
    week: DeliveryWeek
    max_palette_capacity: int
    demi_palette_capacity: int
    n_norvegiennes: int
    norvegienne_capacity: int
    max_stops: int
    max_trips: int
    max_tour_duration: int
    max_tour_duration_with_pickup: int
    wait_at_centres: int
    wait_at_pdrs: int
    fuel_cost: float

    def to_dict(self):
        d = asdict(self)
        d["week"] = self.week.name
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            DeliveryWeek[d["week"]],
            d["max_palette_capacity"],
            d["demi_palette_capacity"],
            d["n_norvegiennes"],
            d["norvegienne_capacity"],
            d["max_stops"],
            d["max_trips"],
            d["max_tour_duration"],
            d["max_tour_duration_with_pickup"],
            d["wait_at_centres"],
            d["wait_at_pdrs"],
            d["fuel_cost"],
        )


class StopType(Enum):
    Livraison = 0
    Ramasse = 1


@dataclass
class Stop:
    index: int
    name: str
    type: StopType
    delivery: tuple[int, int, int] = (0, 0, 0)
    palettes: tuple[int, float, float] = (0, 0, 0)
    norvegiennes: int = 0

    def to_dict(self):
        d = {
            "index": self.index,
            "name": self.name,
            "type": self.type.name,
        }
        if self.type == StopType.Livraison:
            d["delivery"] = self.delivery
            d["palettes"] = self.palettes
            d["norvegiennes"] = self.norvegiennes
        return d

    @classmethod
    def from_dict(cls, dict):
        stoptype = StopType[dict["type"]]
        stop = cls(
            dict["index"],
            dict["name"],
            stoptype,
        )
        if stoptype == StopType.Livraison:
            stop.delivery = dict["delivery"]
            stop.palettes = dict["palettes"]
            stop.norvegiennes = dict["norvegiennes"]
        return stop


@dataclass
class Solution:
    week: int
    total_costs: float
    fixed_costs: float
    variable_costs: float
    total_distance: int
    vehicles_used: list[bool]
    tours: dict[tuple[int, int], list[Stop]]
    tour_durations: dict[tuple[int, int], int]
    tour_durations_adjusted: dict[tuple[int, int], str] | None = None

    def write_as_json(self, outfile):
        sol = {
            "week": self.week,
            "total_costs": round(self.total_costs, 5),
            "fixed_costs": round(self.fixed_costs, 5),
            "variable_costs": round(self.variable_costs, 5),
            "total_distance": self.total_distance,
            "vehicles_used": self.vehicles_used,
        }
        sol["tours"] = {
            str(d) + ", " + str(v): [stop.to_dict() for stop in tour]
            for (d, v), tour in self.tours.items()
        }
        sol["tour_durations"] = {
            str(d) + ", " + str(v): int(td)
            for (d, v), td in self.tour_durations.items()
        }
        if self.tour_durations_adjusted is not None:
            sol["tour_durations_adjusted"] = {
                str(d) + ", " + str(v): td
                for (d, v), td in self.tour_durations_adjusted.items()
            }

        json.dump(sol, open(outfile, "w"))
        print(f"Wrote solution to {outfile}")

    @classmethod
    def read_from_json(cls, file_path):
        data = json.load(open(file_path))

        def key_to_ints(key_str):
            return (
                int(key_str.split(", ")[0]),
                int(key_str.split(", ")[1]),
            )

        tour_durations = {key_to_ints(k): v for k, v in data["tour_durations"].items()}
        if "tour_durations_adjusted" in data:
            tour_durations_adjusted = {
                key_to_ints(k): v for k, v in data["tour_durations_adjusted"].items()
            }
        else:
            tour_durations_adjusted = None

        tours = {
            key_to_ints(k): [Stop.from_dict(stop) for stop in tour]
            for k, tour in data["tours"].items()
        }

        sol = cls(
            data["week"],
            data["total_costs"],
            data["fixed_costs"],
            data["variable_costs"],
            data["total_distance"],
            data["vehicles_used"],
            tours,
            tour_durations,
            tour_durations_adjusted,
        )
        return sol

    def adjust_durations(
        self,
        traffic_matrix: pd.DataFrame,
        no_traffic_matrix: pd.DataFrame,
        params: Params,
    ):
        if self.tour_durations_adjusted is not None:
            return

        self.tour_durations_adjusted = {}
        for (d, v), tour in self.tours.items():
            tour_indices = [stop.index for stop in tour]

            estimated_duration = 0
            for a, b in zip([0] + tour_indices, tour_indices + [0]):
                if estimated_duration <= 90:
                    estimated_duration += traffic_matrix.iloc[a, b] / 60  # type: ignore
                else:
                    estimated_duration += no_traffic_matrix.iloc[a, b] / 60  # type: ignore

            estimated_duration += sum(
                params.wait_at_centres * (stop.type == StopType.Livraison)
                + params.wait_at_pdrs * (stop.type == StopType.Ramasse)
                for stop in tour
            )

            self.tour_durations_adjusted[d, v] = (
                f"{estimated_duration//60:.0f}h{estimated_duration%60:.0f}min"
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
    livraisons_de_ramasses: set[tuple[int, int, int, int]]
    params: Params

    def write_as_json(self, path):
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
            set([(d, v, trip, c) for d, v, trip, c in dict["livraisons_de_ramasses"]]),
            Params.from_dict(dict["params"]),
        )

    @classmethod
    def from_files(
        cls,
        demands_files,
        allowed_days_file,
        week_assignments_file,
        points_de_ramasse_file,
        vehicles_file,
        distance_matrix_file,
        duration_matrix_file,
        no_traffic_duration_matrix_file,
        params_file,
        week,
        custom_assignment=None,
    ):
        """Reads the problem from the usual input files and returns a Problem object"""

        demands_dfs = [
            pd.read_excel(demands_file, index_col=0) for demands_file in demands_files
        ]
        allowed_days_df = pd.read_excel(allowed_days_file)
        week_assignments_df = pd.read_excel(week_assignments_file)
        points_de_ramasse_df = pd.read_excel(points_de_ramasse_file, index_col=0)
        vehicles_df = pd.read_excel(vehicles_file, index_col=0)
        distance_matrix_df = pd.read_excel(distance_matrix_file, index_col=0)
        duration_matrix_df = pd.read_excel(duration_matrix_file, index_col=0)
        no_traffic_duration_matrix_df = pd.read_excel(
            no_traffic_duration_matrix_file, index_col=0
        )
        params_json = json.load(open(params_file, "r"))

        n_centres = len(demands_dfs[0].index)
        n_pdr = len(points_de_ramasse_df.index)
        n_days = 5  # Single week scheduling
        m = len(vehicles_df.index)

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

        demands = {}
        for i, demands_df in enumerate(demands_dfs):
            demands_df *= params_json["robustness_factor"]
            demands_df = demands_df.astype(int)
            for c in range(n_centres):
                demands[i, c] = {
                    ProductType.A: demands_df.iloc[c, 0],
                    ProductType.F: demands_df.iloc[c, 1],
                    ProductType.S: demands_df.iloc[c, 2],
                }

        centres = [
            Centre(
                i,
                week_assignments_df["Centre"][i],
                j_de_livraison_possibles[i],
                DeliveryWeek(week_assignments_df["Week"][i]),
                [demands[scenario, i] for scenario in range(len(demands_dfs))],
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

        index_malepere = 26
        livraisons_de_ramasses = {
            (0, 2, 0, index_malepere),
            (1, 2, 0, index_malepere),
            (2, 2, 0, index_malepere),
            (4, 2, 0, index_malepere),
        }

        params = Params(
            DeliveryWeek(week),
            params_json["max_palette_capacity"],
            params_json["demi_palette_capacity"],
            params_json["n_norvegiennes"],
            params_json["norvegienne_capacity"],
            params_json["max_stops"],
            params_json["max_trips"],
            params_json["max_tour_duration"],
            params_json["max_tour_duration_with_pickup"],
            params_json["wait_at_centres"],
            params_json["wait_at_pdrs"],
            params_json["fuel_cost"],
        )

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

    @property
    def delivered_centres(self):
        return [
            c
            for c in range(1, self.n_centres)
            if self.centres[c].delivery_week in [DeliveryWeek.ANY, self.params.week]
        ]

    @property
    def pdr_nodes(self):
        return range(self.n_centres, self.n_centres + self.n_pdr)

    @property
    def used_nodes(self):
        return [0] + self.delivered_centres + list(self.pdr_nodes)
