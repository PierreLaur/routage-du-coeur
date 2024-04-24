from dataclasses import asdict, dataclass
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
class Demand:
    weight: int
    palettes: float
    norvegiennes: int
    product_type: ProductType

    def to_dict(self):
        d = asdict(self)
        d["product_type"] = self.product_type.name
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["weight"],
            d["palettes"],
            d["norvegiennes"],
            ProductType[d["product_type"]],
        )


@dataclass
class Centre:
    index: int
    name: str
    allowed_days: set[int]
    delivery_week: DeliveryWeek

    def to_dict(self):
        d = asdict(self)
        d["allowed_days"] = list(self.allowed_days)
        d["delivery_week"] = self.delivery_week.name
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["index"],
            d["name"],
            set(d["allowed_days"]),
            DeliveryWeek[d["delivery_week"]],
        )

    def __str__(self) -> str:
        return f"Centre {self.index} {self.name} [allowed_days={self.allowed_days} delivery_week={self.delivery_week.name}"

    def __repr__(self) -> str:
        return f"Centre {self.index} {self.name} [allowed_days={self.allowed_days} delivery_week={self.delivery_week.name}"


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
    wait_between_trips: int
    fuel_cost: float

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(
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
            d["wait_between_trips"],
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
    palettes: tuple[float, float, float] = (0, 0, 0)
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
            stop.delivery = tuple(dict["delivery"])
            stop.palettes = tuple(dict["palettes"])
            stop.norvegiennes = dict["norvegiennes"]
        return stop
