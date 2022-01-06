#! /usr/bin/env python3

import decimal
from typing import Union
import yaml
from pathlib import Path


class Params:
    """Class for input parameters

    Attributes
    -------------


    """

    kb_eV: float = 8.617e-5

    def __init__(self, filename: Union[str, Path] = "") -> None:
        if filename:
            with open(filename) as f:
                input_yaml = yaml.safe_load(f)
            for k, v in input_yaml.items():
                setattr(self, k, v)

    @property
    def temperature_eV(self) -> float:
        return self.temperature * 8.617e-5

    @property
    def atoms_in_BL(self) -> int:
        return int(2 * (self.cell_size_xy) ** 2)

    @property
    def dep_rate_atoms_persec(self) -> decimal.Decimal:
        return decimal.Decimal(self.dep_rate / 60 * self.atoms_in_BL)

    @property
    def total_time(self) -> decimal.Decimal:
        return decimal.Decimal(self.dep_time) * 60

    @property
    def interval(self) -> decimal.Decimal:
        return self.total_time * decimal.Decimal(self.img_per) / 100

    @property
    def total_atoms(self) -> int:
        return int(round(self.dep_rate_atoms_persec * self.total_time))

    @property
    def rec_num_atom_interval(self) -> int:
        return int(round(self.total_atoms * self.img_per / 100))
