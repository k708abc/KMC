﻿#! /usr/bin/env python3

from collections import OrderedDict
import math
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
        """
        self.n_cell_init = 10
        self.z_unit_init = 5
        self.temperature = 550.0
        self.dep_rate = 0.1
        self.dep_time = 20.0
        self.prefactor = 1000000000000.0
        self.binding_energies: OrderedDict[str, float] = OrderedDict()
        self.binding_energies["Si_first"] = 0.2
        self.binding_energies["Si_second"] = 0.2
        self.binding_energies["Si_third"] = 0.2
        self.binding_energies["Si_else"] = 0.2
        self.diffusion_barriers: OrderedDict[str, float] = OrderedDict()
        self.diffusion_barriers["Si_first"] = 1.06
        self.diffusion_barriers["Si_second"] = 1.0
        self.diffusion_barriers["Si_third"] = 1.3
        self.diffusion_barriers["Si_else"] = 1.3
        #
        self.put_first = 10
        self.cut_number = 100000
        self.num_defect = 1
        self.record_name = "Test3"
        self.img_per = 5.0
        self.comments = "Test"
        self.intra_distance = 0.204
        self.inter_distance = 0.612
        self.keep_defect_check = False
        self.first_put_check = True
        self.cut_check = False
        self.limit_check = False
        self.limit_val = 1000
        self.subtract_check = True
        self.transformation = True
        self.trans_val = 5
        self.method = "Rejection_free"

        # for repetition
        self.repeat_combos = (
            "Bond E1: 2nd, E2: multi",
            "Diffuse E1: 2nd, E2: multi",
        )
        self.repeat_val = 1
        self.start_E1 = 1.15
        self.end_E1 = 1.25
        self.diff_E1 = 0.05
        self.start_E2 = 1.15
        self.end_E2 = 1.25
        self.diff_E2 = 0.05
        self.max_workers_val = 3
        """

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

    @property
    def cal_energy(self) -> float:
        return (
            -8.617e-5
            * self.temperature
            * math.log(float(self.limit_num) / float(self.prefactor))
        )