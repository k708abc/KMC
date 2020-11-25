#! /usr/bin/env python3
"""初期値も予め入力したい方式なら、クラス作っちゃった方が楽かも"""

from typing import OrderedDict


class Params:
    """Class for input parameters

    Attributes
    -------------

    """

    kb_eV = 8.617e-5

    def __init__(self) -> None:
        self.n_cell_init = 5
        self.z_unit_init = 5
        self.temperature = 550
        self.dep_rate = 0.4
        self.dep_time = 5
        self.post_anneal = 0
        self.prefactor = "1E+13"
        self.binding_energies: OrderedDict[str, float] = OrderedDict()
        self.binding_energies["AgSi"] = -1.5
        self.binding_energies["Si12"] = -1.5
        self.binding_energies["Si23"] = -1.5
        self.binding_energies["Si34"] = -1.5
        self.binding_energies["Si45"] = -1.5
        self.binding_energies["Si56"] = -1.5
        self.binding_energies["Si_intra"] = -1.5
        self.binding_energies["Si_inter"] = -1.5
        self.binding_energies["Agtop"] = -1.5
        self.transformation = -0.3
        self.record_name = "KMC_Si_rec"
        self.img_per = 10
        self.comments = "No comments"
        self.intra_distance = 0.204
        self.inter_distance = 0.612

    @property
    def temperature_eV(self) -> float:
        return self.temperature * 8.617e-5

    @property
    def atoms_in_BL(self) -> int:
        return 2 * (self.n_cell_init) ** 2

    @property
    def dep_rate_atoms_persec(self) -> float:
        return self.dep_rate / 60 * self.atoms_in_BL

    @property
    def total_time(self) -> float:
        return (self.dep_time + self.post_anneal) * 60

    @property
    def interval(self) -> float:
        return self.total_time * self.img_per / 100
