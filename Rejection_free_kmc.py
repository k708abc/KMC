from typing import List, Dict, Tuple
from kmc_functions import common_functions
import copy
import math
from cal_rates import rate
import time
from lattice_form import lattice_form
from InputParameter import Params
from deposition import deposit_an_atom
from choose_site import choose_atom
from judgement import judge_null
from event_collection import site_events
from normarize_list import normarize_rate
from weighted_choice import choice
from recording import record_data
from atoms_recalculate import recalculate
from rejection_free_choose import rejection_free_choise


class rejetion_free(common_functions):
    def __init__(self) -> None:
        None

    def start(self):
        self.init_value = Params()
        self.start_setting()
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        self.atom_exist: List[Tuple[int, int, int]] = []
        (
            self.lattice,
            self.bonds,
            self.atom_set,
            self.event,
            self.event_time,
            self.event_time_tot,
            self.event_state,
        ) = lattice_form(self.init_value)
        # 最初の二原子を配置
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.put_first)):
                self.rejection_free_deposition()
        #
        self.prev_eve = "dep"