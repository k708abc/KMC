from InputParameter import Params
import time
from typing import List, Dict, Tuple
from Modules.lattice_form import lattice_form
from Modules.deposition import deposit_an_atom
from Modules.atoms_recalculate import recalculate
from Modules.rejection_free_choose import rejection_free_choise
from Modules.event_collection import site_events
import math
from Modules.recording import record_data
import copy

##
from Modules.recording import rec_events_per_dep


class common_functions:
    def __init__(self) -> None:
        self.init_value = Params()

    def start_setting(self):
        self.start_time = time.time()
        self.prog_time = 0
        self.n_atoms = 0
        self.n_events = 0
        self.empty_firstBL = self.init_value.atoms_in_BL
        self.pos_rec: List[dict] = []
        self.time_rec: List[float] = []
        self.cov_rec: List[float] = []
        self.rec_num_atoms = 0
        self.n_events_perdep = 0
        self.n_events_rec: List[int] = []
        self.num_atoms_rec: List[int] = []

    def update_progress(self) -> None:
        self.n_events += 1
        self.n_events_perdep += 1

    def start_rejection_free(self):
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        self.atom_exist: List[Tuple[int, int, int]] = []
        (
            self.lattice,
            self.bonds,
            self.atom_set,
            self.event,
            self.event_time,
            self.event_time_tot,
            _,
        ) = lattice_form(self.init_value)
        self.energy_summarize()

    def energy_summarize(self):
        self.energy_bonding = [
            self.init_value.binding_energies["Si_1st"],
            self.init_value.binding_energies["Si_2nd"],
            self.init_value.binding_energies["Si_3rd"],
        ]
        #
        self.energy_bonding += [
            self.init_value.binding_energies["Si_else"]
            for i in range(self.init_value.z_unit_init * 3 - 2)
        ]
        #
        """
        self.energy_2D3D = [
            (
                self.init_value.binding_energies["Si_1st"]
                + self.init_value.binding_energies["Si_else"]
            )
            / 2,
            (
                self.init_value.binding_energies["Si_2nd"]
                + self.init_value.binding_energies["Si_else"]
            )
            / 2,
            (
                self.init_value.binding_energies["Si_3rd"]
                + self.init_value.binding_energies["Si_else"]
            )
            / 2,
        ]
        self.energy_2D3D += [
            self.init_value.binding_energies["Si_else"]
            for i in range(self.init_value.z_unit_init * 3 - 2)
        ]
        #
        self.energy_3D = [
            self.init_value.binding_energies["Si_else"]
            for i in range(self.init_value.z_unit_init * 3 + 1)
        ]
        """
        #
        self.energy_diffuse = [
            self.init_value.diffusion_barriers["Si_1st"],
            self.init_value.diffusion_barriers["Si_2nd"],
            self.init_value.diffusion_barriers["Si_3rd"],
            self.init_value.diffusion_barriers["Si_else"],
        ]
        #
        self.energy_diffuse += [
            self.init_value.diffusion_barriers["Si_else"]
            for i in range(self.init_value.z_unit_init * 3 - 2)
        ]

    def deposition(self) -> Tuple:
        dep_pos = deposit_an_atom(
            self.atom_set,
            self.bonds,
            self.init_value,
            self.empty_firstBL,
        )
        self.update_after_deposition(dep_pos)
        return dep_pos

    def update_after_deposition(self, dep_pos) -> None:
        ##
        self.n_events_rec.append(self.n_events_perdep)
        self.num_atoms_rec.append(self.n_atoms)
        ##
        self.n_events_perdep = 0
        self.atom_set[dep_pos] = 1
        self.atom_exist.append(dep_pos)
        self.n_atoms += 1
        self.n_events += 1
        self.prog_time += 1 / (self.init_value.dep_rate_atoms_persec)
        if dep_pos[2] in (0, 1):
            self.empty_firstBL -= 1

    def update_events(self):
        self.related_atoms = list(set(self.related_atoms))
        for target_rel in self.related_atoms:
            if self.atom_set[target_rel] == 0:
                events = []
                rates = []

            else:
                events, rates = site_events(
                    self.atom_set,
                    self.bonds,
                    target_rel,
                    self.init_value,
                    self.energy_bonding,
                    self.energy_diffuse,
                )
            self.total_event_time -= self.event_time_tot[target_rel]
            self.event[target_rel] = events
            self.event_time[target_rel] = rates
            self.event_time_tot[target_rel] = sum(rates)
            # self.event_state[target_rel] = states
            self.total_event_time += self.event_time_tot[target_rel]
        self.related_atoms = []
        #
        #
        # self.time_check()

    def time_check(self):
        tot_time = self.init_value.dep_rate_atoms_persec
        for times in self.event_time_tot.values():
            tot_time += times
        if abs(tot_time - self.total_event_time) == 0:
            pass
        else:
            print("recal_time = " + str(tot_time))
            print("tot_time = " + str(self.total_event_time))
            print("diff = " + str(self.total_event_time - tot_time))
            print("events = " + str(self.num_events))
            input()

    def rejection_free_deposition(self):
        print(
            "Progress: "
            + str(self.n_atoms)
            + "/"
            + str(self.init_value.total_atoms)
            + " Event: "
            + str(self.n_events_perdep)
            + "/"
            + str(self.init_value.cut_number)
        )
        dep_pos = self.deposition()
        # 蒸着によりイベントに変化が生じうる原子
        self.related_atoms = recalculate(
            dep_pos, self.bonds, self.atom_set, self.init_value
        )
        # それぞれのイベント等を格納
        self.update_events()
        #
        #
        # self.prev_dep = dep_pos

    def put_first_atoms_rf(self):
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.put_first)):
                self.rejection_free_deposition()

    def rejection_free_event(self):
        self.move_atom = self.event[self.target][self.event_number]
        # self.new_state = self.event_state[self.target][self.event_number]
        self.event_progress()
        self.related_atoms = recalculate(
            self.target, self.bonds, self.atom_set, self.init_value
        )
        self.related_atoms += recalculate(
            self.move_atom, self.bonds, self.atom_set, self.init_value
        )
        """
        if self.new_state == 4:
            for bond in self.bonds[self.target]:
                self.related_atoms += recalculate(
                    bond, self.bonds, self.atom_set, self.init_value
                )
        """
        self.update_events()
        #
        #
        # self.prev_eve = str(self.target) + ":" + str(self.move_atom)

    def rejection_free_loop(self):
        self.target, self.event_number = rejection_free_choise(
            self.total_event_time, self.event_time, self.event_time_tot
        )
        if self.target == (-1, -1, -1):
            self.rejection_free_deposition()
            # self.prev_eve = "dep_nat"
        elif (self.n_events_perdep == int(self.init_value.cut_number)) and (
            self.init_value.cut_check is True
        ):
            self.rejection_free_deposition()
            # self.prev_eve = "dep_lim"
        else:
            self.rejection_free_event()

        # recoding the positions in the middle
        if self.n_atoms >= self.rec_num_atoms:
            self.rec_num_atoms += self.init_value.rec_num_atom_interval
            self.record_position()

    def record_position(self) -> None:
        self.pos_rec.append(copy.copy(self.atom_set))
        self.time_rec.append(self.prog_time)
        self.cov_rec.append(self.n_atoms / self.init_value.atoms_in_BL)

    def middle_check(self):
        # 構造の途中確認用
        if self.n_atoms == 17 and self.total_event_time < self.min_rates:
            from record_for_test import rec_for_test
            from recording import rec_poscar

            rec_for_test(self.atom_set, self.bonds, self.lattice)
            rec_poscar(
                self.atom_set,
                self.lattice,
                self.init_value.n_cell_init,
                self.init_value.z_unit_init,
                "middle_structure.vasp",
            )

        """
        # 構造の途中確認用2
        if self.n_atoms == 80 and self.record_middle == 0:
            from record_for_test import rec_for_test

            self.record_middle = 1
            rec_for_test(self.atom_set, self.bonds, self.lattice)
        """

    def end_of_loop(self) -> None:
        self.record_position()
        self.elapsed_time = time.time() - self.start_time
        self.minute = math.floor(self.elapsed_time / 60)
        self.second = int(self.elapsed_time % 60)
        self.time_per_event = round(self.elapsed_time / self.n_events * 1000, 3)
        rec_events_per_dep(self.n_events_rec, self.num_atoms_rec, self.init_value)
        self.mode_val = record_data(
            self.pos_rec,
            self.time_rec,
            self.cov_rec,
            self.lattice,
            self.init_value,
            self.minute,
            self.second,
            self.time_per_event,
        )

        #
        # self.time_check()

    def defect_check(self):
        if (self.target[2] not in (0, 1)) and (self.move_atom[2] in (0, 1)):
            self.move_atom = self.target
            # self.new_state = self.atom_set[self.target]
            self.n_events -= 1
            self.n_events_perdep -= 1

    def event_progress(self):
        if (self.empty_firstBL == int(self.init_value.num_defect)) and (
            self.init_value.keep_defect_check is True
        ):
            self.defect_check()
        #
        if self.move_atom == self.target:
            """
            if self.new_state == 4:
                # print("clustering")
                # self.prev_eve = "clustering"
                self.atom_set[self.target] = 3
                for bond in self.bonds[self.target]:
                    if self.atom_set[bond] != 0:
                        self.atom_set[bond] = 3
            elif self.new_state == 5:
                # print("de clustering")
                # self.prev_eve = "change to 2D"
                self.atom_set[self.target] = 2
                for bond in self.bonds[self.target]:
                    if self.atom_set[bond] != 0:
                        self.atom_set[bond] = 2

            else:
            """
            self.atom_set[self.target] = 1
            print("selfmove")

        else:

            self.atom_set[self.move_atom] = 1
            self.atom_set[self.target] = 0
            self.atom_exist.remove(self.target)
            self.atom_exist.append(self.move_atom)
            if self.move_atom[2] in (0, 1):
                self.empty_firstBL -= 1
            if self.target[2] in (0, 1):
                self.empty_firstBL += 1
