from InputParameter import Params
import time
from typing import List, Dict, Tuple
from Modules.lattice_form import lattice_form
from Modules.deposition import deposit_an_atom
from Modules.atoms_recalculate import recalculate
from Modules.rejection_free_choose import rejection_free_choise
from Modules.event_collection import site_events
import math
import copy
from Modules.recording import record_data
from Modules.recording import rec_events_per_dep
import os


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
        if os.path.exists("Record") is False:
            os.mkdir("Record")

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
            self.diffuse_candidates,
            self.highest_atom,
        ) = lattice_form(self.init_value)
        self.energy_summarize()
        self.height_change_rem = (0, 0)
        self.height_change_add = (0, 0)
        if self.init_value.transformation is False:
            self.init_value.trans_val = 100

    def energy_summarize(self):
        self.energy_bonding = [
            self.init_value.binding_energies["Si_first"],  # with Ag
            self.init_value.binding_energies["Si_first"],  # in first BL
            (
                self.init_value.binding_energies["Si_first"]
                + self.init_value.binding_energies["Si_second"]
            )
            / 2,  # inter 1st and 2nd
            self.init_value.binding_energies["Si_second"],  # in second
            (
                self.init_value.binding_energies["Si_second"]
                + self.init_value.binding_energies["Si_third"]
            )
            / 2,  # inter second and third
            self.init_value.binding_energies["Si_third"],  # in third
            (
                self.init_value.binding_energies["Si_third"]
                + self.init_value.binding_energies["Si_else"]
            )
            / 2,  # inter third and fourth
        ]
        #
        self.energy_bonding += [
            self.init_value.binding_energies["Si_else"]
            for i in range(self.init_value.z_unit_init * 6 - 2)
        ]
        #

        self.energy_diffuse = [
            self.init_value.diffusion_barriers["Si_first"],
            self.init_value.diffusion_barriers["Si_first"],
            self.init_value.diffusion_barriers["Si_second"],
            self.init_value.diffusion_barriers["Si_second"],
            self.init_value.diffusion_barriers["Si_third"],
            self.init_value.diffusion_barriers["Si_third"],
        ]
        #
        self.energy_diffuse += [
            self.init_value.diffusion_barriers["Si_else"]
            for i in range(self.init_value.z_unit_init * 6 - 2)
        ]
        if self.init_value.subtract_check is True:
            for i in range(len(self.energy_diffuse)):
                self.energy_diffuse[i] -= self.energy_bonding[i]

    def height_check_add(self, pos):
        prev = self.highest_atom[(pos[0], pos[1])]
        if prev < pos[2]:
            self.highest_atom[(pos[0], pos[1])] = pos[2]
            return (prev, pos[2])
        else:
            return (prev, prev)

    def height_check_remove(self, pos):
        prev = self.highest_atom[(pos[0], pos[1])]
        if prev == pos[2]:
            for i in range(pos[2], 0, -1):
                if self.atom_set[(pos[0], pos[1], i)] == 1:
                    self.highest_atom[(pos[0], pos[1])] = i
                    return (i, prev)
            self.highest_atom[(pos[0], pos[1])] = 0
            return (0, prev)
        else:
            return (prev, prev)

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
                    self.diffuse_candidates,
                    self.highest_atom,
                    self.init_value.trans_val,
                )
            self.total_event_time -= self.event_time_tot[target_rel]
            self.event[target_rel] = events
            self.event_time[target_rel] = rates
            self.event_time_tot[target_rel] = sum(rates)
            self.total_event_time += self.event_time_tot[target_rel]
        self.related_atoms = []

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
        if self.num in (1, 0):
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
        if self.init_value.transformation is True:
            self.height_change_add = self.height_check_add(dep_pos)
        # 蒸着によりイベントに変化が生じうる原子
        self.related_atoms = recalculate(
            dep_pos,
            self.atom_set,
            self.diffuse_candidates,
            self.height_change_add,
            self.init_value.trans_val,
        )
        # それぞれのイベント等を格納
        self.update_events()

    def put_first_atoms_rf(self):
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.put_first)):
                self.rejection_free_deposition()

    def rejection_free_event(self):
        self.move_atom = self.event[self.target][self.event_number]
        self.event_progress()
        self.related_atoms = recalculate(
            self.target,
            self.atom_set,
            self.diffuse_candidates,
            self.height_change_rem,
            self.init_value.trans_val,
        )
        self.related_atoms += recalculate(
            self.move_atom,
            self.atom_set,
            self.diffuse_candidates,
            self.height_change_add,
            self.init_value.trans_val,
        )
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
        elif (self.n_events_perdep == int(self.init_value.cut_number)) and (
            self.init_value.cut_check is True
        ):
            self.rejection_free_deposition()
        else:
            self.rejection_free_event()

        # recoding the positions in the middle
        if self.n_atoms >= self.rec_num_atoms:
            self.rec_num_atoms += self.init_value.rec_num_atom_interval
            self.record_position()
        # self.middle_check()

    def record_position(self) -> None:
        self.pos_rec.append(copy.copy(self.atom_set))
        self.time_rec.append(self.prog_time)
        self.cov_rec.append(self.n_atoms / self.init_value.atoms_in_BL)

    def middle_check(self):
        # 構造の途中確認用
        if self.n_atoms == 30 and self.total_event_time < self.min_rates:
            from Test_modules.record_for_test import rec_for_test
            from Modules.recording import rec_poscar

            rec_for_test(
                self.atom_set, self.bonds, self.lattice, self.diffuse_candidates
            )
            rec_poscar(
                self.atom_set,
                self.lattice,
                self.init_value.n_cell_init,
                self.init_value.z_unit_init,
                "middle_structure.vasp",
            )
            print("middle formed")
            input()

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

    def defect_check(self):
        if (self.target[2] not in (0, 1)) and (self.move_atom[2] in (0, 1)):
            self.move_atom = self.target
            self.n_events -= 1
            self.n_events_perdep -= 1

    def event_progress(self):
        if (self.empty_firstBL == int(self.init_value.num_defect)) and (
            self.init_value.keep_defect_check is True
        ):
            self.defect_check()
        #
        self.atom_set[self.target] = 0
        self.atom_set[self.move_atom] = 1
        self.atom_exist.remove(self.target)
        self.atom_exist.append(self.move_atom)
        if self.init_value.transformation is True:
            self.height_change_rem = self.height_check_remove(self.target)
            self.height_change_add = self.height_check_add(self.move_atom)
        if self.move_atom[2] in (0, 1):
            self.empty_firstBL -= 1
        if self.target[2] in (0, 1):
            self.empty_firstBL += 1
