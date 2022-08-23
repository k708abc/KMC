from Modules.InputParameter import Params
import time
from typing import List, Tuple
from Modules.lattice_form import lattice_form
from Modules.deposition import deposit_an_atom
from Modules.atoms_recalculate import recalculate
from Modules.rejection_free_choose import rejection_free_choise
from Modules.event_collection import site_events
import math
import copy

# from Modules.recording import record_data
# from Modules.recording import rec_events_per_dep
import os
import pickle


class common_functions:
    def __init__(self) -> None:
        self.init_value = Params("kmc_input.yml")

    def start_setting(self):
        self.start_time = time.time()
        self.prog_time = 0
        self.n_atoms = 0
        self.n_events = 0
        self.pos_rec: List[dict] = []
        self.time_rec: List[float] = []
        self.cov_rec: List[float] = []
        self.rec_num_atoms = 0
        self.n_events_perdep = 0
        self.n_events_rec: List[int] = []
        self.num_atoms_rec: List[int] = []
        # self.record_middle = 0
        if os.path.exists("Record") is False:
            os.mkdir("Record")

        """
        self.diff_rec = []
        self.eve_num_rec = []
        """

    def update_progress(self) -> None:
        self.n_events += 1
        self.n_events_perdep += 1

        if self.n_events % 100000 == 0:
            self.recalculate_total_rate()
        """      
        if self.n_events % 1000 == 0:
            tot_rate = self.compare_rate()
            diff = self.total_event_time - tot_rate
            self.diff_rec.append(diff)
            self.eve_num_rec.append(self.n_events)
        """

    def start_rejection_free(self):
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        # self.atom_exist: List[Tuple[int, int, int]] = []
        (
            self.lattice,
            self.bonds,
            self.atom_set,
            self.event,
            self.event_time,
            self.event_time_tot,
            self.site_list_correspondance,
            self.list_site_correspondance,
            self.diffuse_candidates,
            self.highest_atom,
        ) = lattice_form(self.init_value)
        self.energy_summarize()
        self.height_change_rem = False
        self.height_change_add = False
        if self.init_value.trans_check is False:
            self.init_value.trans_num = 100

    def energy_summarize(self):
        # interlayer bonding is the average of the intralayer bonding
        self.energy_bonding = [
            self.init_value.energies_binding["first"],  # with Ag
            self.init_value.energies_binding["first"],  # intra first BL
            (
                self.init_value.energies_binding["first"]
                + self.init_value.energies_binding["second"]
            )
            / 2,  # inter 1st and 2nd
            self.init_value.energies_binding["second"],  # intra second
            (
                self.init_value.energies_binding["second"]
                + self.init_value.energies_binding["third"]
            )
            / 2,  # inter second and third
            self.init_value.energies_binding["third"],  # intra third
            (
                self.init_value.energies_binding["third"]
                + self.init_value.energies_binding["upper"]
            )
            / 2,  # inter third and fourth
        ]
        #
        self.energy_bonding += [
            self.init_value.energies_binding["upper"]
            for _ in range(self.init_value.cell_size_z * 6 - 2)
        ]
        #

        self.energy_diffuse = [
            self.init_value.energies_diffusion["first"],
            self.init_value.energies_diffusion["first"],
            self.init_value.energies_diffusion["second"],
            self.init_value.energies_diffusion["second"],
            self.init_value.energies_diffusion["third"],
            self.init_value.energies_diffusion["third"],
        ]
        #
        self.energy_diffuse += [
            self.init_value.energies_diffusion["upper"]
            for _ in range(self.init_value.cell_size_z * 6 - 2)
        ]

        for i in range(len(self.energy_diffuse)):
            self.energy_diffuse[i] -= self.energy_bonding[i]

    def height_check_add(self, pos):
        if pos[2] >= self.init_value.trans_num:
            self.highest_atom[(pos[0], pos[1])] += 1
            if self.highest_atom[(pos[0], pos[1])] == 1:
                return True
            else:
                return False
        else:
            return False

    def height_check_remove(self, pos):
        if pos[2] >= self.init_value.trans_num:
            self.highest_atom[(pos[0], pos[1])] -= 1
            if self.highest_atom[(pos[0], pos[1])] == 0:
                return True
            else:
                return False
        else:
            return False

    def deposition(self) -> Tuple:
        dep_pos = deposit_an_atom(
            self.atom_set,
            self.bonds,
        )
        self.update_after_deposition(dep_pos)
        return dep_pos

    def update_after_deposition(self, dep_pos) -> None:
        #
        self.n_events_rec.append(self.n_events_perdep)
        self.n_events_perdep = 0
        self.num_atoms_rec.append(self.n_atoms)
        self.n_atoms += 1
        #
        self.atom_set[dep_pos] = 1
        # self.atom_exist.append(dep_pos)

        # self.n_events += 1
        self.prog_time += 1 / (self.init_value.dep_rate_atoms_persec)
        print(self.n_atoms)

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
                )

            list_num = self.site_list_correspondance[target_rel]
            self.total_event_time -= self.event_time_tot[list_num]
            self.event[target_rel] = events
            self.event_time[target_rel] = rates
            self.event_time_tot[list_num] = sum(rates)
            self.total_event_time += self.event_time_tot[list_num]
        self.related_atoms = []
        """
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        for _, rate in self.event_time_tot.items():
            self.total_event_time += rate
        """

    def recalculate_total_rate(self):
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        for rate in self.event_time_tot:
            self.total_event_time += rate

    def compare_rate(self):
        total_event_time = self.init_value.dep_rate_atoms_persec
        for rate in self.event_time_tot:
            total_event_time += rate
        return total_event_time

    def rejection_free_deposition(self):
        dep_pos = self.deposition()
        self.move_atom = dep_pos
        self.target = dep_pos
        self.height_change_add = self.height_check_add(dep_pos)
        #
        self.related_atoms = recalculate(
            dep_pos,
            self.atom_set,
            self.bonds,
            self.diffuse_candidates,
            self.height_change_add,
            self.init_value.trans_num,
        )  # Picking up atoms that possibly affectsd by the deposition
        #
        self.update_events()

    def put_first_atoms_rf(self):
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.first_put_num)):
                self.rejection_free_deposition()

    def rejection_free_event(self):
        self.move_atom = self.event[self.target][self.event_number]
        self.event_progress()
        self.related_atoms = recalculate(
            self.target,
            self.atom_set,
            self.bonds,
            self.diffuse_candidates,
            self.height_change_rem,
            self.init_value.trans_num,
        )
        self.related_atoms += recalculate(
            self.move_atom,
            self.atom_set,
            self.bonds,
            self.diffuse_candidates,
            self.height_change_add,
            self.init_value.trans_num,
        )
        self.update_events()

    def event_progress(self):
        self.atom_set[self.target] = 0
        self.atom_set[self.move_atom] = 1
        # self.atom_exist.remove(self.target)
        # self.atom_exist.append(self.move_atom)
        #
        self.height_change_rem = self.height_check_remove(self.target)
        self.height_change_add = self.height_check_add(self.move_atom)

    def pickle_dump(self, obj, path):
        with open(path, mode="wb") as f:
            pickle.dump(obj, f)

    def parameter_record(self):
        dir_name = "Record/" + self.init_value.record_name + "/"
        os.makedirs(dir_name, exist_ok=True)
        file_name = dir_name + "selfdata.pickle"
        self.pickle_dump(self, file_name)

    def rejection_free_loop(self):
        self.target, self.event_number = rejection_free_choise(
            self.total_event_time,
            self.event_time,
            self.event_time_tot,
            self.list_site_correspondance,
            self.init_value.dep_rate_atoms_persec,
        )
        if self.target == (-1, -1, -1):
            self.rejection_free_deposition()
        else:
            self.rejection_free_event()

        if self.n_atoms >= self.rec_num_atoms:
            self.rec_num_atoms += self.init_value.rec_num_atom_interval
            self.record_position()
            if self.setting_value != -1:
                self.parameter_record()

    def record_position(self) -> None:
        self.pos_rec.append(copy.copy(self.atom_set))
        self.time_rec.append(self.prog_time)
        self.cov_rec.append(self.n_atoms / self.init_value.atoms_in_BL)

    def end_of_loop(self) -> None:
        # self.record_position()
        self.elapsed_time = time.time() - self.start_time
        self.minute = math.floor(self.elapsed_time / 60)
        self.second = int(self.elapsed_time % 60)
        self.time_per_event = round(self.elapsed_time / self.n_events * 1000, 3)
        # rec_events_per_dep(self.n_events_rec, self.num_atoms_rec, self.init_value)
        """
        self.mode_val, self.other_modes = record_data(
            self.pos_rec,
            self.time_rec,
            self.cov_rec,
            self.lattice,
            self.init_value,
            self.minute,
            self.second,
            self.time_per_event,
            self.init_value,
        )
        """
        total_time_dir = (
            sum(self.event_time_tot) + self.init_value.dep_rate_atoms_persec
        )
        diff = self.total_event_time - total_time_dir
        print("diff_time = " + str(diff))

    """
    def middle_check(self, val):
        # 構造の途中確認用
        if self.n_atoms == val and self.record_middle == 0:
            from Test_modules.Test_conditions.record_for_test import rec_for_test
            from Modules.recording import rec_poscar

            path = "Test_modules/Test_conditions/"
            rec_for_test(
                self.atom_set, self.bonds, self.lattice, self.diffuse_candidates, path
            )
            rec_poscar(
                self.atom_set,
                self.lattice,
                self.init_value.cell_size_xy,
                self.init_value.cell_size_z,
                "Test_modules/Test_conditions/middle_structure.vasp",
            )
            print("middle formed")
            self.record_middle = 1
            input()

    def trans_check(self):
        for i in range(0, self.init_value.cell_size_xy):
            for k in range(0, self.init_value.cell_size_xy):
                above_trans = 0
                for z in range(
                    self.init_value.trans_num, self.init_value.cell_size_z * 3 - 2
                ):
                    if self.atom_set[(i, k, z)] == 1:
                        above_trans += 1
                if above_trans == self.highest_atom[(i, k)]:
                    pass
                else:
                    print("Transformation counts different")
        print("Trans check finished")

    def isolation_check(self):
        num_bonds = 0
        for i in self.bonds[self.move_atom]:
            if self.atom_set[i] == 1:
                num_bonds += 1
        if num_bonds == 0 and self.move_atom[2] != 0:
            print("Isolate after move")
            print("Before move: " + str(self.target))
            print("After move: " + str(self.move_atom))
            print("Existing atoms: " + str(self.atom_exist))
            input()

        for i in self.bonds[self.target]:
            if self.atom_set[i] == 1:
                bonding = []
                for k in self.bonds[i]:
                    if self.atom_set[k] == 1:
                        bonding.append(k)
                if i[2] == 0:
                    pass
                elif len(bonding) == 0:
                    print("Surrounding atom isolates")
                    print("Position" + str(i))
                    print("Atom befor move: " + str(self.target))
                    print("Atom after move: " + str(self.move_atom))
                    input()
                elif (
                    (len(bonding) == 1)
                    and (bonding[0] == self.move_atom)
                    and (num_bonds == 1)
                    and (self.move_atom[2] != 0)
                ):
                    print("Dimer isolates")
                    print("Position" + str(i))
                    print("Atom befor move: " + str(self.target))
                    print("Atom after move: " + str(self.move_atom))
                    input()


    def isolation_full_check(self):
        checked_atom = []
        clustering_atoms = []
        checking_atom = []
        num_cluster = 0

        for pos in self.atom_exist:
            if pos in checked_atom:
                pass
            else:
                check = False
                clustering_atoms = [pos]
                next_check = [pos]

                while check is False:
                    check = True
                    checking_atom = next_check.copy()
                    next_check = []
                    for atom in checking_atom:
                        for bond in self.bonds[atom]:
                            if (self.atom_set[bond] == 1) and (
                                bond not in clustering_atoms
                            ):
                                next_check.append(bond)
                                clustering_atoms.append(bond)
                                check = False
                        checked_atom.append(atom)
                #
                connected_Ag = 0
                for atom in clustering_atoms:
                    if atom[2] == 0:
                        connected_Ag += 1
                if connected_Ag == 0:
                    print("floating")
                    print("Atom befor move: " + str(self.target))
                    print("Atom after move: " + str(self.move_atom))
                    print("clustering:" + str(clustering_atoms))
                    from Test_modules.Test_conditions.record_for_test import (
                        rec_for_test,
                    )
                    from Modules.recording import rec_poscar
                    import sys

                    path = "Test_modules/middle_rec/"
                    if os.path.exists(path) is False:
                        os.mkdir(path)
                    rec_for_test(
                        self.atom_set,
                        self.bonds,
                        self.lattice,
                        self.diffuse_candidates,
                        path,
                    )
                    rec_poscar(
                        self.atom_set,
                        self.lattice,
                        self.init_value.cell_size_xy,
                        self.init_value.cell_size_z,
                        "Test_modules/middle_rec/middle_structure_"
                        + str(num_cluster)
                        + ".vasp",
                    )
                    num_cluster += 1
                    #
                    #
                    self.atom_set[self.move_atom] = 0
                    self.atom_set[self.target] = 1
                    #
                    events, rates = site_events(
                        self.atom_set,
                        self.bonds,
                        self.target,
                        self.init_value,
                        self.energy_bonding,
                        self.energy_diffuse,
                        self.diffuse_candidates,
                        self.highest_atom,
                    )
                    print("possible diffusion: " + str(events))
                    if self.move_atom in events:
                        print("isolation function fail")
                    else:
                        print("isolation functionis OK")
                    self.atom_set[self.target] = 0
                    self.atom_set[self.move_atom] = 1

        if num_cluster != 0:
            sys.exit()

    def num_atom_check(self):
        num_atom = 0
        for atom_pos, state in self.atom_set.items():
            if state == 1:
                num_atom += 1
        if num_atom != self.n_atoms:
            print("Current number of atoms: " + str(num_atom))
            print("Counted number of atoms: " + str(self.n_atoms))

    def time_check(self):
        tot_time = self.init_value.dep_rate_atoms_persec
        for times in self.event_time_tot.values():
            tot_time += times
        if abs(tot_time - self.total_event_time) < 0.0000001:
            print("Time check OK")
        else:
            print("Time check differt:")
            print("recal_time = " + str(tot_time))
            print("tot_time = " + str(self.total_event_time))
            print("diff = " + str(self.total_event_time - tot_time))
            input()
    """
