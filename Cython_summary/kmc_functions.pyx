# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

import os
from InputParameter cimport Params
import time
from lattice_form cimport (
    lattice_form_lattice, 
    lattice_form_bonds, 
    lattice_form_atom_set, 
    lattice_form_event, 
    lattice_form_event_time, 
    lattice_form_event_time_tot, 
    lattice_form_index_list, 
    lattice_form_diffuse_candidates, 
    lattice_form_highest_atom
)
from deposition cimport deposit_an_atom
from atoms_recalculate cimport recalculate
from rejection_free_choose cimport rejection_free_choise
from event_collection cimport site_events
from Calc_grid_index cimport grid_num
from recording import record_data
from recording import rec_events_per_dep

import pickle

from libcpp.vector cimport vector
from vector_operation cimport search, remove_element, remove_duplicate

cdef class common_functions:

    def __init__(self):
        self.init_value = Params("kmc_input.yml")
        self.unit_length = self.init_value.cell_size_xy
        self.total_time = int(self.init_value.total_time())
        #
        self.z_units = self.init_value.cell_size_z
        self.z_intra = self.init_value.distance_intra
        self.z_inter = self.init_value.distance_inter
        self.prefactor = self.init_value.prefactor
        self.kbt = self.init_value.temperature_eV()
        self.unit_height = self.init_value.unit_height()
        self.z_max = self.init_value.z_max()
        self.num_one_layer = self.init_value.num_one_layer()
        self.num_grids = self.init_value.num_grids()
        self.rec_interval =  self.init_value.rec_num_atom_interval()


    cpdef loop(self):
        while self.prog_time <= self.total_time:
            self.rejection_free_loop()
            self.update_progress()
            #
            # self.middle_check(60)
            # self.isolation_check()
            # self.num_atom_check()
            # self.isolation_full_check()


    cpdef start_setting(self):
        self.start_time = time.time()
        self.prog_time = 0
        self.n_atoms = 0
        self.n_events = 0
        self.pos_rec.clear()
        self.time_rec.clear()
        self.cov_rec.clear()
        self.rec_num_atoms = 0
        self.n_events_perdep = 0
        self.n_events_rec.clear()
        self.num_atoms_rec.clear()
        # self.record_middle = 0
        if os.path.exists("Record") is False:
            os.mkdir("Record")

        """
        self.diff_rec = []
        self.eve_num_rec = []
        """

    cdef update_progress(self):
        self.n_events += 1
        self.n_events_perdep += 1
        if self.n_events % 10000 == 0:
            self.recalculate_total_rate()
        """      
        if self.n_events % 1000 == 0:
            tot_rate = self.compare_rate()
            diff = self.total_event_time - tot_rate
            self.diff_rec.append(diff)
            self.eve_num_rec.append(self.n_events)
        """

    cpdef start_rejection_free(self):
        self.total_event_time = self.init_value.dep_rate_atoms_persec()
        self.lattice = lattice_form_lattice(self.unit_length, self.z_units, self.z_intra, self.z_inter, self.unit_height, self.z_max, self.num_one_layer, self.num_grids)
        self.bonds = lattice_form_bonds(self.unit_length, self.num_grids, self.z_max)
        self.atom_set = lattice_form_atom_set(self.num_grids)
        self.event = lattice_form_event(self.num_grids)
        self.event_time = lattice_form_event_time(self.num_grids)
        self.event_time_tot = lattice_form_event_time_tot(self.num_grids)
        self.index_list = lattice_form_index_list(self.unit_length, self.num_grids, self.z_max)
        self.diffuse_candidates = lattice_form_diffuse_candidates(self.unit_length, self.num_grids, self.z_max, self.bonds, self.index_list)
        self.highest_atom = lattice_form_highest_atom(self.num_one_layer)
        self.energy_summarize()
        self.height_change_rem = False
        self.height_change_add = False
        if self.init_value.trans_check is False:
            self.init_value.trans_num = 100

    cdef energy_summarize(self):
        cdef int _, i
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
        for _ in range(self.init_value.cell_size_z * 6 - 2):
            self.energy_bonding.push_back(self.init_value.energies_binding["upper"])
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
        for _ in range(self.init_value.cell_size_z * 6 - 2):
            self.energy_diffuse.push_back(self.init_value.energies_diffusion["upper"])
            
        for i in range(self.energy_diffuse.size()):
            self.energy_diffuse[i] -= self.energy_bonding[i]

    cdef bint height_check_add(self, int pos):
        cdef int x, y, z, grid_index
        x = self.index_list[pos][0]
        y = self.index_list[pos][1]
        z = self.index_list[pos][2]
        if z >= self.init_value.trans_num:
            grid_index = grid_num(x, y, 0, self.unit_length)
            self.highest_atom[grid_index] += 1
            if self.highest_atom[grid_index] == 1:
                return True
            else:
                return False
        else:
            return False

    cdef bint height_check_remove(self, int pos):
        cdef int x, y, z, grid_index
        x = self.index_list[pos][0]
        y = self.index_list[pos][1]
        z = self.index_list[pos][2]
        if z >= self.init_value.trans_num:
            grid_index = grid_num(x, y, 0, self.unit_length)
            self.highest_atom[grid_index] -= 1
            if self.highest_atom[grid_index] == 0:
                return True
            else:
                return False
        else:
            return False


    cdef update_events(self):
        cdef int target_rel
        cdef vector[int] events
        cdef vector[double] rates
        cdef int list_num 
        cdef double rate, tot
        self.related_atoms = remove_duplicate(self.related_atoms)
        for target_rel in self.related_atoms:
            if self.atom_set[target_rel] == 0:
                events.clear()
                rates.clear()
            else:
                events, rates = site_events(
                    self.atom_set,
                    self.bonds,
                    target_rel,
                    self.unit_length,
                    self.prefactor,
                    self.kbt,
                    self.energy_bonding,
                    self.energy_diffuse,
                    self.diffuse_candidates,
                    self.highest_atom,
                    self.index_list
                )

            self.total_event_time -= self.event_time_tot[target_rel]
            self.event[target_rel] = events
            self.event_time[target_rel] = rates
            tot = 0.0
            for rate in rates:
                tot += rate
            self.event_time_tot[target_rel] = tot

            self.total_event_time += tot

        self.related_atoms.clear()
        """
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        for _, rate in self.event_time_tot.items():
            self.total_event_time += rate
        """

    cdef recalculate_total_rate(self):
        cdef double rate
        self.total_event_time = self.init_value.dep_rate_atoms_persec()
        for rate in self.event_time_tot:
            self.total_event_time += rate

    cdef double compare_rate(self):
        cdef double total_event_time, rate
        total_event_time = self.init_value.dep_rate_atoms_persec()
        for rate in self.event_time_tot:
            total_event_time += rate
        return total_event_time


    cdef update_after_deposition(self, int dep_pos):
        self.n_events_rec.push_back(self.n_events_perdep)
        self.recalculate_total_rate()
        self.n_events_perdep = 0
        self.num_atoms_rec.push_back(self.n_atoms)
        self.n_atoms += 1
        #
        self.atom_set[dep_pos] = 1
        # self.atom_exist.append(dep_pos)
        # self.n_events += 1
        self.prog_time += 1 / (self.init_value.dep_rate_atoms_persec())
        print(self.n_atoms)
        if self.n_atoms >= self.rec_num_atoms:
            self.rec_num_atoms += self.rec_interval
            self.record_position()
            if self.setting_value != -1:
                self.parameter_record()


    cdef int deposition(self):
        cdef int dep_pos
        dep_pos = deposit_an_atom(
            self.atom_set,
            self.bonds,
            self.index_list
        )
        self.update_after_deposition(dep_pos)
        return dep_pos

    cdef rejection_free_deposition(self):
        cdef int dep_pos
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
            self.index_list,
            self.unit_length,
        )  # Picking up atoms that possibly affectsd by the deposition
        #
        self.update_events()

    cpdef put_first_atoms_rf(self):
        cdef int _
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.first_put_num)):
                self.rejection_free_deposition()

    cdef rejection_free_event(self):
        self.move_atom = self.event[self.target][self.event_number]
        self.event_progress()
        self.related_atoms = recalculate(
            self.target,
            self.atom_set,
            self.bonds,
            self.diffuse_candidates,
            self.height_change_rem,
            self.init_value.trans_num,
            self.index_list,
            self.unit_length,
        )
        self.related_atoms_move = recalculate(
            self.move_atom,
            self.atom_set,
            self.bonds,
            self.diffuse_candidates,
            self.height_change_add,
            self.init_value.trans_num,
            self.index_list,
            self.unit_length,
        )

        for index in self.related_atoms_move:
            self.related_atoms.push_back(index)
        self.update_events()


    cdef event_progress(self):
        self.atom_set[self.target] = 0
        self.atom_set[self.move_atom] = 1
        # self.atom_exist.remove(self.target)
        # self.atom_exist.append(self.move_atom)
        #
        self.height_change_rem = self.height_check_remove(self.target)
        self.height_change_add = self.height_check_add(self.move_atom)

    cpdef pickle_dump(self, obj, path):
        with open(path, mode="wb") as f:
            pickle.dump(obj, f)

    cpdef parameter_record(self):
        cdef unicode dir_name, file_name
        dir_name = "Record/" + self.init_value.record_name + "/"
        os.makedirs(dir_name, exist_ok=True)
        file_name = dir_name + "selfdata.pickle"
        self.pickle_dump(self, file_name)

    cdef rejection_free_loop(self):
        self.target, self.event_number = rejection_free_choise(
            self.total_event_time,
            self.event_time,
            self.event_time_tot,
            self.init_value.dep_rate_atoms_persec(),
        )
        if self.target == -1:
            self.rejection_free_deposition()
        else:
            self.rejection_free_event()


    cpdef record_position(self):
        self.pos_rec.push_back(self.atom_set)
        self.time_rec.push_back(self.prog_time)
        self.cov_rec.push_back(self.n_atoms / self.init_value.atoms_in_BL())
        

    cpdef end_of_loop(self):
        cdef double total_time_dir, diff
        self.record_position()
        self.elapsed_time = time.time() - self.start_time
        self.minute = int(self.elapsed_time / 60)
        self.second = int(self.elapsed_time % 60)
        self.time_per_event = round(self.elapsed_time / self.n_events * 1000, 3)
        rec_events_per_dep(self.n_events_rec, self.num_atoms_rec, self.init_value)

        self.mode_val, self.other_modes = record_data(
            self.pos_rec,
            self.time_rec,
            self.cov_rec,
            self.lattice,
            self.init_value.cell_size_z * 6,
            self.unit_length,
            self.init_value.record_name,
            self.init_value.atoms_in_BL(),
            self.init_value.start_from_middle,
            self.z_units,
            self.init_value.temperature,
            self.kbt,
            self.init_value.dep_rate,
            self.init_value.dep_time,
            self.prefactor,
            self.init_value.first_put_check,
            self.init_value.trans_check,
            self.init_value.first_put_num,
            self.init_value.trans_num,
            self.init_value.energies_diffusion,
            self.init_value.energies_binding,
            self.init_value.comments,
            self.minute,
            self.second,
            self.time_per_event,
            self.init_value,
            self.index_list
        )

        total_time_dir = (
            sum(self.event_time_tot) + self.init_value.dep_rate_atoms_persec()
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
