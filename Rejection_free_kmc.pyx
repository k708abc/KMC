#from Modules.kmc_functions cimport common_functions
import datetime
# from Modules.InputParameter cimport Params
import time
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
import pickle

# import decimal
import yaml
from pathlib import Path


cdef class Params:
    cdef int cell_size_xy, cell_size_z, first_put_num, max_workers_val, repeat_val, trans_num
    cdef double dep_rate, dep_time, distance_inter, distance_intra, img_per, prefactor, repeat_E1_diff, repeat_E1_end, repeat_E1_start, repeat_E2_diff, repeat_E2_end, repeat_E2_start, temperature
    cdef list energies_binding, energies_diffusion
    cdef bint first_put_check, start_from_middle, trans_check
    cdef char comments, record_name

    def __cinit__(self, filename) -> None:
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
    def dep_rate_atoms_persec(self) -> float:
        return self.dep_rate / 60 * self.atoms_in_BL

    @property
    def total_time(self) -> float:
        return self.dep_time * 60

    @property
    def interval(self) -> float:
        return self.total_time * self.img_per / 100

    @property
    def total_atoms(self) -> int:
        return int(round(self.dep_rate_atoms_persec * self.total_time))

    @property
    def rec_num_atom_interval(self) -> int:
        return int(round(self.total_atoms * self.img_per / 100))


cdef class common_functions:
    cdef double start_time, total_event_time, elapsed_time, prog_time
    cdef int n_atoms, n_events, rec_num_atoms, n_events_perdep, event_number, setting_value
    cdef list pos_rec, time_rec, cov_rec, n_events_rec, num_atoms_rec, event_time_tot, energy_bonding, energy_diffuse, related_atoms
    cdef dict lattice, bonds, atom_set, event, event_time, site_list_correspondance, diffuse_candidates,highest_atom
    cdef bint height_change_rem, height_change_add
    cdef tuple move_atom, target
    cdef Params init_value

    def __cinit__(self):
        self.init_value = Params("kmc_input.yml")

    cpdef start_setting(self):
        self.start_time = time.time()
        self.prog_time = 0
        self.n_atoms = 0
        self.n_events = 0
        self.pos_rec = []
        self.time_rec = []
        self.cov_rec = []
        self.rec_num_atoms = 0
        self.n_events_perdep = 0
        self.n_events_rec = []
        self.num_atoms_rec = []
        # self.record_middle = 0
        if os.path.exists("Record") is False:
            os.mkdir("Record")

        """
        self.diff_rec = []
        self.eve_num_rec = []
        """

    cpdef update_progress(self):
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

    cpdef start_rejection_free(self):
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

    cpdef energy_summarize(self):
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

    cpdef height_check_add(self, tuple pos):
        if pos[2] >= self.init_value.trans_num:
            self.highest_atom[(pos[0], pos[1])] += 1
            if self.highest_atom[(pos[0], pos[1])] == 1:
                return True
            else:
                return False
        else:
            return False

    cpdef bint height_check_remove(self, tuple pos):
        if pos[2] >= self.init_value.trans_num:
            self.highest_atom[(pos[0], pos[1])] -= 1
            if self.highest_atom[(pos[0], pos[1])] == 0:
                return True
            else:
                return False
        else:
            return False

    cpdef tuple deposition(self):
        cdef tuple dep_pos
        dep_pos = deposit_an_atom(
            self.atom_set,
            self.bonds,
        )
        self.update_after_deposition(dep_pos)
        return dep_pos

    cpdef update_after_deposition(self, tuple dep_pos):
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

    cpdef update_events(self):
        cdef tuple target_rel
        cdef list events, rates
        cdef int list_num 
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

    cpdef recalculate_total_rate(self):
        cdef double rate
        self.total_event_time = self.init_value.dep_rate_atoms_persec
        for rate in self.event_time_tot:
            self.total_event_time += rate

    cpdef double compare_rate(self):
        cdef double total_event_time, rate
        total_event_time = self.init_value.dep_rate_atoms_persec
        for rate in self.event_time_tot:
            total_event_time += rate
        return total_event_time

    cpdef rejection_free_deposition(self):
        cdef tuple dep_pos
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

    cpdef put_first_atoms_rf(self):
        cdef int _
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.first_put_num)):
                self.rejection_free_deposition()

    cpdef rejection_free_event(self):
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

    cpdef event_progress(self):
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

    cpdef rejection_free_loop(self):
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

    cpdef record_position(self):
        self.pos_rec.append(copy.copy(self.atom_set))
        self.time_rec.append(self.prog_time)
        self.cov_rec.append(self.n_atoms / self.init_value.atoms_in_BL)

    cpdef end_of_loop(self):
        cdef double total_time_dir, diff
        self.record_position()
        self.elapsed_time = time.time() - self.start_time
        self.minute = math.floor(self.elapsed_time / 60)
        self.second = int(self.elapsed_time % 60)
        self.time_per_event = round(self.elapsed_time / self.n_events * 1000, 3)
        rec_events_per_dep(self.n_events_rec, self.num_atoms_rec, self.init_value)

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

        total_time_dir = (
            sum(self.event_time_tot) + self.init_value.dep_rate_atoms_persec
        )
        diff = self.total_event_time - total_time_dir
        print("diff_time = " + str(diff))


cdef class rejection_free(common_functions):

    def __cinit__(self, int setting_value):
        common_functions.__init__(self)
        self.setting_value = setting_value
        """
        setting_value: executed by
            0   Rejection_free_kmc.py
            1   Rejection_free_repetetion.py
            2   RFKMC_multi.py
        """

    cpdef loop(self):
        while int(self.prog_time) <= int(self.init_value.total_time):
            self.rejection_free_loop()
            self.update_progress()
            #
            # self.middle_check(60)
            # self.isolation_check()
            # self.num_atom_check()
            # self.isolation_full_check()

    cpdef end(self):
        if self.setting_value in (1, 0):
            print("Recording")
        self.end_of_loop()
        # self.trans_check()
        # self.time_check()
        if self.setting_value in (1, 0):
            print("Finished: " + str(self.minute) + " min " + str(self.second) + " sec")
            print("Time/event: " + str(self.time_per_event) + " ms")
        if self.setting_value == 0:
            pass

        """
        fig = plt.figure()
        plt.plot(self.eve_num_rec, self.diff_rec)
        plt.show()
        """

    cpdef start(self):
        if self.setting_value in (1, 0):
            print("Calculation start")
        self.start_setting()
        self.start_rejection_free()
        # Put first atom
        self.put_first_atoms_rf()
        if self.setting_value in (1, 0):
            print("Loop start")
        self.loop()
        self.end()

    cpdef start_from_middle(self):
        print("current number of atoms" + str(self.n_atoms))
        print("started time" + str(datetime.datetime.fromtimestamp(self.start_time)))
        self.loop()
        self.end()


if __name__ == "__main__":
    rf_class = rejection_free(0)
    rf_class.start()
