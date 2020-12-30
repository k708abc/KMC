#!/usr/bin/env python3

from typing import List, Dict, Tuple
import tkinter as tk
import tkinter.ttk as ttk

# import copy
import math
from cal_rates import rate
import time
from lattice_form import lattice_form

# from InputParameter import Params
# from deposition import deposit_an_atom
from choose_site import choose_atom
from judgement import judge_null
from event_collection import site_events
from normarize_list import normarize_rate
from weighted_choice import choice

# from recording import record_data

# from atoms_recalculate import recalculate
# from rejection_free_choose import rejection_free_choise
from kmc_functions import common_functions

##
# from recording import rec_events_per_dep


class Window(ttk.Frame, common_functions):
    kb_eV = 8.617e-5
    padWE: Dict = dict(sticky=(tk.W, tk.E), padx=15, pady=2)

    def __init__(self, master) -> None:
        super().__init__(master, padding=2)
        # self.init_value = Params()
        common_functions.__init__(self)
        self.create_frame_basics()
        self.create_frame_energies()
        self.create_frame_checks()
        self.create_frame_memos()
        self.create_method()
        self.create_frame_buttons()
        self.create_frame_progress()
        self.create_frame_bar()
        master.title("KMC_Si")
        self.update_values()

    def create_frame_basics(self) -> None:
        self.frame_basics = ttk.Frame()
        self.create_widgets_basics()
        self.create_layout_basic()
        self.frame_basics.pack()

    def create_frame_energies(self) -> None:
        self.frame_energies = ttk.Frame()
        self.create_widgets_energies()
        self.create_layout_energies()
        self.frame_energies.pack()

    def create_frame_checks(self) -> None:
        self.frame_checks = ttk.Frame()
        self.create_widgets_checks()
        self.create_layout_checks()
        self.frame_checks.pack()

    def create_frame_memos(self) -> None:
        self.frame_records = ttk.Frame()
        self.create_widgets_records()
        self.create_layout_records()
        self.frame_records.pack()

    def create_method(self) -> None:
        self.frame_method = ttk.Frame()
        self.create_widgets_method()
        self.create_layout_method()
        self.frame_method.pack()

    def create_frame_buttons(self) -> None:
        self.frame_buttons = ttk.Frame()
        self.create_widgets_buttons()
        self.create_layout_buttons()
        self.frame_buttons.pack()

    def create_frame_bar(self) -> None:
        self.frame_bars = ttk.Frame()
        self.create_widgets_bar()
        self.create_layout_bar()
        self.frame_bars.pack()

    def create_frame_progress(self) -> None:
        self.frame_progress = ttk.Frame()
        self.create_widgets_progress()
        self.create_layout_progress()
        self.frame_progress.pack()

    def create_widgets_basics(self) -> None:
        # The first row
        self.n_cell_label = ttk.Label(self.frame_basics, text="Number of cell")
        self.n_cell = ttk.Entry(self.frame_basics, width=7)
        self.n_cell.insert(tk.END, self.init_value.n_cell_init)
        self.n_cell.bind("<Return>", self.update_click)
        self.z_unit_label = ttk.Label(self.frame_basics, text="Z unit")
        self.z_unit = ttk.Entry(self.frame_basics, width=7)
        self.z_unit.insert(tk.END, self.init_value.z_unit_init)
        self.z_unit.bind("<Return>", self.update_click)
        self.temperature_label = ttk.Label(self.frame_basics, text="T (K)")
        self.temperature = ttk.Entry(self.frame_basics, width=7)
        self.temperature.insert(tk.END, self.init_value.temperature)
        self.temperature.bind("<Return>", self.update_click)
        self.temperautre_energy_label = ttk.Label(self.frame_basics, text="kbT")
        self.temperautre_energy = ttk.Label(
            self.frame_basics, text="{:.3g}".format(self.init_value.temperature_eV)
        )
        # The second row
        self.deposition_rate_label = ttk.Label(
            self.frame_basics, text="Dep. rate(ML/min)"
        )
        self.deposition_rate = ttk.Entry(self.frame_basics, width=7)
        self.deposition_rate.insert(tk.END, self.init_value.dep_rate)
        self.deposition_rate.bind("<Return>", self.update_click)

        self.dep_rate_conv_label = ttk.Label(
            self.frame_basics, text="Dep. rate(atoms/s)"
        )
        self.dep_rate_conv_val = ttk.Label(
            self.frame_basics,
            text=("{:.3f}".format(self.init_value.dep_rate_atoms_persec)),
        )

        self.deposition_time_label = ttk.Label(
            self.frame_basics, text="Dep. time (min)"
        )
        self.deposition_time = ttk.Entry(self.frame_basics, width=7)
        self.deposition_time.insert(tk.END, self.init_value.dep_time)
        self.deposition_time.bind("<Return>", self.update_click)
        self.prefactor_label = ttk.Label(self.frame_basics, text="Prefactor (1/s)")
        self.prefactor = ttk.Entry(self.frame_basics, width=7)
        self.prefactor.insert(tk.END, self.init_value.prefactor)
        self.prefactor.bind("<Return>", self.update_click)

    def create_layout_basic(self) -> None:
        """Frame for number of cell, Z unit, ..."""
        # the first row
        self.n_cell_label.grid(row=0, column=0, **self.padWE)
        self.n_cell.grid(row=0, column=1, **self.padWE)
        self.z_unit_label.grid(row=0, column=2, **self.padWE)
        self.z_unit.grid(row=0, column=3, **self.padWE)
        self.temperature_label.grid(row=0, column=4, **self.padWE)
        self.temperature.grid(row=0, column=5, **self.padWE)
        self.temperautre_energy_label.grid(row=0, column=6, **self.padWE)
        self.temperautre_energy.grid(row=0, column=7, **self.padWE)
        # the second row
        self.deposition_rate_label.grid(row=1, column=0, **self.padWE)
        self.deposition_rate.grid(row=1, column=1, **self.padWE)

        self.dep_rate_conv_label.grid(row=1, column=2, **self.padWE)
        self.dep_rate_conv_val.grid(row=1, column=3, **self.padWE)

        self.deposition_time_label.grid(row=1, column=4, **self.padWE)
        self.deposition_time.grid(row=1, column=5, **self.padWE)

        self.prefactor_label.grid(row=1, column=6, **self.padWE)
        self.prefactor.grid(row=1, column=7, **self.padWE)

    def create_widgets_energies(self) -> None:
        self.labels: List[str] = [
            "  ",
            "Base",
            "AgSi",
            "Si12",
            "Si23",
            "Si34",
            "Si45",
            "Si56",
            "Si_intra",
            "Si_inter",
        ]

        self.energylabels = []
        self.energies = []
        self.rate_labels = []
        for _ in range(9):
            self.energies.append(ttk.Entry(self.frame_energies, width=7))
            self.rate_labels.append(ttk.Label(self.frame_energies, text="0"))
        for (energy, (key, val)) in zip(
            self.energies, self.init_value.binding_energies.items()
        ):
            energy.insert(tk.END, val)
            energy.bind("<Return>", self.update_click)
        for label in self.labels:
            self.energylabels.append(ttk.Label(self.frame_energies, text=label))

    def create_layout_energies(self) -> None:
        # self.update_values()
        for i, energylabel in enumerate(self.energylabels):
            energylabel.grid(row=0, column=i, **self.padWE)
        self.energy_label0 = ttk.Label(self.frame_energies, text="Energy (eV)")
        self.energy_label0.grid(row=1, column=0, **self.padWE)
        for i, energy_entry in enumerate(self.energies):
            energy_entry.grid(row=1, column=i + 1)
        self.rate_bond_label = ttk.Label(self.frame_energies, text="Rates/bond")
        self.rate_bond_label.grid(row=2, **self.padWE)
        for i, ratelabel in enumerate(self.rate_labels):
            ratelabel.grid(row=2, column=i + 1)

    def create_widgets_checks(self) -> None:
        self.bln_tr = tk.BooleanVar()
        self.bln_tr.set(self.init_value.trans_check)
        self.chk_transformation_label = ttk.Label(
            self.frame_checks, text="Transformation"
        )
        self.chk_transformation = ttk.Checkbutton(
            self.frame_checks, variable=self.bln_tr
        )
        self.transformation = ttk.Entry(self.frame_checks, width=7)
        self.transformation.insert(tk.END, self.init_value.transformation)
        self.chk_transformation_label2 = ttk.Label(self.frame_checks, text="")
        #
        self.bln_defect = tk.BooleanVar()
        self.bln_defect.set(self.init_value.keep_defect_check)
        self.chk_defect_label = ttk.Label(
            self.frame_checks, text="Keep defect in the first layer"
        )
        self.chk_defect = ttk.Checkbutton(self.frame_checks, variable=self.bln_defect)
        self.num_defect = ttk.Entry(self.frame_checks, width=7)
        self.num_defect.insert(tk.END, self.init_value.num_defect)
        self.chk_defect_label2 = ttk.Label(self.frame_checks, text="")
        #
        self.bln_first = tk.BooleanVar()
        self.bln_first.set(self.init_value.first_put_check)
        self.chk_first_put_label = ttk.Label(
            self.frame_checks, text="Put atom at first"
        )
        self.chk_first_put = ttk.Checkbutton(self.frame_checks, variable=self.bln_first)
        self.put_first = ttk.Entry(self.frame_checks, width=7)
        self.put_first.insert(tk.END, self.init_value.put_first)
        self.chk_first_put_label2 = ttk.Label(self.frame_checks, text="")
        #
        self.bln_cut = tk.BooleanVar()
        self.bln_cut.set(self.init_value.cut_check)
        self.chk_cut_label = ttk.Label(self.frame_checks, text="Cut event")
        self.chk_cut = ttk.Checkbutton(self.frame_checks, variable=self.bln_cut)
        self.cut = ttk.Entry(self.frame_checks, width=7)
        self.cut.insert(tk.END, self.init_value.cut_number)
        self.chk_cut_label2 = ttk.Label(self.frame_checks, text="")
        #
        self.bln_limit = tk.BooleanVar()
        self.bln_limit.set(self.init_value.limit_check)
        self.chk_limit_label = ttk.Label(self.frame_checks, text="Rate limit")
        self.chk_limit = ttk.Checkbutton(self.frame_checks, variable=self.bln_limit)
        self.limit = ttk.Entry(self.frame_checks, width=7)
        self.limit.insert(tk.END, self.init_value.limit_val)
        self.limit.bind("<Return>", self.update_click)
        self.chk_limit_label2 = ttk.Label(
            self.frame_checks,
            text="(" + "{:.3g}".format(self.init_value.cal_energy) + " eV)",
        )

    def create_layout_checks(self) -> None:
        self.chk_transformation_label.grid(row=0, column=0, **self.padWE)
        self.chk_transformation.grid(row=0, column=1, **self.padWE)
        self.transformation.grid(row=0, column=2, **self.padWE)
        self.chk_transformation_label2.grid(row=0, column=3, **self.padWE)
        #
        self.chk_defect_label.grid(row=1, column=0, **self.padWE)
        self.chk_defect.grid(row=1, column=1, **self.padWE)
        self.num_defect.grid(row=1, column=2, **self.padWE)
        self.chk_defect_label2.grid(row=1, column=3, **self.padWE)
        #
        self.chk_first_put_label.grid(row=2, column=0, **self.padWE)
        self.chk_first_put.grid(row=2, column=1, **self.padWE)
        self.put_first.grid(row=2, column=2, **self.padWE)
        self.chk_first_put_label2.grid(row=2, column=3, **self.padWE)
        #
        self.chk_cut_label.grid(row=3, column=0, **self.padWE)
        self.chk_cut.grid(row=3, column=1, **self.padWE)
        self.cut.grid(row=3, column=2, **self.padWE)
        self.chk_cut_label2.grid(row=3, column=3, **self.padWE)
        #
        self.chk_limit_label.grid(row=4, column=0, **self.padWE)
        self.chk_limit.grid(row=4, column=1, **self.padWE)
        self.limit.grid(row=4, column=2, **self.padWE)
        self.chk_limit_label2.grid(row=4, column=3, **self.padWE)

    def create_widgets_records(self) -> None:
        self.record_label = tk.Label(self.frame_records, text="Record")
        self.record = tk.Entry(self.frame_records, width=50)
        self.record.insert(tk.END, self.init_value.record_name)
        self.image_rec_label = tk.Label(self.frame_records, text="Image rec. (%)")
        self.image_rec = tk.Entry(self.frame_records, width=10)
        self.image_rec.insert(tk.END, self.init_value.img_per)
        self.comments_label = tk.Label(self.frame_records, text="Comments")
        self.comments = tk.Entry(self.frame_records, width=80)
        self.comments.insert(tk.END, self.init_value.comments)

    def create_layout_records(self) -> None:
        self.record_label.grid(row=0, column=0, **self.padWE)
        self.record.grid(row=0, column=1, **self.padWE)
        self.image_rec_label.grid(row=0, column=2, **self.padWE)
        self.image_rec.grid(row=0, column=3, **self.padWE)
        self.comments_label.grid(row=1, column=0, **self.padWE)
        self.comments.grid(row=1, column=1, columnspan=3, **self.padWE)

    def create_widgets_method(self) -> None:
        self.var_method = tk.StringVar()
        methods = ("Null event", "Rejection free", "Simple 1D")
        self.method_label = tk.Label(self.frame_method, text="Method")
        self.method_cb = ttk.Combobox(
            self.frame_method,
            textvariable=self.var_method,
            values=methods,
            state="readonly",
        )
        self.method_cb.current(methods.index(self.init_value.method))

    def create_layout_method(self) -> None:
        self.method_label.grid(row=0, column=0, **self.padWE)
        self.method_cb.grid(row=0, column=1, **self.padWE)

    def create_widgets_buttons(self) -> None:
        self.start = tk.Button(
            self.frame_buttons, text="Start", command=self.start_function, width=20
        )
        self.close = tk.Button(
            self.frame_buttons, text="close", command=self.close_function, width=20
        )

    def create_layout_buttons(self) -> None:
        self.start.grid(row=0, column=0)
        self.close.grid(row=0, column=1)

    def create_widgets_progress(self) -> None:
        self.progress_label = ttk.Label(self.frame_progress, text="Progress (%)")
        self.progress_time = ttk.Label(self.frame_progress, text="Sim. time")
        self.progress_coverage = ttk.Label(self.frame_progress, text="Coverage")
        self.progress_atoms = ttk.Label(self.frame_progress, text="Num. atoms")
        self.progress_events = ttk.Label(self.frame_progress, text="Num. events")
        self.progress_expectation = ttk.Label(
            self.frame_progress, text="Expected Num. events"
        )
        self.expectd_cal_time = ttk.Label(
            self.frame_progress, text="Expected cal. time"
        )

    def create_layout_progress(self) -> None:
        self.progress_label.grid(row=0, column=0, padx=4, pady=10)
        self.progress_time.grid(row=0, column=1, padx=4, pady=10)
        self.progress_coverage.grid(row=0, column=2, padx=4, pady=10)
        self.progress_atoms.grid(row=0, column=3, padx=4, pady=10)
        self.progress_events.grid(row=0, column=4, padx=4, pady=10)
        self.progress_expectation.grid(row=0, column=5, padx=4, pady=10)
        self.expectd_cal_time.grid(row=0, column=6, padx=4, pady=10)

    def create_widgets_bar(self) -> None:
        self.pbval = 0
        self.progress_bar = ttk.Progressbar(
            self.frame_bars,
            orient=tk.HORIZONTAL,
            length=350,
            mode="determinate",
        )
        self.progress_bar.configure(maximum=100, value=self.pbval)

    def create_layout_bar(self) -> None:
        self.progress_bar.grid(pady=20)

    def close_function(self) -> None:
        self.quit()

    def run(self) -> None:
        self.mainloop()

    def update_values(self) -> None:
        self.init_value.n_cell_init = int(self.n_cell.get())
        self.init_value.z_unit_init = int(self.z_unit.get())
        self.init_value.temperature = float(self.temperature.get())
        self.init_value.dep_rate = float(self.deposition_rate.get())
        self.init_value.dep_time = float(self.deposition_time.get())

        self.init_value.prefactor = float(self.prefactor.get())
        for energy_entry, energy_key in zip(
            self.energies, self.init_value.binding_energies
        ):
            self.init_value.binding_energies[energy_key] = float(energy_entry.get())
        self.init_value.transformation = float(self.transformation.get())
        self.init_value.record_name = str(self.record.get())
        self.init_value.img_per = float(self.image_rec.get())
        self.init_value.comments = str(self.comments.get())
        self.init_value.keep_defect_check = self.bln_defect.get()
        self.init_value.trans_check = self.bln_tr.get()
        self.init_value.first_put_check = self.bln_first.get()
        self.init_value.cut_check = self.bln_cut.get()
        self.init_value.limit_check = self.bln_limit.get()

        self.init_value.num_defect = self.num_defect.get()
        self.init_value.put_first = self.put_first.get()
        self.init_value.cut_number = self.cut.get()
        self.init_value.limit_val = self.limit.get()

        self.init_value.method = self.var_method.get()
        self.record_middle = 0
        #
        kbt = self.init_value.temperature_eV
        self.temperautre_energy["text"] = "{:.3g}".format(kbt)
        for energy, rate_label in zip(self.energies, self.rate_labels):
            rate_label["text"] = str(
                "{:.3g}".format(
                    rate(
                        float(self.prefactor.get()),
                        kbt,
                        float(energy.get())
                        + float(self.init_value.binding_energies["Base"]),
                    )
                )
            )
        self.dep_rate_conv_val["text"] = "{:.3f}".format(
            self.init_value.dep_rate_atoms_persec
        )

        self.chk_limit_label2["text"] = (
            "(" + "{:.3g}".format(self.init_value.cal_energy) + " eV)"
        )

        self.update()

    def update_click(self, event) -> None:
        self.update_values()

    def start_setting_tk(self):
        self.progress_label["text"] = "Started"
        self.progress_time["text"] = "0 s"
        self.progress_coverage["text"] = "0 ML"
        self.progress_atoms["text"] = "0 atoms"
        self.progress_events["text"] = "0"
        self.pbval = 0
        self.progress_bar.configure(value=self.pbval)
        self.update()

    """
    def start_setting(self) -> None:
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
        #
        #
        self.n_events_rec: List[int] = []
        self.num_atoms_rec: List[int] = []
    """

    def update_progress_tk(self):
        self.pbval = int(self.prog_time / self.init_value.total_time * 100)
        self.progress_bar.configure(value=self.pbval)
        self.progress_label["text"] = str(self.pbval) + " %"
        self.progress_time["text"] = str("{:.2f}".format(self.prog_time)) + " s"
        self.progress_coverage["text"] = (
            str("{:.2f}".format(self.n_atoms / self.init_value.atoms_in_BL)) + " ML"
        )
        self.progress_atoms["text"] = str(self.n_atoms) + " atoms"
        self.progress_events["text"] = str(self.n_events) + " events"
        if (self.n_events == 100) and (self.init_value.method == "Null event"):
            time_middle = time.time() - self.start_time
            expected_time = self.num_events / 100 * time_middle
            self.expectd_cal_time["text"] = (
                str(math.floor(expected_time / 3600))
                + " h"
                + str(int((expected_time % 3600) / 60))
                + " min"
                + str(int(expected_time % 60))
                + " sec"
            )
        elif self.n_events == 100:
            self.expectd_cal_time["text"] = "---"

        self.update()

    """
    def update_progress(self) -> None:
        self.n_events += 1
        self.n_events_perdep += 1
    """

    def det_normarize(self) -> None:
        kbt = self.init_value.temperature_eV
        fast_event = max(
            [
                rate(float(self.prefactor.get()), kbt, float(energy.get()))
                for energy in self.energies
            ]
        )
        if self.init_value.limit_check is True:
            if fast_event > float(self.init_value.limit_val):
                fast_event = float(self.init_value.limit_val)

        if self.bln_tr.get() is True:
            self.normarize = 6 * fast_event + rate(
                float(self.prefactor.get()), kbt, float(self.init_value.transformation)
            )
        else:
            self.normarize = 6 * fast_event

    """
    def record_position(self) -> None:
        self.pos_rec.append(copy.copy(self.atom_set))
        self.time_rec.append(self.prog_time)
        self.cov_rec.append(self.n_atoms / self.init_value.atoms_in_BL)
    """

    def cal_expected_events(self) -> None:
        dep_success = self.init_value.dep_rate_atoms_persec / self.normarize
        num_atom_total = self.init_value.total_atoms + 1
        self.num_events = 1 / dep_success * num_atom_total * num_atom_total / 2
        self.progress_expectation["text"] = (
            "Expected: " + str(int(self.num_events)) + " events"
        )

    """
    def update_after_deposition(self, dep_pos, atom_type) -> None:
        ##
        self.n_events_rec.append(self.n_events_perdep)
        self.num_atoms_rec.append(self.n_atoms)
        ##
        self.n_events_perdep = 0
        self.atom_set[dep_pos] = atom_type
        self.atom_exist.append(dep_pos)
        self.n_atoms += 1
        self.n_events += 1
        self.prog_time += 1 / (self.init_value.dep_rate_atoms_persec)
        if dep_pos[2] in (0, 1):
            self.empty_firstBL -= 1
    """

    def try_deposition(self) -> None:
        # deposition
        judge = judge_null(self.init_value.dep_rate_atoms_persec / self.normarize)
        if judge:
            self.deposition()

    """
    def deposition(self) -> Tuple:
        # self.n_events_perdep = 0
        dep_pos, atom_type = deposit_an_atom(
            self.atom_set,
            self.bonds,
            self.init_value,
        )
        self.update_after_deposition(dep_pos, atom_type)
        return dep_pos
    """
    """
    def defect_check(self):
        if (self.target[2] not in (0, 1)) and (self.move_atom[2] in (0, 1)):
            self.move_atom = self.target
            self.new_state = self.atom_set[self.target]
            self.n_events -= 1
            self.n_events_perdep -= 1
    """

    """
    def event_progress(self):
        if (self.empty_firstBL == int(self.init_value.num_defect)) and (
            self.init_value.keep_defect_check is True
        ):
            self.defect_check()
        #
        if self.move_atom == self.target:
            if self.new_state == 4:
                # クラスターで3次元化
                # print("clustering")
                # print(self.target)
                self.prev_eve = "clustering"

                self.atom_set[self.target] = 3
                for bond in self.bonds[self.target]:
                    if self.atom_set[bond] != 0:
                        self.atom_set[bond] = 3
                        # print(bond)
                # input()
            else:
                self.prev_eve = (
                    "state change: "
                    + str(self.target)
                    + " : "
                    + str(self.atom_set[self.target])
                    + "→"
                    + str(self.new_state)
                )
                self.atom_set[self.target] = self.new_state
    """
    """
            elif self.new_state != self.atom_set[self.target]:
                if self.new_state == 2:
                    print("3D→2D")
                    print("target: " + str(self.target))
                    print("state: " + str(self.atom_set[self.target]))
                    print("new state: " + str(self.new_state))
                    print(self.event_time[self.target])
                    print(self.event[self.target])
                    num_b2 = 0
                    num_b3 = 0
                    for bond in self.bonds[self.target]:
                        if self.atom_set[bond] == 2:
                            num_b2 += 1
                        elif self.atom_set[bond] == 3:
                            num_b3 += 1
                    print("bond2: " + str(num_b2))
                    print("bond3: " + str(num_b3))

                else:
                    print("2D→3D")
                    print("target: " + str(self.target))
                    print("state: " + str(self.atom_set[self.target]))
                    print("new state: " + str(self.new_state))
                    print(self.event_time[self.target])
                    print(self.event[self.target])
                    num_b2 = 0
                    num_b3 = 0
                    for bond in self.bonds[self.target]:
                        if self.atom_set[bond] == 2:
                            num_b2 += 1
                        elif self.atom_set[bond] == 3:
                            num_b3 += 1
                    print("bond2: " + str(num_b2))
                    print("bond3: " + str(num_b3))

                self.atom_set[self.target] = self.new_state
    """
    """
        else:
            self.prev_eve = (
                "move : "
                + str(self.target)
                + " : "
                + str(self.atom_set[self.target])
                + " → "
                + str(self.move_atom)
                + " ("
                + str(self.atom_set[self.move_atom])
                + ") : "
                + str(self.new_state)
            )
            self.atom_set[self.move_atom] = self.new_state
            self.atom_set[self.target] = 0
            self.atom_exist.remove(self.target)
            self.atom_exist.append(self.move_atom)
            if self.move_atom[2] in (0, 1):
                self.empty_firstBL -= 1
            if self.target[2] in (0, 1):
                self.empty_firstBL += 1
    """

    def try_events(self) -> None:
        events, rates, states = site_events(
            self.atom_set,
            self.bonds,
            self.target,
            self.init_value,
            self.bln_defect.get(),
            self.empty_firstBL,
            self.bln_tr.get(),
        )
        # Normarize rates
        norm_rates = normarize_rate(rates, self.normarize)
        # add null event
        events.append(self.target)
        states.append(self.atom_set[self.target])
        # choose an event
        self.move_atom, self.new_state = choice(events, norm_rates, states)
        if (self.target != self.move_atom) or (
            self.atom_set[self.target] != self.new_state
        ):
            self.n_events_perdep += 1

        # event progress
        self.event_progress()

    def end_of_loop_tk(self):
        self.progress_label["text"] = (
            "Finished: " + str(self.minute) + " min " + str(self.second) + " s"
        )
        self.update()

    """
    def end_of_loop(self) -> None:
        self.record_position()
        elapsed_time = time.time() - self.start_time
        self.minute = math.floor(elapsed_time / 60)
        self.second = int(elapsed_time % 60)
        record_data(
            self.pos_rec,
            self.time_rec,
            self.cov_rec,
            self.lattice,
            self.init_value,
            self.minute,
            self.second,
            self.bln_defect.get(),
        )
        #
        #
        rec_events_per_dep(self.n_events_rec, self.num_atoms_rec)
        #
        #
    """

    def null_event_kmc(self) -> None:  # 長過ぎ！　Helper function 作ってコンパクトにしないと見通し悪い。
        self.start_setting_tk()
        self.start_setting()
        self.atom_exist: List[Tuple[int, int, int]] = [(-1, -1, -1)]
        self.lattice, self.bonds, self.atom_set, _, _, _, _ = lattice_form(
            self.init_value
        )
        # return lattice, bonds, atom_set, event, event_time, event_time_tot,event_state
        self.det_normarize()
        self.cal_expected_events()
        """
        最初に2原子置くのは、1原子のみが拡散する時間が無駄なためでしたが、
        Null event methodでは微々たる差なので消しました
        """
        if self.bln_first.get() is True:
            for _ in range(int(self.put_first.get())):
                _ = self.deposition()

        while int(self.prog_time) <= int(self.init_value.total_time):
            self.target = choose_atom(self.atom_exist)
            if (self.n_events_perdep >= int(self.init_value.cut_number)) and (
                self.init_value.cut_check is True
            ):
                _ = self.deposition()

            elif self.target == (-1, -1, -1):
                self.try_deposition()
            else:
                self.try_events()
            # end of an event
            self.update_progress()
            self.update_progress_tk()
            # recoding the positions in the middle
            if self.n_atoms >= self.rec_num_atoms:
                self.rec_num_atoms += self.init_value.rec_num_atom_interval
                self.record_position()

            """
            # 構造の途中確認用
            if self.n_atoms == 80 and self.record_middle == 0:
                from record_for_test import rec_for_test

                self.record_middle = 1
                rec_for_test(self.atom_set, self.bonds, self.lattice)
            """

        # end of the loop
        self.end_of_loop()
        self.end_of_loop_tk()

    """
    def update_events(self):
        self.related_atoms = list(set(self.related_atoms))
        for target_rel in self.related_atoms:
            if self.atom_set[target_rel] == 0:
                events = []
                rates = []
                states = []

            else:
                events, rates, states = site_events(
                    self.atom_set,
                    self.bonds,
                    target_rel,
                    self.init_value,
                    self.bln_defect.get(),
                    self.empty_firstBL,
                    self.bln_tr.get(),
                )

            self.total_event_time -= self.event_time_tot[target_rel]
            self.event[target_rel] = events
            self.event_time[target_rel] = rates
            self.event_time_tot[target_rel] = sum(rates)
            self.event_state[target_rel] = states
            self.total_event_time += self.event_time_tot[target_rel]
        self.related_atoms = []
    """

    """
    def rejection_free_deposition(self):
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
    """

    """
    def rejection_free_event(self):

        if self.target is None or self.event_number is None:
            print("None in rejection free")
            print(self.target)
            print(self.event_number)

        self.move_atom = self.event[self.target][self.event_number]
        self.new_state = self.event_state[self.target][self.event_number]
        self.event_progress()
        self.related_atoms = recalculate(
            self.target, self.bonds, self.atom_set, self.init_value
        )
        self.related_atoms.extend(
            recalculate(self.move_atom, self.bonds, self.atom_set, self.init_value)
        )
        if self.new_state == 4:
            for bond in self.bonds[self.target]:
                self.related_atoms.extend(
                    recalculate(bond, self.bonds, self.atom_set, self.init_value)
                )
        self.update_events()
        #
        #
        # self.prev_eve = str(self.target) + ":" + str(self.move_atom)

    """

    def num_atom_check(self):
        num = 0
        for state in self.atom_set.values():
            if state != 0:
                num += 1
        num_ex = len(self.atom_exist)
        #
        if self.n_atoms != num:
            print("Atom number not match")
            print("Rec. num : " + str(self.n_atoms))
            print("Real num : " + str(num))
            print("Exist num :" + str(num_ex))
            print(self.prev_eve)
            input()
        if num_ex != num:
            print("Atom_number not much with existing list")
            print("Rec. num : " + str(self.n_atoms))
            print("Real num : " + str(num))
            print("Exist num :" + str(num_ex))
            print(self.prev_eve)
            input()

    """
    def atom_count(self):
        if self.empty_firstBL < int(self.init_value.num_defect):
            count = 0
            for key, val in self.atom_set.items():
                if val != 0 and key[2] in (0, 1):
                    count += 1
            print("defect error")
            print(self.prev_eve)
            print(self.empty_firstBL)
            print(int(self.init_value.num_defect))
            print(count)
            print(self.prev_dep)
            input()
        count = 0
        for key, val in self.atom_set.items():
            if val != 0:
                count += 1
        if count != self.n_atoms:
            print("count error")
            print(self.prev_eve)
            input()
    """

    def rejection_free_kmc(self) -> None:
        #
        # self.min_rates = 10000000000
        #
        self.start_setting_tk()
        self.start_setting()
        self.progress_expectation["text"] = "---"
        self.start_rejection_free()
        # 最初の二原子を配置
        self.put_first_atoms_rf()
        #
        # self.prev_eve = "dep"
        while int(self.prog_time) <= int(self.init_value.total_time):
            # self.num_atom_check()
            # self.atom_count()
            self.rejection_free_loop()
            self.update_progress()
            self.update_progress_tk()

            # self.middle_check()
        self.end_of_loop()
        self.end_of_loop_tk()

    def start_function(self) -> None:
        self.update_values()
        if self.var_method.get() == "Null event":
            self.null_event_kmc()
        elif self.var_method.get() == "Rejection free":
            self.rejection_free_kmc()
        elif self.var_method.get() == "Simple 1D":
            pass


if __name__ == "__main__":
    application = tk.Tk()
    application.title("kMC_Si")
    Window(application)
    application.mainloop()
