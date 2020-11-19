#!/usr/bin/env python3

from typing import List, Dict
import tkinter as tk
import tkinter.ttk as ttk
from collections import OrderedDict
from cal_rates import rate
from calculation import cal_start
import time
from lattice_form import lattice_form
from lattice_form_check import check
from InputParameter import Params
from depsoition import deposit_an_atom


class Window(ttk.Frame):
    kb_eV = 8.617e-5
    padWE: Dict = dict(sticky=(tk.W, tk.E), padx=15, pady=2)

    def __init__(self, master):
        super().__init__(master, padding=2)
        self.init_value = Params()
        self.create_variables()
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

    def create_variables(self):
        pass

    def create_frame_basics(self):
        self.frame_basics = ttk.Frame()
        self.create_widgets_basics()
        self.create_layout_basic()
        self.frame_basics.pack()

    def create_frame_energies(self):
        self.frame_energies = ttk.Frame()
        self.create_widgets_energies()
        self.create_layout_energies()
        self.frame_energies.pack()

    def create_frame_checks(self):
        self.frame_checks = ttk.Frame()
        self.create_widgets_checks()
        self.create_layout_checks()
        self.frame_checks.pack()

    def create_frame_memos(self):
        self.frame_records = ttk.Frame()
        self.create_widgets_records()
        self.create_layout_records()
        self.frame_records.pack()

    def create_method(self):
        self.frame_method = ttk.Frame()
        self.create_widgets_method()
        self.create_layout_method()
        self.frame_method.pack()

    def create_frame_buttons(self):
        self.frame_buttons = ttk.Frame()
        self.create_widgets_buttons()
        self.create_layout_buttons()
        self.frame_buttons.pack()

    def create_frame_bar(self):
        self.frame_bars = ttk.Frame()
        self.create_widgets_bar()
        self.create_layout_bar()
        self.frame_bars.pack()

    def create_frame_progress(self):
        self.frame_progress = ttk.Frame()
        self.create_widgets_progress()
        self.create_layout_progress()
        self.frame_progress.pack()

    def create_widgets_basics(self):
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
        self.deposition_time_label = ttk.Label(self.frame_basics, text="Dep. time(min")
        self.deposition_time = ttk.Entry(self.frame_basics, width=7)
        self.deposition_time.insert(tk.END, self.init_value.dep_time)
        self.deposition_time.bind("<Return>", self.update_click)
        self.postanneal_time_label = ttk.Label(
            self.frame_basics, text="Post aneeal time (min)"
        )
        self.postanneal_time = ttk.Entry(self.frame_basics, width=7)
        self.postanneal_time.insert(tk.END, self.init_value.post_anneal)
        self.postanneal_time.bind("<Return>", self.update_click)
        self.prefactor_label = ttk.Label(self.frame_basics, text="Prefactor(1/s)")
        self.prefactor = ttk.Entry(self.frame_basics, width=7)
        self.prefactor.insert(tk.END, self.init_value.prefactor)
        self.prefactor.bind("<Return>", self.update_click)

    def create_layout_basic(self):
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
        self.deposition_time_label.grid(row=1, column=2, **self.padWE)
        self.deposition_time.grid(row=1, column=3, **self.padWE)
        self.postanneal_time_label.grid(row=1, column=4, **self.padWE)
        self.postanneal_time.grid(row=1, column=5, **self.padWE)
        self.prefactor_label.grid(row=1, column=6, **self.padWE)
        self.prefactor.grid(row=1, column=7, **self.padWE)

    def create_widgets_energies(self):
        self.labels: List[str] = [
            "  ",
            "AgSi",
            "Si12",
            "Si23",
            "Si34",
            "Si45",
            "Si56",
            "Si_intra",
            "Si_inter",
            "Ag_top",
        ]

        self.energylabels = []
        self.energies = []
        self.rate_labels = []
        for _ in range(9):
            self.energies.append(
                ttk.Entry(self.frame_energies, width=7)
            )  # entryを9つ生成。_が変数だが、使われていない。
            self.rate_labels.append(ttk.Label(self.frame_energies, text="0"))
        for (energy, (key, val)) in zip(
            self.energies, self.init_value.binding_energies.items()
        ):  # 上で作成したentryそれぞれに値を格納
            energy.insert(tk.END, val)
            energy.bind("<Return>", self.update_click)
        for label in self.labels:
            self.energylabels.append(ttk.Label(self.frame_energies, text=label))

    def create_layout_energies(self):
        # self.update_values()
        for i, energylabel in enumerate(self.energylabels):  # インデックスとリストの要素を同時に取得しループ
            energylabel.grid(row=0, column=i, **self.padWE)
        self.energy_label0 = ttk.Label(self.frame_energies, text="Energy (eV)")
        self.energy_label0.grid(row=1, column=0, **self.padWE)
        for i, energy_entry in enumerate(self.energies):
            energy_entry.grid(row=1, column=i + 1)
        self.rate_bond_label = ttk.Label(self.frame_energies, text="Rates/bond")
        self.rate_bond_label.grid(row=2, **self.padWE)
        for i, ratelabel in enumerate(self.rate_labels):
            ratelabel.grid(row=2, column=i + 1)

    def create_widgets_checks(self):
        self.bln_tr = tk.BooleanVar()
        self.bln_tr.set(True)
        self.chk_transformation_label = ttk.Label(
            self.frame_checks, text="Transformation"
        )
        self.chk_transformation = ttk.Checkbutton(
            self.frame_checks, variable=self.bln_tr
        )
        self.transformation = ttk.Entry(self.frame_checks, width=7)
        self.transformation.insert(tk.END, self.init_value.transformation)
        self.bln_defect = tk.BooleanVar()
        self.bln_defect.set(True)
        self.chk_defect_label = ttk.Label(
            self.frame_checks, text="Keep defect in the first layer"
        )
        self.chk_defect = ttk.Checkbutton(self.frame_checks, variable=self.bln_defect)

    def create_layout_checks(self):
        self.chk_transformation_label.grid(row=0, column=0, **self.padWE)
        self.chk_transformation.grid(row=0, column=1, **self.padWE)
        self.transformation.grid(row=0, column=2, **self.padWE)
        self.chk_defect_label.grid(row=1, column=0, **self.padWE)
        self.chk_defect.grid(row=1, column=1, **self.padWE)

    def create_widgets_records(self):
        self.record_label = tk.Label(self.frame_records, text="Record")
        self.record = tk.Entry(self.frame_records, width=50)
        self.record.insert(tk.END, self.init_value.record_name)
        self.image_rec_label = tk.Label(self.frame_records, text="Image rec. (%)")
        self.image_rec = tk.Entry(self.frame_records, width=10)
        self.image_rec.insert(tk.END, self.init_value.img_per)
        self.comments_label = tk.Label(self.frame_records, text="Comments")
        self.comments = tk.Entry(self.frame_records, width=80)
        self.comments.insert(tk.END, self.init_value.comments)

    def create_layout_records(self):
        self.record_label.grid(row=0, column=0, **self.padWE)
        self.record.grid(row=0, column=1, **self.padWE)
        self.image_rec_label.grid(row=0, column=2, **self.padWE)
        self.image_rec.grid(row=0, column=3, **self.padWE)
        self.comments_label.grid(row=1, column=0, **self.padWE)
        self.comments.grid(row=1, column=1, columnspan=3, **self.padWE)

    def create_widgets_method(self):
        self.var_method = tk.StringVar()
        method_list = ["Check lattice", "Null event", "Rejection free"]
        self.method_label = tk.Label(self.frame_method, text="Method")
        self.method_cb = ttk.Combobox(
            self.frame_method,
            textvariable=self.var_method,
            values=method_list,
            state="readonly",
        )
        self.method_cb.current(0)

    def create_layout_method(self):
        self.method_label.grid(row=0, column=0, **self.padWE)
        self.method_cb.grid(row=0, column=1, **self.padWE)

    def create_widgets_buttons(self):
        self.start = tk.Button(
            self.frame_buttons, text="Start", command=self.start_function, width=20
        )
        self.close = tk.Button(
            self.frame_buttons, text="close", command=self.close_function, width=20
        )

    def create_layout_buttons(self):
        self.start.grid(row=0, column=0)
        self.close.grid(row=0, column=1)

    def create_widgets_progress(self):
        self.progress_label = ttk.Label(self.frame_progress, text="Progress (%)")
        self.progress_time = ttk.Label(self.frame_progress, text="Sim. time")
        self.progress_coverage = ttk.Label(self.frame_progress, text="Coverage")
        self.progress_atoms = ttk.Label(self.frame_progress, text="Num. atoms")
        self.progress_events = ttk.Label(self.frame_progress, text="Num. events")

    def create_layout_progress(self):
        self.progress_label.grid(row=0, column=0, padx=4, pady=10)
        self.progress_time.grid(row=0, column=1, padx=4, pady=10)
        self.progress_coverage.grid(row=0, column=2, padx=4, pady=10)
        self.progress_atoms.grid(row=0, column=3, padx=4, pady=10)
        self.progress_events.grid(row=0, column=4, padx=4, pady=10)

    def create_widgets_bar(self):
        self.pbval = 0
        self.progress_bar = ttk.Progressbar(
            self.frame_bars,
            orient=tk.HORIZONTAL,
            length=350,
            mode="determinate",
        )
        self.progress_bar.configure(maximum=100, value=self.pbval)

    def create_layout_bar(self):
        self.progress_bar.grid(pady=20)

    def close_function(self):
        self.quit()

    def run(self):
        self.mainloop()

    def update_values(self):
        self.init_value.n_cell_init = int(self.n_cell.get())
        self.init_value.z_unit_init = int(self.z_unit.get())
        self.init_value.temperature = float(self.temperature.get())
        self.init_value.dep_rate = float(self.deposition_rate.get())
        self.init_value.dep_time = float(self.deposition_time.get())
        self.init_value.post_anneal = float(self.postanneal_time.get())
        self.init_value.prefactor = float(self.prefactor.get())
        for energy_entry, key in zip(self.energies, self.init_value.binding_energies):
            self.init_value.binding_energies[key] = float(energy_entry.get())
        self.init_value.transformation = float(self.transformation.get())
        self.init_value.record_name = str(self.record.get())
        self.init_value.img_per = float(self.image_rec.get())
        self.init_value.comments = str(self.comments.get())
        #
        kbt = self.init_value.temperature_eV
        self.temperautre_energy["text"] = "{:.3g}".format(kbt)
        for energy, rate_label in zip(self.energies, self.rate_labels):
            rate_label["text"] = str(
                "{:.3g}".format(
                    rate(float(self.prefactor.get()), kbt, float(energy.get()))
                )
            )
        self.update()

    def update_click(self, event):
        self.update_values()

    def lattice_check(self):
        lattice_formed = lattice_form(self.init_value)
        check(self.init_value, lattice_formed[0], lattice_formed[1])

    def start_setting(self):
        self.progress_label["text"] = "Started"
        self.progress_time["text"] = "0 s"
        self.progress_coverage["text"] = "0 ML"
        self.progress_atoms["text"] = "0 atoms"
        self.progress_events["text"] = "0"
        self.pbval = 0
        self.progress_bar.configure(value=self.pbval)
        self.update()
        self.start_time = time.time()
        self.prog_time = 0
        self.n_atoms = 0
        self.n_events = 0

    def update_progress(self):
        self.pbval = int(self.prog_time / self.init_value.total_time * 100)
        self.progress_bar.configure(value=self.pbval)
        self.progress_label["text"] = str(self.pbval) + " %"
        self.progress_time["text"] = str("{:.2f}".format(self.prog_time)) + " s"
        self.progress_coverage["text"] = (
            str("{:.2f}".format(self.n_atoms / self.init_value.atoms_in_BL)) + " ML"
        )
        self.progress_atoms["text"] = str(self.n_atoms) + "atoms"
        self.progress_events["text"] = str(self.n_events) + "events"
        self.update()

    def start_function(self):
        self.update_values()
        if self.var_method.get() == "Check lattice":
            self.lattice_check()

        elif self.var_method.get() == "Null event":
            self.start_setting()
            lattice, bonds, atom_set, _, _, _ = lattice_form(self.init_value)
            # return lattice, bonds, atom_set, event, event_time, event_time_tot
            # put first and second atom
            for _ in 2:
                dep_pos, atom_type = deposit_an_atom(atom_set, bonds)
                atom_set[dep_pos] = atom_type
                self.n_atoms += 1
                self.n_events += 1
            self.update_progress()

        elif self.var_method.get() == "Rejection free":
            self.start_setting()


if __name__ == "__main__":
    application = tk.Tk()
    application.title("kMC_Si")
    Window(application)
    application.mainloop()
