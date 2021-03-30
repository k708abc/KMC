#!/usr/bin/env python3

from typing import List, Dict, Tuple
import tkinter as tk
import tkinter.ttk as ttk
import math
from Modules.cal_rates import rate
import time
from Modules.kmc_functions import common_functions
from Modules.Null_kmc_functions import null_functions


class Window(ttk.Frame, common_functions, null_functions):
    kb_eV = 8.617e-5
    padWE: Dict = dict(sticky=(tk.W, tk.E), padx=15, pady=2)

    def __init__(self, master) -> None:
        super().__init__(master, padding=2)
        # self.init_value = Params()
        common_functions.__init__(self)
        null_functions.__init__(self)
        self.create_frame_basics()
        self.create_frame_energies()
        self.create_frame_checks()
        self.create_frame_memos()
        self.create_method()
        self.create_frame_buttons()
        self.create_frame_progress()
        self.create_frame_bar()
        self.cerate_frame_repeat()
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

    def cerate_frame_repeat(self):
        self.frame_repeat = ttk.Frame()
        self.create_widgets_repeat()
        self.create_layout_repeat()
        self.frame_repeat.pack()

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
        self.prefactor.insert(tk.END, "{:.1e}".format(self.init_value.prefactor))
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
            "Si_1st BL",
            "Si_2nd BL",
            "Si_3rd BL",
            "Si_else",
        ]

        self.energylabels = []
        self.diffuse = []
        self.energies = []
        self.rate_labels = []
        for _ in range(4):
            self.diffuse.append(ttk.Entry(self.frame_energies, width=7))
            self.energies.append(ttk.Entry(self.frame_energies, width=7))
            self.rate_labels.append(ttk.Label(self.frame_energies, text="0"))
        for (energy, (key, val)) in zip(
            self.energies, self.init_value.binding_energies.items()
        ):
            energy.insert(tk.END, val)
            energy.bind("<Return>", self.update_click)
        for (diffuse, (key, val)) in zip(
            self.diffuse, self.init_value.diffusion_barriers.items()
        ):
            diffuse.insert(tk.END, val)
            diffuse.bind("<Return>", self.update_click)

        for label in self.labels:
            self.energylabels.append(ttk.Label(self.frame_energies, text=label))

    def create_layout_energies(self) -> None:
        # self.update_values()
        for i, energylabel in enumerate(self.energylabels):
            energylabel.grid(row=0, column=i, **self.padWE)
        self.energy_label_diffuse = ttk.Label(
            self.frame_energies, text="Diffusion Energy (eV)"
        )
        self.energy_label_diffuse.grid(row=1, column=0, **self.padWE)
        for i, diffuse_entry in enumerate(self.diffuse):
            diffuse_entry.grid(row=1, column=i + 1)

        self.energy_label0 = ttk.Label(self.frame_energies, text="Bond Energy (eV)")
        self.energy_label0.grid(row=2, column=0, **self.padWE)

        for i, energy_entry in enumerate(self.energies):
            energy_entry.grid(row=2, column=i + 1)
        self.rate_bond_label = ttk.Label(self.frame_energies, text="Diffuse rate")
        self.rate_bond_label.grid(row=3, **self.padWE)
        for i, ratelabel in enumerate(self.rate_labels):
            ratelabel.grid(row=3, column=i + 1)

    def create_widgets_checks(self) -> None:
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
        self.chk_defect_label.grid(row=0, column=0, **self.padWE)
        self.chk_defect.grid(row=0, column=1, **self.padWE)
        self.num_defect.grid(row=0, column=2, **self.padWE)
        self.chk_defect_label2.grid(row=0, column=3, **self.padWE)
        #
        self.chk_first_put_label.grid(row=1, column=0, **self.padWE)
        self.chk_first_put.grid(row=1, column=1, **self.padWE)
        self.put_first.grid(row=1, column=2, **self.padWE)
        self.chk_first_put_label2.grid(row=1, column=3, **self.padWE)
        #
        self.chk_cut_label.grid(row=2, column=0, **self.padWE)
        self.chk_cut.grid(row=2, column=1, **self.padWE)
        self.cut.grid(row=2, column=2, **self.padWE)
        self.chk_cut_label2.grid(row=2, column=3, **self.padWE)
        #
        self.chk_limit_label.grid(row=3, column=0, **self.padWE)
        self.chk_limit.grid(row=3, column=1, **self.padWE)
        self.limit.grid(row=3, column=2, **self.padWE)
        self.chk_limit_label2.grid(row=3, column=3, **self.padWE)

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
        methods = ("Null_event", "Rejection_free")
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
        self.save = tk.Button(
            self.frame_buttons, text="Save values", command=self.save_function, width=20
        )

    def create_layout_buttons(self) -> None:
        self.start.grid(row=0, column=0)
        self.close.grid(row=0, column=1)
        self.save.grid(row=0, column=2)

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

    def create_widgets_repeat(self) -> None:
        self.repeat_label = ttk.Label(self.frame_repeat, text="Repetition setting")
        self.start_label = ttk.Label(self.frame_repeat, text="Start")
        self.end_label = ttk.Label(self.frame_repeat, text="End")
        self.diff_label = ttk.Label(self.frame_repeat, text="Diff")

        self.repeat_label_1 = ttk.Label(self.frame_repeat, text="E1 (eV)")
        self.repeat_label_2 = ttk.Label(self.frame_repeat, text="E2 (eV)")

        self.start_E1 = ttk.Entry(self.frame_repeat, width=7)
        self.start_E1.insert(tk.END, self.init_value.start_E1)

        self.end_E1 = ttk.Entry(self.frame_repeat, width=7)
        self.end_E1.insert(tk.END, self.init_value.end_E1)

        self.diff_E1 = ttk.Entry(self.frame_repeat, width=7)
        self.diff_E1.insert(tk.END, self.init_value.diff_E1)

        self.start_E2 = ttk.Entry(self.frame_repeat, width=7)
        self.start_E2.insert(tk.END, self.init_value.start_E2)

        self.end_E2 = ttk.Entry(self.frame_repeat, width=7)
        self.end_E2.insert(tk.END, self.init_value.end_E2)

        self.diff_E2 = ttk.Entry(self.frame_repeat, width=7)
        self.diff_E2.insert(tk.END, self.init_value.diff_E2)

    def create_layout_repeat(self) -> None:
        self.repeat_label.grid(row=0, column=0, **self.padWE)
        self.start_label.grid(row=0, column=1, **self.padWE)
        self.end_label.grid(row=0, column=2, **self.padWE)
        self.diff_label.grid(row=0, column=3, **self.padWE)

        self.repeat_label_1.grid(row=1, column=0, **self.padWE)
        self.start_E1.grid(row=1, column=1, **self.padWE)
        self.end_E1.grid(row=1, column=2, **self.padWE)
        self.diff_E1.grid(row=1, column=3, **self.padWE)
        self.repeat_label_2.grid(row=2, column=0, **self.padWE)
        self.start_E2.grid(row=2, column=1, **self.padWE)
        self.end_E2.grid(row=2, column=2, **self.padWE)
        self.diff_E2.grid(row=2, column=3, **self.padWE)

    def close_function(self) -> None:
        self.quit()

    def save_function(self):
        self.update_values()
        self.rewrite_input()

    def run(self) -> None:
        self.mainloop()

    def update_values(self) -> None:
        self.init_value.n_cell_init = int(self.n_cell.get())
        self.init_value.z_unit_init = int(self.z_unit.get())
        self.init_value.temperature = float(self.temperature.get())
        self.init_value.dep_rate = float(self.deposition_rate.get())
        self.init_value.dep_time = float(self.deposition_time.get())

        self.init_value.prefactor = float(self.prefactor.get())
        #
        for diffuse_entry, diffuse_key in zip(
            self.diffuse, self.init_value.diffusion_barriers
        ):
            self.init_value.diffusion_barriers[diffuse_key] = float(diffuse_entry.get())
        #
        for energy_entry, energy_key in zip(
            self.energies, self.init_value.binding_energies
        ):
            self.init_value.binding_energies[energy_key] = float(energy_entry.get())
        self.init_value.record_name = str(self.record.get())
        self.init_value.img_per = float(self.image_rec.get())
        self.init_value.comments = str(self.comments.get())
        self.init_value.keep_defect_check = self.bln_defect.get()
        self.init_value.first_put_check = self.bln_first.get()
        self.init_value.cut_check = self.bln_cut.get()
        self.init_value.limit_check = self.bln_limit.get()
        self.init_value.num_defect = self.num_defect.get()
        self.init_value.put_first = self.put_first.get()
        self.init_value.cut_number = self.cut.get()
        self.init_value.limit_val = self.limit.get()
        self.init_value.method = self.var_method.get()
        self.init_value.start_E1 = float(self.start_E1.get())
        self.init_value.end_E1 = float(self.end_E1.get())
        self.init_value.diff_E1 = float(self.diff_E1.get())
        self.init_value.start_E2 = float(self.start_E2.get())
        self.init_value.end_E2 = float(self.end_E2.get())
        self.init_value.diff_E2 = float(self.diff_E2.get())

        self.record_middle = 0
        #
        kbt = self.init_value.temperature_eV
        self.temperautre_energy["text"] = "{:.3g}".format(kbt)
        for diffuse_energy, rate_label in zip(self.diffuse, self.rate_labels):
            E_value = float(diffuse_energy.get())

            rate_label["text"] = str(
                "{:.3g}".format(rate(float(self.prefactor.get()), kbt, E_value))
            )
        self.dep_rate_conv_val["text"] = "{:.3f}".format(
            self.init_value.dep_rate_atoms_persec
        )

        self.chk_limit_label2["text"] = (
            "(" + "{:.3g}".format(self.init_value.cal_energy) + " eV)"
        )

        self.update()

    def rewrite_input(self):
        file_name = "InputParameter.py"
        temp_rec = []
        check_val = 0
        with open(file_name, "r", encoding="utf-8_sig") as f:
            for row in f:
                row_sp = row.split()
                if check_val == 0:
                    if "@property" in row_sp:
                        check_val = 1
                        temp_rec.append(row)
                    elif "self.n_cell_init" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.n_cell_init))
                        )
                    elif "self.z_unit_init" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.z_unit_init))
                        )
                    elif "self.temperature" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.temperature))
                        )
                    elif "self.dep_rate" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.dep_rate))
                        )
                    elif "self.dep_time" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.dep_time))
                        )
                    elif "self.prefactor" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.prefactor))
                        )
                    elif 'self.binding_energies["Si_first"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.binding_energies["Si_first"]),
                            )
                        )
                    elif 'self.binding_energies["Si_second"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.binding_energies["Si_second"]),
                            )
                        )
                    elif 'self.binding_energies["Si_third"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.binding_energies["Si_third"]),
                            )
                        )
                    elif 'self.binding_energies["Si_else"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.binding_energies["Si_else"]),
                            )
                        )
                    elif 'self.diffusion_barriers["Si_first"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.diffusion_barriers["Si_first"]),
                            )
                        )
                    elif 'self.diffusion_barriers["Si_second"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.diffusion_barriers["Si_second"]),
                            )
                        )
                    elif 'self.diffusion_barriers["Si_third"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.diffusion_barriers["Si_third"]),
                            )
                        )
                    elif 'self.diffusion_barriers["Si_else"]' in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2],
                                str(self.init_value.diffusion_barriers["Si_else"]),
                            )
                        )
                    elif "self.put_first" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.put_first))
                        )
                    elif "self.cut_number" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.cut_number))
                        )
                    elif "self.num_defect" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.num_defect))
                        )
                    elif "self.record_name" in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2], '"' + str(self.init_value.record_name) + '"'
                            )
                        )
                    elif "self.img_per" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.img_per))
                        )
                    elif "self.comments" in row_sp:
                        comment = ""
                        for i in range(2, len(row_sp)):
                            comment += row_sp[i]
                            if i != len(row_sp) - 1:
                                comment += " "

                        temp_rec.append(
                            row.replace(
                                comment, '"' + str(self.init_value.comments) + '"'
                            )
                        )

                    elif "self.keep_defect_check" in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2], str(self.init_value.keep_defect_check)
                            )
                        )
                    elif "self.first_put_check" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.first_put_check))
                        )
                    elif "self.cut_check" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.cut_check))
                        )
                    elif "self.limit_check" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.limit_check))
                        )
                    elif "self.limit_val" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.limit_val))
                        )
                    elif "self.method" in row_sp:
                        temp_rec.append(
                            row.replace(
                                row_sp[2], '"' + str(self.init_value.method) + '"'
                            )
                        )
                    elif "self.start_E1" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.start_E1))
                        )
                    elif "self.end_E1" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.end_E1))
                        )
                    elif "self.diff_E1" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.diff_E1))
                        )

                    elif "self.start_E2" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.start_E2))
                        )
                    elif "self.end_E2" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.end_E2))
                        )
                    elif "self.diff_E2" in row_sp:
                        temp_rec.append(
                            row.replace(row_sp[2], str(self.init_value.diff_E2))
                        )

                    else:
                        temp_rec.append(row)
                elif check_val == 1:
                    temp_rec.append(row)

        file_name = "InputParameter.py"
        with open(file_name, "w", encoding="utf-8_sig") as f:
            for rep in temp_rec:
                f.write(rep)

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
            expected_time = int(self.expected_num_events) / 100 * time_middle
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

    def cal_expected_events(self) -> None:
        dep_success = self.init_value.dep_rate_atoms_persec / self.normarize
        if dep_success > 1:
            dep_success = 1
        num_atom_total = self.init_value.total_atoms + 1
        self.expected_num_events = 1 / dep_success * num_atom_total * num_atom_total / 2
        self.progress_expectation["text"] = (
            "Expected: " + str(int(self.expected_num_events)) + " events"
        )

    def end_of_loop_tk(self):
        self.progress_label["text"] = (
            "Finished: " + str(self.minute) + " min " + str(self.second) + " s"
        )
        self.update()

    def null_event_kmc(self) -> None:  # 長過ぎ！　Helper function 作ってコンパクトにしないと見通し悪い。
        self.start_setting_tk()
        self.start_setting()
        self.start_null()
        self.cal_expected_events()
        self.put_first_atoms_null()

        while int(self.prog_time) <= int(self.init_value.total_time):
            self.null_event_loop()
            # end of an event
            self.update_progress()
            self.update_progress_tk()
        # end of the loop
        self.end_of_loop()
        self.end_of_loop_tk()

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
        if self.var_method.get() == "Null_event":
            self.null_event_kmc()
        elif self.var_method.get() == "Rejection_free":
            self.rejection_free_kmc()


if __name__ == "__main__":
    application = tk.Tk()
    application.title("kMC_Si")
    Window(application)
    application.mainloop()
