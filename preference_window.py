#!/usr/bin/env python3

from typing import List, Dict
import tkinter as tk
import tkinter.ttk as ttk
from collections import OrderedDict
from cal_rates import rate
from calculation import cal_start

init_values: Dict = dict(
    n_cell_init=5,
    z_unit_init=5,
    temperature=550,
    dep_rate=0.4,
    dep_time=5,
    post_anneal=0,
    prefactor="1E+13",
    AgSi=-1.4,
    Si12=-1.3,
    Si23=-1.4,
    Si34=-1.4,
    Si45=-1.2,
    Si56=-1.4,
    Si_intra=-1.4,
    Si_inter=-1.4,
    Agtop=-1.4,
    transformation=-0.3,
    record_name="KMC_Si_rec",
    img_per=10,
    comments="No comments",
)


class Window(ttk.Frame):
    kb_eV = 8.617e-5
    padWE: Dict = dict(sticky=(tk.W, tk.E), padx=10, pady=2)

    def __init__(self, master):
        super().__init__(master, padding=2)
        self.create_variables()
        self.create_frame_basics()
        self.create_frame_energies()
        self.create_frame_checks()
        self.create_frame_memos()
        self.create_frame_buttons()
        self.create_frame_bar()
        master.title("KMC_Si")

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

    def create_frame_buttons(self):
        self.frame_buttons = ttk.Frame()
        self.create_widgets_buttons()
        self.create_layout_buttons()
        self.frame_buttons.pack()

    def create_frame_bar(self):
        self.frame_bars = ttk.Frame()
        self.create_widgets_progresses()
        self.create_layout_progresses()
        self.frame_bars.pack()

    def create_widgets_basics(self):
        # The first row
        self.n_cell_label = ttk.Label(self.frame_basics, text="Number of cell")
        self.n_cell = ttk.Entry(self.frame_basics, width=7)
        self.n_cell.insert(tk.END, init_values["n_cell_init"])
        self.n_cell.bind("<Return>", self.update_click)
        self.z_unit_label = ttk.Label(self.frame_basics, text="Z unit")
        self.z_unit = ttk.Entry(self.frame_basics, width=7)
        self.z_unit.insert(tk.END, init_values["z_unit_init"])
        self.z_unit.bind("<Return>", self.update_click)
        self.temperature_label = ttk.Label(self.frame_basics, text="T (K)")
        self.temperature = ttk.Entry(self.frame_basics, width=7)
        self.temperature.insert(tk.END, init_values["temperature"])
        self.temperature.bind("<Return>", self.update_click)
        self.temperautre_energy_label = ttk.Label(self.frame_basics, text="kbT")
        temperature = self.temperature.get()
        self.temperautre_energy = ttk.Label(
            self.frame_basics, text="{:.3g}".format(float(temperature) * self.kb_eV)
        )
        # The second row
        self.deposition_rate_label = ttk.Label(
            self.frame_basics, text="Dep. rate(ML/min)"
        )
        self.deposition_rate = ttk.Entry(self.frame_basics, width=7)
        self.deposition_rate.insert(tk.END, init_values["dep_rate"])
        self.deposition_rate.bind("<Return>", self.update_click)
        self.deposition_time_label = ttk.Label(self.frame_basics, text="Dep. time(min")
        self.deposition_time = ttk.Entry(self.frame_basics, width=7)
        self.deposition_time.insert(tk.END, init_values["dep_time"])
        self.deposition_time.bind("<Return>", self.update_click)
        self.postanneal_time_label = ttk.Label(
            self.frame_basics, text="Post aneeal time (min)"
        )
        self.postanneal_time = ttk.Entry(self.frame_basics, width=7)
        self.postanneal_time.insert(tk.END, init_values["post_anneal"])
        self.postanneal_time.bind("<Return>", self.update_click)
        self.prefactor_label = ttk.Label(self.frame_basics, text="Prefactor(1/s)")
        self.prefactor = ttk.Entry(self.frame_basics, width=7)
        self.prefactor.insert(tk.END, init_values["prefactor"])
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
        labels: List[str] = [
            "  ",
            "Ag-Si",
            "Si(1-2)",
            "Si(2-3)",
            "Si(3-4)",
            "Si(4-5)",
            "Si(5-6)",
            "Si(intra)",
            "Si(inter)",
            "Ag(top)",
        ]
        e_values: List[float] = [
            init_values["AgSi"],
            init_values["Si12"],
            init_values["Si23"],
            init_values["Si34"],
            init_values["Si45"],
            init_values["Si56"],
            init_values["Si_intra"],
            init_values["Si_inter"],
            init_values["Agtop"],
        ]
        self.energylabels = []
        self.energies = []
        self.rate_labels = []

        for _ in range(9):
            self.energies.append(
                ttk.Entry(self.frame_energies, width=7)
            )  # entryを9つ生成。_が変数だが、使われていない。
            self.rate_labels.append(ttk.Label(self.frame_energies, text="0"))
        for i, energy in enumerate(self.energies):  # 上で作成したentryそれぞれに値を格納
            energy.insert(tk.END, e_values[i])
            energy.bind("<Return>", self.update_click)
        for label in labels:
            self.energylabels.append(ttk.Label(self.frame_energies, text=label))

    def create_layout_energies(self):
        self.update_values()
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
        self.transformation.insert(tk.END, init_values["transformation"])
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
        self.record.insert(tk.END, init_values["record_name"])
        self.image_rec_label = tk.Label(self.frame_records, text="Image rec. (%)")
        self.image_rec = tk.Entry(self.frame_records, width=10)
        self.image_rec.insert(tk.END, init_values["img_per"])
        self.comments_label = tk.Label(self.frame_records, text="Comments")
        self.comments = tk.Entry(self.frame_records, width=80)
        self.comments.insert(tk.END, init_values["comments"])

    def create_layout_records(self):
        self.record_label.grid(row=0, column=0, **self.padWE)
        self.record.grid(row=0, column=1, **self.padWE)
        self.image_rec_label.grid(row=0, column=2, **self.padWE)
        self.image_rec.grid(row=0, column=3, **self.padWE)
        self.comments_label.grid(row=1, column=0, **self.padWE)
        self.comments.grid(row=1, column=1, columnspan=3, **self.padWE)

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

    def create_widgets_progresses(self):
        self.pbval = 0
        self.progress_label = ttk.Label(self.frame_bars, text="Progress...")
        self.progress_bar = ttk.Progressbar(
            self.frame_bars,
            orient=tk.HORIZONTAL,
            length=350,
            mode="determinate",
        )
        self.progress_bar.configure(maximum=100, value=self.pbval)

    def create_layout_progresses(self):
        self.progress_label.grid(row=0, pady=15)
        self.progress_bar.grid(pady=20)

    def close_function(self):
        # plt.close("all")  <<< FIXME
        self.quit()

    def run(self):
        self.mainloop()

    def update_values(self):
        temperature: str = self.temperature.get()
        kbt = float(temperature) * self.kb_eV
        for energy, rate_label in zip(self.energies, self.rate_labels):

            rate_label["text"] = str(
                "{:.3g}".format(
                    rate(float(self.prefactor.get()), kbt, float(energy.get()))
                )
            )
        self.update()

    def update_click(self, event):
        self.update_values()

    def start_function(self):
        self.update_values()
        cal_start()


if __name__ == "__main__":
    application = tk.Tk()
    application.title("kMC_Si")
    Window(application)
    application.mainloop()
