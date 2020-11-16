from typing import List, Dict
import tkinter as tk
from tkinter import ttk

def update_values():
    input_params["Unit_length"] = float(entry_unit_length.get())
    input_params["Number_of_z_units"] = float(entry_zunit.get())
    input_params["Temperature"] = float(entry_temperature.get())
    input_params["kbt"] = float(entry_temperature.get())*8.617e-5
    input_params["deposition_rate"] = float(entry_dep_rate.get())
    input_params["deposition_time"] = float(entry_dep_time.get())
    input_params["post_anneal"] = float(entry_post_anneal.get())
    input_params["prefactor"] = float(entry_prefactor.get())
    # calculate deposition speed in atoms/sec
    unit_length = int(entry_unit_length.get())
    deposition_rate = float(entry_dep_rate.get())  # ML/min
    deposition_rate = deposition_rate/60  # ML/s
    atoms_in_BL = 2 *unit_length**2  # atoms/ML
    deposition_rate = deposition_rate*atoms_in_BL  # atoms/s
    input_params["Atoms_in_BL"] = float(atoms_in_BL)
    input_params["dep_rate_(atom/s)"] = float(deposition_rate)
    text_dep_atom_persec["text"] = str("{:.5f}".format(deposition_rate))
    #update rates
    for i in range(len(rates_list)):
        rates_list[i]["text"] = str("{:.5f}".format(cal_rate(float(bonding_entry_list[i].get()))))
    #
    root.update()
def update(event):
    update_values()

def GUI_set():

    #For adjusting the position of entries
    arx: List[int] = [0, 20, 140, 200, 300, 380, 500, 560, 660]
    ary: List[int] = [0, 20, 50, 100, 130, 160, 0, 190, 220, 250, 280, 310]
    #
    root = tk.Tk()
    root.title("kMC_Si")
    root.geometry("800x550")
    #
    text_unit_length = tk.Label(root, text="Number of cells")
    text_unit_length.place(x=arx[1], y=ary[1])
    entry_unit_length = tk.Entry(root, text="Number of cells", width=7)
    entry_unit_length.place(x=arx[2], y=ary[1])
    entry_unit_length.delete(0, tk.END)
    entry_unit_length.insert(tk.END, input_params["Unit_length"] )
    entry_unit_length.bind("<Return>", update)
    #
    text_zunit = tk.Label(root, text="Z units")
    text_zunit.place(x=arx[3], y=ary[1])
    entry_zunit = tk.Entry(root, text="Z units", width=7)
    entry_zunit.place(x=arx[4], y=ary[1])
    entry_zunit.delete(0, tk.END)
    entry_zunit.insert(tk.END, input_params["Number_of_z_units"])
    entry_zunit.bind("<Return>", update)
    #
    text_temperature = tk.Label(root, text="T (K)")
    text_temperature.place(x=arx[5], y=ary[1])
    entry_temperature = tk.Entry(root, text="T (K)", width=7)
    entry_temperature.place(x=arx[6], y=ary[1])
    entry_temperature.delete(0, tk.END)
    entry_temperature.insert(tk.END, input_params["Temperature"])
    entry_temperature.bind("<Return>", update)
    #
    text_kbt = tk.Label(root, text="kbT")
    text_kbt.place(x=arx[7], y=ary[1])
    text_kbt = tk.Label(root, text=str("{:.3g}".format(input_params["kbt"])))
    text_kbt.place(x=arx[8], y=ary[1])
    #
    text_dep_rate = tk.Label(root, text="dep_rate (ML/min)")
    text_dep_rate.place(x=arx[1], y=ary[2])
    entry_dep_rate = tk.Entry(root, text="dep_rate", width=7)
    entry_dep_rate.place(x=arx[2], y=ary[2])
    entry_dep_rate.delete(0, tk.END)
    entry_dep_rate.insert(tk.END, input_params["deposition_rate"])
    entry_dep_rate.bind("<Return>", update)
    #
    text_dep_atom_persec = tk.Label(root, text="0")
    text_dep_atom_persec.place(x=arx[2], y=ary[2] + 30)
    #
    text_dep_time = tk.Label(root, text="Dep.time (min)")
    text_dep_time.place(x=arx[3], y=ary[2])
    entry_dep_time = tk.Entry(root, text="Dep.time", width=7)
    entry_dep_time.place(x=arx[4], y=ary[2])
    entry_dep_time.delete(0, tk.END)
    entry_dep_time.insert(tk.END, input_params["deposition_time"])
    entry_dep_time.bind("<Return>", update)
    #
    text_post_anneal = tk.Label(root, text="Post anneal (min)")
    text_post_anneal.place(x=arx[5], y=ary[2])
    entry_post_anneal = tk.Entry(root, text="Post annealing", width=7)
    entry_post_anneal.place(x=arx[6], y=ary[2])
    entry_post_anneal.delete(0, tk.END)
    entry_post_anneal.insert(tk.END, input_params["post_anneal"] )
    entry_post_anneal.bind("<Return>", update)
    #
    text_prefactor = tk.Label(root, text="prefactor (1/s)")
    text_prefactor.place(x=arx[7], y=ary[2])
    entry_prefactor = tk.Entry(root, text="Prefactor (1/s)", width=7)
    entry_prefactor.place(x=arx[8], y=ary[2])
    entry_prefactor.delete(0, tk.END)
    entry_prefactor.insert(tk.END, '{:.1e}'.format(input_params["prefactor"]))
    entry_prefactor.bind("<Return>", update)    
    #set bonding parameters
    text_bonding = tk.Label(root, text="Energy")
    text_bonding.place(x=20, y=ary[4])
    # Ag-Si
    text_AgSi = tk.Label(root, text="Ag-Si")
    text_AgSi.place(x=100, y=ary[3])
    entry_AgSi = tk.Entry(root, text="Ag-Si", width=7)
    entry_AgSi.place(x=100, y=ary[4])
    entry_AgSi.delete(0, tk.END)
    entry_AgSi.insert(tk.END, input_params["bond_Ag_Si"] )
    entry_AgSi.bind("<Return>", update)
    # Si1-2
    text_Si12 = tk.Label(root, text="Si(1-2)")
    text_Si12.place(x=180, y=ary[3])
    entry_Si12 = tk.Entry(root, text="Si(1-2)", width=7)
    entry_Si12.place(x=180, y=ary[4])
    entry_Si12.delete(0, tk.END)
    entry_Si12.insert(tk.END, input_params["bond_Si1_Si2"])
    entry_Si12.bind("<Return>", update)
    # Si2-3
    text_Si23 = tk.Label(root, text="Si(2-3)")
    text_Si23.place(x=260, y=ary[3])
    entry_Si23 = tk.Entry(root, text="Si(2-3)", width=7)
    entry_Si23.place(x=260, y=ary[4])
    entry_Si23.delete(0, tk.END)
    entry_Si23.insert(tk.END, input_params["bond_Si2_Si3"])
    entry_Si23.bind("<Return>", update)
    # Si3-4
    text_Si34 = tk.Label(root, text="Si(3-4)")
    text_Si34.place(x=340, y=ary[3])
    entry_Si34 = tk.Entry(root, text="Si(3-4)", width=7)
    entry_Si34.place(x=340, y=ary[4])
    entry_Si34.delete(0, tk.END)
    entry_Si34.insert(tk.END, input_params["bond_Si3_Si4"])
    entry_Si34.bind("<Return>", update)
    # Si4-5
    text_Si45 = tk.Label(root, text="Si(4-5)")
    text_Si45.place(x=420, y=ary[3])
    entry_Si45 = tk.Entry(root, text="Si(4-5)", width=7)
    entry_Si45.place(x=420, y=ary[4])
    entry_Si45.delete(0, tk.END)
    entry_Si45.insert(tk.END, input_params["bond_Si4_Si5"])
    entry_Si45.bind("<Return>", update)
    # Si5-6
    text_Si56 = tk.Label(root, text="Si(5-6)")
    text_Si56.place(x=500, y=ary[3])
    entry_Si56 = tk.Entry(root, text="Si(5-6)", width=7)
    entry_Si56.place(x=500, y=ary[4])
    entry_Si56.delete(0, tk.END)
    entry_Si56.insert(tk.END, input_params["bond_Si5_Si6"])
    entry_Si56.bind("<Return>", update)
    # else intra layers
    text_Si_intra = tk.Label(root, text="Si(intra)")
    text_Si_intra.place(x=580, y=ary[3])
    entry_Si_intra = tk.Entry(root, text="Si(intra)", width=7)
    entry_Si_intra.place(x=580, y=ary[4])
    entry_Si_intra.delete(0, tk.END)
    entry_Si_intra.insert(tk.END, input_params["bond_Si_intra"])
    entry_Si_intra.bind("<Return>", update)
    # else inter layers
    text_Si_inter = tk.Label(root, text="Si(inter)")
    text_Si_inter.place(x=660, y=ary[3])
    entry_Si_inter = tk.Entry(root, text="Si(inter)", width=7)
    entry_Si_inter.place(x=660, y=ary[4])
    entry_Si_inter.delete(0, tk.END)
    entry_Si_inter.insert(tk.END, input_params["bond_Si_inter"])
    entry_Si_inter.bind("<Return>", update)
    # Agtop
    text_Ag_top = tk.Label(root, text="Ag(top)")
    text_Ag_top.place(x=740, y=ary[3])
    entry_Ag_top = tk.Entry(root, text="Ag(top)", width=7)
    entry_Ag_top.place(x=740, y=ary[4])
    entry_Ag_top.delete(0, tk.END)
    entry_Ag_top.insert(tk.END, input_params["bond_Ag_top"])
    entry_Ag_top.bind("<Return>", update)
    #bond entry lists
    bonding_entry_list = [
        entry_AgSi, entry_Si12, entry_Si23, entry_Si34, entry_Si45,
        entry_Si56, entry_Si_intra, entry_Si_inter
        ]
    # Reference rates
    text_rates = tk.Label(root, text="Rates/bond")
    text_rates.place(x=20, y=ary[5])
    # Ag-Si
    text_AgSi_rates = tk.Label(root, text="0")
    text_AgSi_rates.place(x=100, y=ary[5])
    # Si1-2
    text_Si12_rates = tk.Label(root, text="0")
    text_Si12_rates.place(x=180, y=ary[5])
    # Si2-3
    text_Si23_rates = tk.Label(root, text="0")
    text_Si23_rates.place(x=260, y=ary[5])
    # Si3-4
    text_Si34_rates = tk.Label(root, text="0")
    text_Si34_rates.place(x=340, y=ary[5])
    # Si4-5
    text_Si45_rates = tk.Label(root, text="0")
    text_Si45_rates.place(x=420, y=ary[5])
    # Si5-6
    text_Si56_rates = tk.Label(root, text="0")
    text_Si56_rates.place(x=500, y=ary[5])
    # else between layers
    text_Si_intra_rates = tk.Label(root, text="0")
    text_Si_intra_rates.place(x=580, y=ary[5])
    # else inter layers
    text_Si_inter_rates = tk.Label(root, text="0")
    text_Si_inter_rates.place(x=660, y=ary[5])
    # Agtop
    text_Agtp_rates = tk.Label(root, text="0")
    text_Agtp_rates.place(x=740, y=ary[5])
    #rates list
    rates_list = [
        text_AgSi_rates, text_Si12_rates, text_Si23_rates, text_Si34_rates,
        text_Si45_rates, text_Si56_rates, text_Si_intra_rates, text_Si_inter_rates
    ]
    #set transformation
    bln_transformation = tk.BooleanVar()
    bln_transformation.set(True)
    chk_transformation = tk.Checkbutton(root, variable=bln_transformation, text="Transformation")
    chk_transformation.place(x=20, y=ary[7])
    entry_transformation = tk.Entry(root, text="transformation", width=7)
    entry_transformation.place(x=200, y=ary[7] + 5)
    entry_transformation.delete(0, tk.END)
    entry_transformation.insert(tk.END, input_params["transformation_energy"])
    entry_transformation.bind("<Return>", update)
    #set keep defects
    bln_keep_defect = tk.BooleanVar()
    bln_keep_defect.set(True)
    chk_keep_defect = tk.Checkbutton(
        root, variable=bln_keep_defect, text="Keep a defect in first layer"
    )
    chk_keep_defect.place(x=20, y=ary[8])
    #Set recording settings
    text_record = tk.Label(root, text="Record")
    text_record.place(x=20, y=ary[9])
    entry_record = tk.Entry(root, text="Name", width=50)
    entry_record.place(x=100, y=ary[9])
    entry_record.delete(0, tk.END)
    entry_record.insert(tk.END, input_params["record_name"])
    #Image recording setting
    text_img_every = tk.Label(root, text="Image rec. (%) :  ")
    text_img_every.place(x=450, y=ary[9])
    entry_img_every = tk.Entry(root, text="img", width=10)
    entry_img_every.place(x=570, y=ary[9])
    entry_img_every.delete(0, tk.END)
    entry_img_every.insert(tk.END, input_params["record_image_every"])
    #Set comments
    text_comment = tk.Label(root, text="Comments")
    text_comment.place(x=20, y=ary[10])
    entry_comment = tk.Entry(root, text="Comments", width=110)
    entry_comment.place(x=100, y=ary[10])
    entry_comment.delete(0, tk.END)
    entry_comment.insert(tk.END, input_params["comments"])
    #Set progress bar
    input_params["progress_value"] = 0
    progress_bar = ttk.Progressbar(root, orient=tk.HORIZONTAL, length=350, mode="determinate")
    progress_bar.configure(maximum=100, value=input_params["progress_value"])
    progress_bar.place(x=350, y=ary[11] + 5)
    #Other statuses
    text_count_progress = tk.Label(root, text="Waiting")
    text_count_progress.place(x=720, y=ary[11] + 5)
    #
    text_time_progress = tk.Label(root, text="time (s)")
    text_time_progress.place(x=370, y=ary[11] + 35)
    #
    text_event = tk.Label(root, text="events")
    text_event.place(x=470, y=ary[11] + 35)
    #
    text_number_atoms = tk.Label(root, text="Num. atoms")
    text_number_atoms.place(x=570, y=ary[11] + 35)
    #
    text_coverage = tk.Label(root, text="Coverage")
    text_coverage.place(x=670, y=ary[11] + 35)







