#!/usr/bin/env python3

import tkinter as tk
from preference_window import Window
from InputParameter import Params
from Rejection_free_kmc import rejection_free


class App(Window):
    def __init__(self, master):
        super().__init__(master)


"""
#use tk
if __name__ == "__main__":
    application = tk.Tk()
    app = App(application)
    app.run()
"""

# Direct start

if __name__ == "__main__":
    parameters = Params()
    if parameters.method == "Rejection free":
        rf_class = rejection_free
        rf_class.start()

    elif parameters.method == "Null method":
        print("n")
