#!/usr/bin/env python3

import math
import tkinter as tk
from preference_window import Window
from typing import List


class App(Window):
    def __init__(self, master):
        super().__init__(master)


if __name__ == "__main__":
    unit_x: List[float] = [1, 0, 0]
    unit_y: List[float] = [0.5, math.sqrt(3) / 2, 0]
    unit_z: List[float] = [0, 0, 1]
    application = tk.Tk()
    app = App(application)
    app.run()
