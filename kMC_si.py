#!/usr/bin/env python3

import tkinter as tk
from preference_window import Window


class App(Window):
    def __init__(self, master):
        super().__init__(master)


if __name__ == "__main__":
    application = tk.Tk()
    app = App(application)
    app.run()
