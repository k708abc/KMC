from Modules.kmc_functions import common_functions
from Modules.Null_kmc_functions import null_functions


class null_event(common_functions, null_functions):
    def __init__(self) -> None:
        common_functions.__init__(self)
        null_functions.__init__(self)

    def start(self):
        print("Calculation start")
        self.start_setting()
        self.start_null()
        self.put_first_atoms_null()
        print("Loop start")
        while int(self.prog_time) <= int(self.init_value.total_time):
            self.null_event_loop()
            # end of an event
            self.update_progress()
        # end of the loop
        print("Recording")
        self.end_of_loop()
        print("Finished: " + str(self.minute) + " min " + str(self.second) + " sec")
        print(
            "Time/event: "
            + str(round(self.elapsed_time / self.n_events * 1000, 3))
            + " ms"
        )
        input()


if __name__ == "__main__":
    null_class = null_event()
    null_class.start()
