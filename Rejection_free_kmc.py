from Modules.kmc_functions import common_functions


class rejection_free(common_functions):
    def __init__(self, setting_value) -> None:
        common_functions.__init__(self)
        self.setting_value = setting_value
        """
        setting_value: executed by
            0   Rejection_free_kmc.py
            1   Rejection_free_repetetion.py
            2   RFKMC_multi.py
        """

    def start(self):
        if self.setting_value in (1, 0):
            print("Calculation start")
        self.start_setting()
        self.start_rejection_free()
        # Put first atom
        self.put_first_atoms_rf()
        if self.setting_value in (1, 0):
            print("Loop start")

        while int(self.prog_time) <= int(self.init_value.total_time):
            self.rejection_free_loop()
            self.update_progress()
            #
            # self.middle_check(30)
            # self.isolation_check()
            # self.num_atom_check()

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
            # input()


if __name__ == "__main__":
    rf_class = rejection_free(0)
    rf_class.start()
