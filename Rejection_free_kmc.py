from kmc_functions import common_functions


class rejection_free(common_functions):
    def __init__(self) -> None:
        common_functions.__init__(self)

    def start(self):
        print("cal start")
        self.start_setting()
        self.start_rejection_free()
        # 最初の二原子を配置
        self.put_first_atoms_rf()
        #
        # self.prev_eve = "dep"
        print("loop start")
        while int(self.prog_time) <= int(self.init_value.total_time):
            # self.num_atom_check()
            # self.atom_count()
            self.rejection_free_loop()
            self.update_progress()
            # self.middle_check()
        print("recording")
        self.end_of_loop()
        print("Finished: " + str(self.minute) + " min " + str(self.second) + " sec")


if __name__ == "__main__":
    rf_class = rejection_free()
    rf_class.start()
