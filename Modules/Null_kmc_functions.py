from typing import List, Dict, Tuple
from Modules.cal_rates import rate
from Modules.lattice_form import lattice_form
from Modules.choose_site import choose_atom
from Modules.judgement import judge_null
from Modules.event_collection import site_events
from Modules.normarize_list import normarize_rate
from Modules.weighted_choice import choice
import decimal


class null_functions:
    def __init__(self) -> None:
        pass

    def det_normarize(self) -> None:
        kbt = self.init_value.temperature_eV
        pre = float(self.init_value.prefactor)
        # E1 = float(self.init_value.binding_energies["AgSi"])
        E2 = float(self.init_value.binding_energies["Si12"]) + float(
            self.init_value.binding_energies["Base"]
        )
        E3 = float(self.init_value.binding_energies["Si23"]) + float(
            self.init_value.binding_energies["Base"]
        )
        E4 = float(self.init_value.binding_energies["Si34"]) + float(
            self.init_value.binding_energies["Base"]
        )
        E5 = float(self.init_value.binding_energies["Si45"]) + float(
            self.init_value.binding_energies["Base"]
        )
        E6 = float(self.init_value.binding_energies["Si56"]) + float(
            self.init_value.binding_energies["Base"]
        )
        E7 = float(self.init_value.binding_energies["Si_intra"]) + float(
            self.init_value.binding_energies["Base"]
        )
        E8 = float(self.init_value.binding_energies["Si_inter"]) + float(
            self.init_value.binding_energies["Base"]
        )
        fast_event = max(
            [
                rate(pre, kbt, E2),
                rate(pre, kbt, E3),
                rate(pre, kbt, E4),
                rate(pre, kbt, E5),
                rate(pre, kbt, E6),
                rate(pre, kbt, E7),
                rate(pre, kbt, E8),
            ]
        )
        if self.init_value.limit_check is True:
            if fast_event > float(self.init_value.limit_val):
                fast_event = float(self.init_value.limit_val)

        if self.init_value.trans_check is True:
            self.normarize = 6 * fast_event + rate(
                pre, kbt, float(self.init_value.transformation)
            )
        else:
            self.normarize = decimal.Decimal(6 * fast_event)

    def try_deposition(self) -> None:
        # deposition
        judge = judge_null(self.init_value.dep_rate_atoms_persec / self.normarize)
        if judge:
            print(
                "Progress: "
                + str(self.n_atoms)
                + "/"
                + str(self.init_value.total_atoms)
            )
            self.deposition()

    def try_events(self) -> None:
        events, rates, states = site_events(
            self.atom_set,
            self.bonds,
            self.target,
            self.init_value,
        )
        # Normarize rates
        norm_rates = normarize_rate(rates, self.normarize)
        # add null event
        events.append(self.target)
        states.append(self.atom_set[self.target])
        # choose an event
        self.move_atom, self.new_state = choice(events, norm_rates, states)
        if (self.target != self.move_atom) or (
            self.atom_set[self.target] != self.new_state
        ):
            self.n_events_perdep += 1

        # event progress
        self.event_progress()

    def start_null(self):
        self.atom_exist: List[Tuple[int, int, int]] = [(-1, -1, -1)]
        self.lattice, self.bonds, self.atom_set, _, _, _, _ = lattice_form(
            self.init_value
        )
        # return lattice, bonds, atom_set, event, event_time, event_time_tot,event_state
        self.det_normarize()

    def put_first_atoms_null(self):
        if self.init_value.first_put_check is True:
            for _ in range(int(self.init_value.put_first)):
                _ = self.deposition()

    def null_event_loop(self):
        self.target = choose_atom(self.atom_exist)
        if (self.n_events_perdep >= int(self.init_value.cut_number)) and (
            self.init_value.cut_check is True
        ):
            _ = self.deposition()

        elif self.target == (-1, -1, -1):
            self.try_deposition()
        else:
            self.try_events()

        # recoding the positions in the middle
        if self.n_atoms >= self.rec_num_atoms:
            self.rec_num_atoms += self.init_value.rec_num_atom_interval
            self.record_position()

        """
        # 構造の途中確認用
        if self.n_atoms == 80 and self.record_middle == 0:
            from record_for_test import rec_for_test

            self.record_middle = 1
            rec_for_test(self.atom_set, self.bonds, self.lattice)
        """