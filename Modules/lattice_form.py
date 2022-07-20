#! /usr/bin/env python3
from typing import List, Dict, Tuple
from Modules.InputParameter import Params
from Modules.find_candidates import find_candidates


lattice_first: Dict[Tuple, List] = {}
lattice: Dict[Tuple, List] = {}
atom_set: Dict[Tuple, int] = {}
bonds: Dict[Tuple, List] = {}
event: Dict[Tuple, List] = {}
event_time: Dict[Tuple, List] = {}
event_time_tot: List[float] = []
site_list_correspondance: Dict[int, Tuple] = {}
list_site_correspondance: Dict[int, Tuple] = {}
diffuse_candidates: Dict[Tuple, List] = {}
highest_atom: Dict[Tuple, int] = {}


def reset_dicts() -> None:
    lattice_first.clear()
    lattice.clear()
    atom_set.clear()
    bonds.clear()
    event.clear()
    event_time.clear()
    event_time_tot.clear()
    site_list_correspondance.clear()
    list_site_correspondance.clear()
    diffuse_candidates.clear()
    highest_atom.clear()


def form_first_3BL(z_intra: float, z_inter: float):
    for site in lattice_first:
        x, y = site
        lattice_first[site] = [
            [x, y, 0],
            [x + 1 / 3.0, y + 1 / 3.0, z_intra],
            [x + 1 / 3.0, y + 1 / 3.0, z_intra + z_inter],
            [x + 2 / 3.0, y + 2 / 3.0, 2 * z_intra + z_inter],
            [x + 2 / 3.0, y + 2 / 3.0, 2 * (z_inter + z_intra)],
            [x, y, 2 * (z_inter + z_intra) + z_intra],  # for accuracy
        ]


def lattice_full_layers(unit_height: int):
    for site in lattice:
        site_xy: Tuple[int, int] = (site[0], site[1])
        lattice[site] = [
            lattice_first[site_xy][site[2] % 6][0],
            lattice_first[site_xy][site[2] % 6][1],
            lattice_first[site_xy][site[2] % 6][2] + unit_height * (site[2] // 6),
        ]


def neighbor_points(
    site_index: Tuple[int, int, int], z_judge: int, unit_length: int, z_max: int
) -> List[Tuple[int, int, int]]:
    neighbors: List[Tuple[int, int, int]]
    x, y, z = site_index
    if z == z_max:
        neighbors = [
            ((x - 1) % unit_length, y, z - 1),
            (x, (y - 1) % unit_length, z - 1),
            ((x - 1) % unit_length, (y - 1) % unit_length, z - 1),
        ]
    elif z_judge in (0, 2):
        neighbors = [
            ((x - 1) % unit_length, y, z + 1),
            (x, (y - 1) % unit_length, z + 1),
            (x, y, z + 1),
        ]
        if z != 0:
            neighbors.append((x, y, z - 1))
    elif z_judge in (1, 3):
        neighbors = [
            ((x + 1) % unit_length, y, z - 1),
            (x, (y + 1) % unit_length, z - 1),
            (x, y, z - 1),
            (x, y, z + 1),
        ]
    elif z_judge == 4:
        neighbors = [
            ((x + 1) % unit_length, y, z + 1),
            ((x + 1) % unit_length, (y + 1) % unit_length, z + 1),
            (x, (y + 1) % unit_length, z + 1),
            (x, y, z - 1),
        ]
    elif z_judge == 5:
        neighbors = [
            ((x - 1) % unit_length, y, z - 1),
            (x, (y - 1) % unit_length, z - 1),
            ((x - 1) % unit_length, (y - 1) % unit_length, z - 1),
            (x, y, z + 1),
        ]
    else:
        raise RuntimeError("Something wrong happens. check z_judge value")
    return neighbors


def search_bond(unit_length: int, z_max: int):
    # Search for bonding atoms for all the atoms
    for bond_site in bonds:
        z_judge = bond_site[2] % 6
        bonds[bond_site] = neighbor_points(bond_site, z_judge, unit_length, z_max)


def lattice_form(input_params: Params):
    unit_length: int = input_params.cell_size_xy
    z_units: int = input_params.cell_size_z
    z_intra: float = float(input_params.distance_intra)
    z_inter: float = float(input_params.distance_inter)
    unit_height = 3 * (z_intra + z_inter)
    reset_dicts()
    z_max = z_units * 6 - 1
    #
    for i in range(unit_length):
        for j in range(unit_length):
            lattice_first[(i, j)] = []
            highest_atom[(i, j)] = 0
            for k in range(z_units * 6):
                lattice[(i, j, k)] = []
                atom_set[(i, j, k)] = 0
                bonds[(i, j, k)] = []
                event[(i, j, k)] = []
                event_time[(i, j, k)] = []
                event_time_tot.append(0)
                list_site_correspondance[len(event_time_tot)] = (i, j, k)
                site_list_correspondance[(i, j, k)] = len(event_time_tot)
                diffuse_candidates[(i, j, k)] = []
    #
    form_first_3BL(z_intra, z_inter)
    lattice_full_layers(unit_height)
    search_bond(unit_length, z_max)
    for index in diffuse_candidates:
        diffuse_candidates[index] = find_candidates(bonds, index, unit_length, z_max)

    return (
        lattice,
        bonds,
        atom_set,
        event,
        event_time,
        event_time_tot,
        site_list_correspondance,
        list_site_correspondance,
        diffuse_candidates,
        highest_atom,
    )
