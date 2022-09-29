# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True

from cal_rates cimport rate
from Calc_grid_index cimport grid_num
from vector_operation cimport search, remove_element, remove_duplicate
from libcpp.vector cimport vector


cdef double total_energy(
    vector[int] atom_set, vector[vector[int]] bonds, int target, vector[double] energy_bonding, vector[double] energy_diffuse, vector[int] highest_atom, vector[vector[int]] index_list, int unit_length
):
    cdef double bond_energy = 0.0
    cdef int x ,y, z, bond, b_z
    x = index_list[target][0]
    y = index_list[target][1]
    z = index_list[target][2]
    if z == 0:
        bond_energy += energy_bonding[0]
    # highest_atom >= 1 if sites above the transformation layer is occupied
    if highest_atom[grid_num(x, y, 0, unit_length)] == 0:
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                b_z = index_list[bond][2]
                bond_energy += energy_bonding[max(z, b_z)]
        return energy_diffuse[z] + bond_energy
    elif highest_atom[grid_num(x, y, 0, unit_length)] >= 1:  # transformtion is activated
        for bond in bonds[target]:
            if atom_set[bond] == 1:
                bond_energy += energy_bonding[energy_bonding.size() - 1]  # bulk bonding energy
        return energy_diffuse[energy_diffuse.size() - 1] + bond_energy  # bulk diffusion energy


cdef vector[int] find_filled_sites(vector[int] atom_set, vector[int] indexes):
    cdef int atom_nn
    cdef vector[int] filled_list
    for atom_nn in indexes:
        if atom_set[atom_nn] == 1:
            filled_list.push_back(atom_nn)
    return filled_list


cdef vector[int] find_empty_sites(vector[int] atom_set, vector[int] indexes):
    cdef int atom_nn
    cdef vector[int] empty_list
    for atom_nn in indexes:
        if atom_set[atom_nn] == 0:
            empty_list.push_back(atom_nn)
    return empty_list


cdef (vector[int], bint) check_cluster(int c_target, vector[int] atom_set, vector[vector[int]] bonds, vector[int] cluster_start, vector[vector[int]] index_list):
    cdef vector[int] cluster, next_atoms, check_atoms, empty_list
    cdef int bond, atom, c_z, b_z

    cluster.insert(cluster.end(), cluster_start.begin(), cluster_start.end())
    cluster.push_back(c_target)
    next_atoms.push_back(c_target)
    c_z = index_list[c_target][2]
    if c_z == 0:
        return empty_list, False
    while next_atoms.size() != 0:
        check_atoms.clear()
        check_atoms.insert(check_atoms.end(), next_atoms.begin(), next_atoms.end())
        next_atoms.clear()
        for atom in check_atoms:
            for bond in bonds[atom]:
                if (atom_set[bond] == 1) and (search(cluster, bond) is False):
                    cluster.push_back(bond)
                    next_atoms.push_back(bond)
                    b_z = index_list[bond][2]
                    if (b_z == 0) or (cluster.size() >= 5):
                        return cluster, False
    return cluster, True


# Remove events which cause isolated atoms
# eventsのうち、拡散後に結合がなくなるものを排除する
"""
cdef vector[int] judge_isolation(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    vector[int] nn_atoms,
    vector[int] events,
    vector[vector[int]] index_list
):
    cdef vector[int] checked, remove, cluster, nn_nn_list, common, nonused, empty
    cdef vector[vector[int]] clusters
    cdef int nn_atom, event, nn_nn, i, e_z
    cdef bint z_judge
    #拡散後にはtargetから原子は無くなる
    atom_set[target] = 0
    #nn_atomは、targetの結合サイトで占有されているもの。targetの消失で浮く可能性あり
    for nn_atom in nn_atoms:
        #クラスターに属するか未確認のサイトについて
        if search(checked, nn_atom) is False:
            cluster, z_judge = check_cluster(nn_atom, atom_set, bonds, empty, index_list)
            #checkedはクラスターに属するか確認済みのサイト
            checked.insert(checked.end(), cluster.begin(), cluster.end())
            #浮いてるクラスターがclustersに
            if z_judge is True:
                clusters.push_back(cluster)
        else:
            pass
    #
    if clusters.size() == 0:
        for event in events:
            nn_nn_list.clear()
            e_z = index_list[event][2]
            for nn_nn in bonds[event]:
                if atom_set[nn_nn] == 1:
                    nn_nn_list.push_back(nn_nn)
            if (nn_nn_list.size() == 0) and (e_z != 0):
                remove.push_back(event)
    else:
        for event in events:
            atom_set[event] = 1
            for cluster in clusters:
                common.clear()
                for bond in bonds[event]:
                    if search(cluster, bond):
                        common.push_back(bond)
                if common.size() == 0:
                    remove.push_back(event)
                    break
                else:
                    nonused, z_judge = check_cluster(event, atom_set, bonds, cluster, index_list)
                    if z_judge:
                        remove.push_back(event)
            atom_set[event] = 0
    remove = remove_duplicate(remove)
    events = remove_element(events, remove)
    atom_set[target] = 1
    return events
"""

cdef vector[vector[int]] clusters_form(vector[int] atom_set, vector[int] nn_atoms, vector[vector[int]] bonds, vector[vector[int]] index_list):
    cdef vector[int] checked, cluster, empty
    cdef vector[vector[int]] clusters
    cdef int nn_atom
    cdef bint z_judge
    #拡散後にはtargetから原子は無くなる
    atom_set[target] = 0
    #nn_atomは、targetの結合サイトで占有されているもの。targetの消失で浮く可能性あり
    for nn_atom in nn_atoms:
        #クラスターに属するか未確認のサイトについて
        if search(checked, nn_atom) is False:
            cluster, z_judge = check_cluster(nn_atom, atom_set, bonds, empty, index_list)
            #checkedはクラスターに属するか確認済みのサイト
            checked.insert(checked.end(), cluster.begin(), cluster.end())
            #浮いてるクラスターがclustersに
            if z_judge is True:
                clusters.push_back(cluster)
        else:
            pass
    return clusters


cdef bint judge_isolation2(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    int cand,
    vector[vector[int]] index_list,
    vector[vector[int]] clusters
):
    cdef vector[int]  remove, cluster, nn_nn_list, common, nonused
    cdef int event, nn_nn, i, e_z
    cdef bint z_judge
    #拡散後にはtargetから原子は無くなる
    atom_set[target] = 0
    #
    if clusters.size() == 0:
        e_z = index_list[cand][2]
        for nn_nn in bonds[cand]:
            if atom_set[nn_nn] == 1:
                nn_nn_list.push_back(nn_nn)
        if (nn_nn_list.size() == 0) and (e_z != 0):
            #孤立する
            return False
        else:
            #孤立しない
            return True
    else:
        atom_set[cand] = 1
        for cluster in clusters:
            common.clear()
            for bond in bonds[cand]:
                if search(cluster, bond):
                    common.push_back(cand)
            if common.size() == 0:
                #孤立する
                return False
            else:
                nonused, z_judge = check_cluster(cand, atom_set, bonds, cluster, index_list)
                if z_judge:
                    #孤立する
                    return False
        return True


# targetで起こりうるベントを求める
cdef vector[double] possible_events(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    double pre,
    double kbt,
    double rearange_rate,
    vector[int] diffuse_candidates,
    vector[vector[int]] index_list
):
    cdef vector[double] eve_rate = [0]*diffuse_candidates.size()
    cdef vector[int] nn_atom, event_f, event_test
    cdef vector[double] rates_f
    cdef vector[vector[int]] clusters
    cdef int cand, nonused
    cdef int i, k, j
    nn_atom = find_filled_sites(atom_set, bonds[target])
    #
    if rearrange_rate == 0:
        return eve_rate
    else:
        clusters = clusters_form(atom_set, nn_atoms, bonds, index_list)
        i = 0
        for cand in diffuse_candidates:
            if atom_set[cand] == 1:
                pass
            else:
                if judge_isolation2(atom_set, bonds, target, cand, index_list, clusters):
                    eve_rate[i] = rearange_rate
                else:
                    pass
            i += 1
    return eve_rate


cdef vector[double] site_events(
    vector[int] atom_set,
    vector[vector[int]] bonds,
    int target,
    int unit_length,
    double pre,
    double kbt,
    vector[double] energy_bonding,
    vector[double] energy_diffuse,
    vector[vector[int]] diffuse_candidates,
    vector[int] highest_atom,
    vector[vector[int]] index_list
):
    cdef vector[double] event_list
    # cdef vector[double] rate_list
    cdef double energy, rearange_rate = 0
    cdef int bond, atom, i
    if atom_set[target] == 1:
        energy = total_energy(
            atom_set,
            bonds,
            target,
            energy_bonding,
            energy_diffuse,
            highest_atom,
            index_list,
            unit_length
        )
        rearange_rate = rate(pre, kbt, energy)
    else:
        rearange_rate = 0
    # calculate possible events
    event_list = possible_events(
        atom_set, bonds, target, pre, kbt, rearange_rate, diffuse_candidates[target], index_list
    )
    return event_list
