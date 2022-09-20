# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from typing import ValuesView
from Calc_grid_index cimport grid_num

def height_check(pos_x, pos_y, pos, maxz, unit_length):
    max_pos = -1
    for z in range(maxz + 1):
        index = grid_num(pos_x, pos_y, z, unit_length)
        if pos[index] != 0:
            max_pos = z
    return max_pos


def growth_check(pos, unit_length, maxz, atom_BL, index_list):
    num_ag = 0
    num_1st = 0
    num_2nd = 0
    num_multi = 0
    for i in range(unit_length):
        for k in range(unit_length):
            height = height_check(i, k, pos, maxz, unit_length)
            if height == -1:
                num_ag += 2
            elif height == 0:
                num_ag += 1
                num_1st += 1
            elif height == 1:
                num_1st += 2
            elif height == 2:
                num_1st += 1
                num_2nd += 1
            elif height == 3:
                num_2nd += 2
            elif height == 4:
                num_2nd += 1
                num_multi += 1
            else:
                num_multi += 2
    return (
        round(num_ag / atom_BL * 100, 2),
        round(num_1st / atom_BL * 100, 2),
        round(num_2nd / atom_BL * 100, 2),
        round(num_multi / atom_BL * 100, 2),
    )


def num_check(coverage, num):
    value = 100
    number = 0
    for i in range(len(coverage)):
        cf_val = abs(num - coverage[i])
        if value > cf_val:
            number = i
            value = cf_val
    return number


#
# 以前まで用いていた方法
def growth_val_prev(growth_mode, coverage, num_1ML, num_2ML):
    First_at_1ML = growth_mode[num_1ML][1]
    Second_at_2ML = growth_mode[num_2ML][2]
    if growth_mode[num_2ML][3] == 100:
        First_at_second = 0
    else:
        First_at_second = (
            growth_mode[num_2ML][1] / (100 - growth_mode[num_2ML][3]) * 100
        )
    if First_at_1ML < 50:
        if Second_at_2ML < 50:
            mode = "VW"
        else:
            mode = "BL"
    else:
        if Second_at_2ML >= 50:
            mode = "FM"
        else:
            if First_at_second >= 50:
                mode = "SK"
            else:
                mode = "DW"

    return [
        First_at_1ML,
        Second_at_2ML,
        First_at_second,
        mode,
        coverage[num_1ML],
        coverage[num_2ML],
    ]


def mode_check_prev(growth_mode, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    return growth_val_prev(growth_mode, coverage, num_1ML, num_2ML)


#
#
# 分母に2層目+3層目を使用(二層目以上を除いたときの1層目の割合)
def growth_val_1(growth_mode, coverage, num_1ML, num_2ML):
    First_at_1ML = growth_mode[num_1ML][1]
    Second_at_2ML = growth_mode[num_2ML][2]
    if (growth_mode[num_2ML][2] + growth_mode[num_2ML][3]) == 100:
        First_at_second = 0
    else:
        First_at_second = (
            growth_mode[num_2ML][1]
            / (100 - growth_mode[num_2ML][3] - growth_mode[num_2ML][2])
            * 100
        )
    if First_at_1ML < 50:
        if Second_at_2ML < 50:
            mode = "VW"
        else:
            mode = "BL"
    else:
        if Second_at_2ML >= 50:
            mode = "FM"
        else:
            if First_at_second >= 50:
                mode = "SK"
            else:
                mode = "DW"

    return [
        First_at_1ML,
        Second_at_2ML,
        First_at_second,
        mode,
        coverage[num_1ML],
        coverage[num_2ML],
    ]


def mode_check_1(growth_mode, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    return growth_val_1(growth_mode, coverage, num_1ML, num_2ML)


#
#
# 分子にAg部分を使用
def growth_val_2(growth_mode, coverage, num_1ML, num_2ML):
    First_at_1ML = growth_mode[num_1ML][1]
    Second_at_2ML = growth_mode[num_2ML][2]
    if growth_mode[num_2ML][3] == 100:
        First_at_second = 0
    else:
        First_at_second = (
            growth_mode[num_2ML][0] / (100 - growth_mode[num_2ML][3]) * 100
        )
    if First_at_1ML < 50:
        if Second_at_2ML < 50:
            mode = "VW"
        else:
            mode = "BL"
    else:
        if Second_at_2ML >= 50:
            mode = "FM"
        else:
            if First_at_second < 50:
                mode = "SK"
            else:
                mode = "DW"

    return [
        First_at_1ML,
        Second_at_2ML,
        First_at_second,
        mode,
        coverage[num_1ML],
        coverage[num_2ML],
    ]


def mode_check_2(growth_mode, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    return growth_val_2(growth_mode, coverage, num_1ML, num_2ML)


#
#
#
def get_theta1(occupation):
    return occupation[0] - occupation[1]


def get_theta2(occupation):
    return occupation[1] - occupation[2]


def get_theta3_1(occupation):
    if occupation[2] == 100:
        return 0
    else:
        return (occupation[0] - occupation[1]) / (100 - occupation[2]) * 100


def get_theta3_2(occupation):
    if occupation[1] == 100:
        return 0
    else:
        return (occupation[0] - occupation[1]) / (100 - occupation[1]) * 100


def get_theta3_3(occupation):
    if occupation[2] == 100:
        return 0
    else:
        return (100 - occupation[0]) / (100 - occupation[2]) * 100


def get_theta3_4(occupation):
    if occupation[1] == 100:
        return 0
    else:
        return (100 - occupation[0]) / (100 - occupation[1]) * 100


def mode_judge(t1, t2, t3):
    if t1 < 50 and t2 < 50:
        return "VW"
    elif t1 >= 50 and t2 >= 50:
        return "FM"
    elif t1 < 50 and t2 >= 50:
        return "BL"
    elif t1 >= 50 and t2 <= 50:
        if t3 >= 50:
            return "SK"
        else:
            return "DW"
    else:
        return "Error"


def mode_judge2(t1, t2, t3):
    if t1 < 50 and t2 < 50:
        return "VW"
    elif t1 >= 50 and t2 >= 50:
        return "FM"
    elif t1 < 50 and t2 >= 50:
        return "BL"
    elif t1 >= 50 and t2 <= 50:
        if t3 < 50:
            return "SK"
        else:
            return "DW"
    else:
        return "Error"


#
#
# theta3で分母に3層目を使用
def mode_check_3(occupation, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    theta_1 = get_theta1(occupation[num_1ML])
    theta_2 = get_theta2(occupation[num_2ML])
    theta_3 = get_theta3_1(occupation[num_2ML])
    mode = mode_judge(theta_1, theta_2, theta_3)
    return [theta_1, theta_2, theta_3, mode]


#
#
# theta3で分母に2層目を使用
def mode_check_4(occupation, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    theta_1 = get_theta1(occupation[num_1ML])
    theta_2 = get_theta2(occupation[num_2ML])
    theta_3 = get_theta3_2(occupation[num_2ML])
    mode = mode_judge(theta_1, theta_2, theta_3)
    return [theta_1, theta_2, theta_3, mode]


#
#
# theta3で分母に3層目を使用、分子にAg領域を使用
def mode_check_5(occupation, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    theta_1 = get_theta1(occupation[num_1ML])
    theta_2 = get_theta2(occupation[num_2ML])
    theta_3 = get_theta3_3(occupation[num_2ML])
    mode = mode_judge2(theta_1, theta_2, theta_3)
    return [theta_1, theta_2, theta_3, mode]


#
#
# theta3で分母に2層目を使用、分子にAg領域を使用
def mode_check_6(occupation, coverage):
    num_1ML = num_check(coverage, 1)
    num_2ML = num_check(coverage, 2)
    theta_1 = get_theta1(occupation[num_1ML])
    theta_2 = get_theta2(occupation[num_2ML])
    theta_3 = get_theta3_4(occupation[num_2ML])
    mode = mode_judge2(theta_1, theta_2, theta_3)
    return [theta_1, theta_2, theta_3, mode]
