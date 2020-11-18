import math

def cal_rate(pre, kbt, E):
    rate = pre * math.exp(E/kbt)
    return rate