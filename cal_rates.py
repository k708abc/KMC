import math


def rate(pre: float, kbt: float, energy: float) -> float:
    return pre * math.exp(energy / kbt)
