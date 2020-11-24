from typing import List


def normarize_rate(rates: List, n_value: float) -> List:
    n_rates = [i / n_value for i in rates]
    n_rates.append(1 - sum(n_rates))  # rate of null event
    return n_rates
