from typing import List


def normarize_rate(rates: List, n_value: float) -> List:
    n_rates = [i / n_value for i in rates]
    n_rates.append(1 - sum(n_rates))  # rate of null event
    return n_rates


if __name__ == "__main__":
    rates = [0.1, 10, 100]
    n_value = 1000
    print("Test: normarize list")
    n_rates = normarize_rate(rates, n_value)

    print("rates = " + str(rates))
    print("normarize = " + str(n_value))
    print("After normarize = " + str(n_rates))
    print("Sum = " + str(sum(n_rates)))
