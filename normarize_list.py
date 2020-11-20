def normarize_rate(rates, n_value):
    n_rates = [i / n_value for i in rates]
    n_rates.append(1 - sum(n_rates))
    return n_rates
