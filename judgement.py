import random


def judge_null(success):
    if success > 1:
        print("judge over 1")
    null = 1 - success
    result = ["success", "failure"]
    weight = (success, null)
    judge = random.choices(result, weights=weight)
    return judge
