import random


def judge_null(success: float) -> str:
    if success > 1:
        print("Error: probablity over 1")
    null = 1 - success
    result = ["success", "failure"]
    weight = [success, null]
    judge = random.choices(result, weights=weight)
    return judge[0]
