import random


def judge_null(success: float) -> bool:
    # assert success <= 1, success
    if success > 1:
        success = 1
    null = 1 - success
    result = [True, False]
    weight = [success, null]
    judge = random.choices(result, weights=weight)
    return judge[0]


if __name__ == "__main__":
    success = 0.7
    success_i = 0
    null_i = 0
    repetition = 100
    for _ in range(repetition):
        judge = judge_null(success)
        if judge:
            success_i += 1
        else:
            null_i += 1

    print("Test : judge null event")
    print("success probability = " + str(success))
    print("repetition = " + str(repetition))
    print("success times = " + str(success_i))
    print("null times = " + str(null_i))
