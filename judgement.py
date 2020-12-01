import random


def judge_null(success: float) -> str:
    if success > 1:
        print("Error: probablity over 1")
    null = 1 - success
    result = ["success", "failure"]
    weight = [success, null]
    judge = random.choices(result, weights=weight)
    return judge[0]


if __name__ == "__main__":
    success = 1.2
    success_i = 0
    null_i = 0
    for _ in range(100):
        judge = judge_null(success)
        if judge == "success":
            success_i += 1
        elif judge == "failure":
            null_i += 1

    print("judge null test")
    print("success probability = " + str(success))
    print("success times = " + str(success_i))
    print("null times = " + str(null_i))
