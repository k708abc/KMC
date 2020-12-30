import time

num = 10000
# method1
start = time.time()
list_1 = []
list_2 = []

for i in range(num):
    list_1.append(i)
    list_2.append(i)
time_1 = time.time() - start

# method2
start = time.time()
list_1 = [i for i in range(num)]
list_2 = [i for i in range(num)]
time_2 = time.time() - start
# method3
start = time.time()
list_3 = [[i, i] for i in range(num)]
time_3 = time.time() - start

print("method 1: " + str(time_1))
print("method 2: " + str(time_2))
print("method 3: " + str(time_3))