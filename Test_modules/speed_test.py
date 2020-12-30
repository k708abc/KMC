import time

start = time.time()
test_list = []
for i in range(100000):
    test_list.append(i)

for i in range(100000):
    test_list.append(i)

for i in range(100000):
    test_list.append(i)

time_for = time.time() - start

start = time.time()
test_1 = [i for i in range(100000)]
test_2 = [i for i in range(100000)]
test_3 = [i for i in range(100000)]

test_in = test_1 + test_2 + test_3

time_in = time.time() - start

#
start = time.time()
test_1 = [i for i in range(100000)]
test_1 += [i for i in range(100000)]
test_1 += [i for i in range(100000)]
time_new = time.time() - start

print("for :" + str(time_for))
print("in: " + str(time_in))
print("new: " + str(time_new))