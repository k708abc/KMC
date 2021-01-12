import time
import math
import numpy as np

num = 100000
start = time.time()
value = 0.2225

for i in range(num):
    k = math.exp(value)

end_time = time.time() - start
print(end_time)

#
start = time.time()

for i in range(num):
    k = np.exp(value)

end_time = time.time() - start
print(end_time)
