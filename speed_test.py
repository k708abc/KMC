import time

list_val = []
dict_val = {}

for i in range(0, 10000000):
    list_val.append(i)
    dict_val[i] = i

start_time_list = time.time()
val = 0
for v in list_val:
    val += v
elapsed_time_list = time.time() - start_time_list
#
start_time_list_sum = time.time()
val = sum(list_val)
elapsed_time_list_sum = time.time() - start_time_list_sum
#

#
start_time_dict = time.time()
val = 0
for _, v in dict_val.items():
    val += v
elapsed_time_dict = time.time() - start_time_dict
#
print("list_time:" + str(elapsed_time_list))
print("list_time_sum:" + str(elapsed_time_list_sum))
print("dict_time:" + str(elapsed_time_dict))
