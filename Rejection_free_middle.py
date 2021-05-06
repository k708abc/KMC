import pickle
from Rejection_free_kmc import rejection_free

rf_class = pickle.load(open("selfdata.pickle", "rb"))
print("pickle read")
rf_class.init_value.start_from_middle = True
rf_class.start_from_middle()
