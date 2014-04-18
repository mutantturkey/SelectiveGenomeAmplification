import time
import filter_max_consecutive_binding as fmax

fmax.test()

t = time.time()
for i in range(100000):
		fmax.max_consecutive_binding('ACCATATCA', 'TACATATA')
print time.time() - t
