import numpy as np

a = np.random.random((3, 3))

b = a.view()

b[0][0] = -100
print(a)
print(b)