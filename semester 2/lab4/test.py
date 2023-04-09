import numpy as np


x = np.linspace(0,1,10)
y = np.linspace(0,1,20)
X, Y = np.meshgrid(x, y)
print(X)
print(Y)

print(X+Y)

