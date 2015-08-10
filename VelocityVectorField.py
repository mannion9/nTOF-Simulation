import math
import numpy as np
import matplotlib.pyplot as plt

def erf(x):
    return math.exp(x**2)

x = np.linspace(0,1,100)
y = [erf(i) for i in x]

plt.plot(x,y)
plt.show()
