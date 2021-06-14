# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

def general_erf(t, mean, std_dev):
    return (erf((t - mean) / np.sqrt(2 * std_dev**2)) - erf((- mean) / np.sqrt(2 * std_dev**2))) / 2

mean = 800
std_dev = 1

x = np.linspace(mean - 3*std_dev, mean + 3*std_dev, 1000)
y = np.array([1 - general_erf(xi, mean, std_dev) for xi in x])

plt.clf()
plt.style.use('seaborn-whitegrid')
plt.plot(x, y)
plt.grid(linestyle="--")
plt.show()