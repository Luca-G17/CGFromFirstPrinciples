import psutil
import numpy as np
import matplotlib.pyplot as plt

for i in range(10000):
    t = i
    temp = psutil.sensors_temperatures()['k10temp'][0].current
    plt.scatter(t, temp, c="red")
    plt.pause(1)

plt.show()

