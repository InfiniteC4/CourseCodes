import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

with open('data\\sea2020.csv', 'r') as f:
    reader = csv.reader(f)
    x = []
    y = []
    for row in reader:
        x.append(row[0])
        y.append(row[1])
f.close()
print(x)
print(y)
# x[0] = 0.0
# y[0] = 0.0
del x[0]
del y[0]
x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
print(x)
print(y)

new_x = np.arange(0, 5000, 10, dtype=float)

ipo3 = spi.splrep(x, y, k=3)
iy3 = spi.splev(new_x, ipo3)

plt.plot(x, y, 'o', label="样本点")
plt.plot(new_x, iy3, label="插值点")
plt.ylabel('Y')
plt.xlabel('X')
plt.show()
