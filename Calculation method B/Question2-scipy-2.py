import csv
import numpy.matlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

with open('data\\alltime_world.csv', 'r') as f:
    reader = csv.reader(f)
    No = []
    date = []
    confirm = []
    for row in reader:
        No.append(row[0])
        date.append(row[1])
        confirm.append(row[2])
f.close()
del No[0]
del date[0]
del confirm[0]
No = np.array(No, dtype=float) + 1.0
date = np.array(date, dtype=float)
confirm = np.array(confirm, dtype=float)

N = 2  # 拟合多项式的最高次数
M = 307  # 日期的个数


def error(p, x, y):
    return p[0] * x**2 + p[1] * x + p[2] - y


p0 = [5, 2, 10]

ret = leastsq(error, p0, args=(No, confirm))

a, b, c = ret[0]

x = np.arange(1, 307, 0.1, dtype=float)
y = a * x**2 + b * x + c

plt.plot(No, confirm, 'o', label="样本点")
plt.plot(x, y, label="插值点")
plt.ylabel('CONFIRM')
plt.xlabel('DATE')
plt.show()
