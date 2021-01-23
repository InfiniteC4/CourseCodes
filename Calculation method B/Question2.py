import csv
import math
import numpy.matlib
import numpy as np
import matplotlib.pyplot as plt

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

N = 8  # 拟合多项式的最高次数
M = 307
num = 0

# 生成G矩阵
G = np.matlib.ones((307, 1))
for i in range(1, N):
    g = No**i
    g = np.mat(g).T
    G = np.hstack((G, g))
G = np.hstack((G, np.mat(confirm).T))

# 生成Q矩阵
for i in range(0, N):
    if G[i, i] > 0:
        sgn = 1
    elif G[i, i] == 0:
        sgn = 0
    else:
        sgn = -1
    for j in range(i, M):
        num = G[j, i]**2 + num
    sigma = -sgn * math.sqrt(num)
    num = 0
    omega = np.zeros(N + 1)
    omega[i] = G[i, i] - sigma
    for k in range(i + 1, M):
        omega[k] = G[k, i]
    beta = sigma * omega[i]
    # 变换Gi-1到Gi
    G[i, i] = sigma
    for j in range(i + 1, N + 1):
        for k in range(i, M):
            num = omega[k] * G[k, j] + num
        x = num / beta
        num = 0
        for k in range(1, M):
            G[k, j] = G[k, j] + x * omega[k]

# 解三角方程求多项式系数X
X = np.zeros(N + 1)
X[N] = G[N, N + 1] / G[N, N]
for i in range(N - 1, 0, -1):
    for j in range(i + 1, N + 1):
        num = G[i, j] * X[j] + num
    X[i] = (G[i, N + 1] - num) / G[i, i]
    num = 0

# 计算误差
E = 0
for i in range(N + 1, M - 1):
    E = E + G[i, N + 2]**2

# 输出参数
y = np.zeros(M)
x = np.arange(1, M - 1, 0.1, dtype=float)
for i in range(0, M):
    for j in range(0, N + 1):
        y[i] = X[j] * x[i]**(j - 1) + y[i]

# plt.plot([0,1,2,3,4])  # 方式1：list参数，默认为y值，自动生成x轴值
plt.plot(date, y)  # 方式2：list参数，传入x和y值,
plt.ylabel('confirm')
plt.xlabel('date')
plt.show()
