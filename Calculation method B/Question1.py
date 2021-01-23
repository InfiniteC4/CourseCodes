import csv
import numpy as np
import matplotlib.pyplot as plt

with open('data\\sea2020.csv', 'r') as f:
    reader = csv.reader(f)
    x = []
    y = []
    for row in reader:
        x.append(row[0])
        y.append(row[1])
f.close()
x[0] = 0.0
y[0] = 0.0
# del x[0]
# del y[0]
x = np.array(x, dtype=float)
y = np.array(y, dtype=float)

# 题目条件
N = 51

# 计算三次样条的参数Mi和d
M = np.zeros(N + 1)
M = y
d = np.zeros(N + 1)
c = np.zeros(N + 1)
b = np.zeros(N + 1)
a = np.zeros(N + 1)
h = np.zeros(N + 1)
for k in range(1, 3):
    for i in range(N, k + 2, -1):
        M[i] = (M[i] - M[i - 1]) / (x[i] - x[i - k])
h[1] = x[2] - x[1]
for i in range(2, N):
    h[i] = x[i + 1] - x[i]
    c[i] = h[i] / (h[i] + h[i - 1])
    a[i] = 1 - c[i]
    b[i] = 2
    d[i] = 6 * M[i + 1]

# 边界条件1得到
d[1] = 0
d[N] = 0
lam = 0
miu = 0
M[1] = 0
M[N] = 0

c[1] = lam
b[1] = 2
a[N] = miu
b[N] = 2

# TSS算法,解三对角线性方程组
u = np.zeros(N + 1)
ty = np.zeros(N + 1)
s = np.zeros(N + 1)
u[1] = b[1]
ty[1] = d[1]
for i in range(2, N + 1):
    s[i] = a[i] / u[i - 1]
    u[i] = b[i] - s[i] * c[i - 1]
    ty[i] = d[i] - s[i] * ty[i - 1]
M[N] = ty[N] / u[N]
for i in range(N - 1, 0, -1):
    M[i] = (ty[i] - c[i] * M[i + 1]) / u[i]

# 找出tx所在的区间(x[k-1]<=tx<=x[k])，以x[k]下标值K作为结果
# 用数据点{(x[i],y[i])}作三次插值样条函数,然后求x处的值y=s(x)
tx = np.arange(0, 5000, 10, dtype=float)
tn = len(tx)
s = np.zeros(tn)
for i in range(1, tn):
    k = 1
    for j in range(2, N):
        if tx[i] <= x[j]:
            k = j - 1
            break
        else:
            k = j
    H = x[k + 1] - x[k]
    x1 = x[k + 1] - tx[i - 1]
    x2 = tx[i - 1] - x[k]
    s[i - 1] = (M[k] * x1**3 / 6 + M[k + 1] * x2**3 / 6 +
                (y[k] - M[k] * H * H) * x1 +
                (y[k + 1] - M[k + 1] * H * H / 6) * x2) / H

print(M)

plt.plot(tx, s)
plt.ylabel('深度/m')
plt.xlabel('探测点/相距m')
plt.show()
