import csv
import numpy as np
import matplotlib.pyplot as plt


def DD(x=[], y=[], n=int):
    # 给定插值点,计算差商表d d顺序为y[n-1,n] y[n-2,n-1] ... y[0,1] y[n-2,n-1,n] ... y[0,1,...,n]
    k = 0
    d = []
    while (k < n):
        k += 1
        i = n + 1
        while (i > k):
            i -= 1
            y[i] = (y[i] - y[i - 1]) / (x[i] - x[i - k])
            d.append(y[i])
    return d


def TSS(a=[], b=[], c=[], d=[], n=int, x=[]):
    # 解三对角线性方程组,系数存于a,b,c中,d是右端向量
    u = [0] * (n + 1)
    y = [0] * (n + 1)
    s = [0] * (n + 1)
    u[0] = b[0]
    y[0] = d[0]
    for i in range(1, n):
        s[i] = a[i - 1] / u[i - 1]
        u[i] = b[i - 1] - s[i] * c[i - 1]
        y[i] = d[i] - s[i] * y[i - 1]
    x[n - 1] = y[n - 1] / u[n - 1]
    k = n - 1
    while (k > 0):
        k -= 1
        x[k] = (y[k] - c[k] * x[k + 1]) / u[k]
    return x


def FINDK(x=[], y=[], n=int, xb=float):
    # 找出xb所在的区间(x[k-1]<=xb<=x[k])，以x[k]下标值K作为结果
    K = 1
    for i in range(1, n):
        if xb <= x[i]:
            K = i
            break
        else:
            K = i + 1
    return K


def SPLINEM(x=[], y=[], M=[], n=int, lam=float, d0=float, miu=float, dn=float):
    # 计算三次样条的参数Mi,存放于M中，参数lam,d0,miu,dn由边界条件给定
    for i in range(0, n + 1):
        M[i] = y[i]
    for k in range(1, 3):
        for i in range(n, k - 1, -1):
            M[i] = (M[i] - M[i - 1]) / (x[i] - x[i - k])
    c = [0.0] * (n + 1)
    b = [0.0] * (n + 1)
    a = [0.0] * (n + 1)
    h = [0.0] * (n + 1)
    h[1] = x[1] - x[0]
    for i in range(1, n):
        h[i + 1] = x[i + 1] - x[i]
        c[i] = h[i + 1] / (h[i] + h[i + 1])
        a[i] = 1 - c[i]
        b[i] = 2
        M[i] = 6 * M[i + 1]
    M[0] = d0
    M[n] = dn
    c[0] = lam
    b[0] = 2
    a[n] = miu
    b[n] = 2
    S = TSS(a, b, c, M, n, M)
    return S


def EVASPLINE(x=[],
              y=[],
              n=int,
              xb=float,
              lam=float,
              d0=float,
              miu=float,
              dn=float):
    # 用数据点{(x[i],y[i])}作三次插值样条函数,然后求xb处的值yb=s(xb)
    M = [0.0] * (n + 1)
    k = 1
    M = np.array(SPLINEM(x, y, M, n, lam, d0, miu, dn))
    print(M)
    k = np.array(FINDK(x, y, n, xb))
    h = np.array(x[k] - x[k - 1])
    xbb = np.array(x[k] - xb)
    xba = np.array(xb - x[k - 1])
    yb = [
        M[k - 1] * xbb * xbb * xbb / 6 +
        (y[k - 1] - M[k - 1] * h * h / 6) * xbb +
        (y[k] - M[k] * h * h / 6) * xba
    ] / h
    return yb


with open('data\\sea2020.csv', 'r') as f:
    reader = csv.reader(f)
    x = []
    y = []
    for row in reader:
        x.append(row[0])
        y.append(row[1])
f.close()
del x[0]
del y[0]
x = list(map(int, x))
y = list(map(float, y))

# 边界条件1得到
d0 = 0
d50 = 0
lam = 0
miu = 0
M0 = 0
M50 = 0

d = DD(x, y, 50)

ybb = [0.0] * 51
ybb = np.array(ybb)

for i in range(0, 51):
    ybb[i] = EVASPLINE(x, y, 50, x[i], lam, d0, miu, d50)

# plt.plot([0,1,2,3,4])  # 方式1：list参数，默认为y值，自动生成x轴值
plt.plot(x, ybb)  # 方式2：list参数，传入x和y值,
plt.ylabel('深度/m')
plt.xlabel('探测点/相距m')
plt.show()
