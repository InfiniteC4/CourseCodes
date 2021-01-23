import struct
import numpy as np
np.set_printoptions(linewidth=400)  # 终端矩阵随窗口换行

# 文件列表
FName1 = "data\\data20201.dat"
FName2 = "data\\data20202.dat"
FName3 = "data\\data20203.dat"
FName4 = "data\\data20204.dat"

# 打开文件
FName = FName2
print("\n\n数据文件:", FName)
f = open(FName, 'rb')

# 读入文件头标志部分
id, ver, id1 = struct.unpack('<3I', f.read(12))

# 读入系数矩阵结构信息
n, q, p = struct.unpack('<3I', f.read(12))

typename = {0x102: '非压缩格式', 0x202: '压缩格式'}
print('id:', hex(id), '  ver:', hex(ver), '  id1:', hex(id1), "   ",
      typename[ver])
print('阶数:', n, ' 上带宽:', q, ' 下带宽:', p)

# 判断格式,读取系数矩阵
if ver == 0x102:
    A = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            A[i, j], = struct.unpack('f', f.read(4))
elif ver == 0x202:
    A = np.zeros([n, p + q + 1])
    for i in range(n):
        for j in range(p + q + 1):
            A[i, j], = struct.unpack('f', f.read(4))
    # 整理系数矩阵
    for i in range(p):
        for j in range(i + q + 1):
            A[i, j] = A[i, p - i + j]
            A[i, p - i + j] = 0.0
else:
    print("error: this ver is unknow!")

# 定义右端常量
B = np.zeros([n, 1])

# 读取右端常量
for i in range(n):
    B[i, 0], = struct.unpack('f', f.read(4))
f.close()

# 判断格式，Guess消元
if ver == 0x102:
    # 获取增广矩阵
    M = np.hstack((A, B))
    # 增广矩阵列数
    col = M.shape[1]
    # 消元
    for i in range(col - 2):
        # 计算系数
        p = M[i + 1:, i] / M[i, i]
        # 计算消元时减去的矩阵
        m = np.tile(M[i, :], (p.shape[0], 1)) * np.tile(p, (col, 1)).T
        # 消元
        M[i + 1:, :] = M[i + 1:, :] - m
    # 回代
    x = np.zeros(col - 1)
    for i in range(col - 2, -1, -1):
        x[i] = (M[i, -1] - np.dot(M[i, :-1], x.T)) / M[i, i]
elif ver == 0x202:
    for i in range(1, n):
        if i <= p:
            for j in range(i):
                B[i, 0] = B[i, 0] - B[j, 0] * A[i, 0] / A[j, 0]
                A[i, :] = A[i, :] - A[j, :] * A[i, 0] / A[j, 0]
                for u in range(p + q):
                    A[i, u] = A[i, u + 1]
                    A[i, u + 1] = 0.0
        else:
            for j in range(p - 1, -1, -1):
                B[i, 0] = B[i, 0] - B[i - j - 1, 0] * A[i, 0] / A[i - j - 1, 0]
                A[i, :] = A[i, :] - A[i - j - 1, :] * A[i, 0] / A[i - j - 1, 0]
                for u in range(p + q):
                    A[i, u] = A[i, u + 1]
                    A[i, u + 1] = 0.0
    print(A, '\n', B)
    # 回代
    x = np.zeros(n)
    x[n - 1] = B[n - 1] / A[n - 1, 0]
    for i in range(n - 2, -1, -1):
        if i >= n - q:
            for j in range(n - i - 1):
                B[i] = B[i] - x[i + 1 + j] * A[i, 1 + j]
        else:
            for j in range(q):
                B[i] = B[i] - x[i + 1 + j] * A[i, 1 + j]
        x[i] = B[i] / A[i, 0]
else:
    print("error: this ver is unknow!")

print(x)
