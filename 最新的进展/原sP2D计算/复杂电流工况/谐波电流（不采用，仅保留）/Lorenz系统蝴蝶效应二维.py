import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 定义Lorenz系统的微分方程
def lorenz(state, t, sigma, rho, beta):
    x, y, z = state
    dx_dt = sigma * (y - x)
    dy_dt = x * (rho - z) - y
    dz_dt = x * y - beta * z
    return [dx_dt, dy_dt, dz_dt]

# 参数设置
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0

# 时间点
t = np.linspace(0, 50, 10000)

# 初始条件
state0 = [1.0, 1.0, 1.0]

# 数值积分求解Lorenz系统
states = odeint(lorenz, state0, t, args=(sigma, rho, beta))

# 提取x, y, z坐标
x = states[:, 0]
y = states[:, 1]
z = states[:, 2]

# 可视化
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

# 绘制轨迹（改动部分）
ax.scatter(x,z, c=t, cmap='hot', s=1, alpha=0.9)

# 设置长宽比
ax.set_aspect(0.5)  # 等比例，X轴和Z轴单位长度相等

# 添加网格
ax.grid(True, linestyle='--', color='gray', alpha=0.5)

# 设置标签
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_title('Lorenz Attractor')

plt.colorbar(ax.scatter(x,z, c=t, cmap='hot', s=0.5, alpha=0.9), label='Time',shrink=0.5,aspect =10)

# 显示图形
plt.show()
