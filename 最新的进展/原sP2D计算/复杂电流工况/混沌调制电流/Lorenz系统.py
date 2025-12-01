import numpy as np
import matplotlib.pyplot as plt

# ======================
# 1. 定义Lorenz系统微分方程
# ======================
def lorenz_system(t, state, sigma, rho, beta):
    x, y, z = state
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return np.array([dx, dy, dz])

# ======================
# 2. 四阶龙格-库塔法求解
# ======================
def runge_kutta_4(func, t_span, y0, dt, params):
    t = np.arange(t_span[0], t_span[1], dt)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(1, len(t)):
        k1 = func(t[i-1], y[i-1], *params) * dt
        k2 = func(t[i-1] + dt/2, y[i-1] + k1/2, *params) * dt
        k3 = func(t[i-1] + dt/2, y[i-1] + k2/2, *params) * dt
        k4 = func(t[i-1] + dt, y[i-1] + k3, *params) * dt
        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
    return t, y

# ======================
# 3. 参数设置与求解
# ======================
# Lorenz系统参数 (经典混沌参数)
sigma, rho, beta = 10.0, 28.0, 8.0/3.0

# 初始条件和时间设置
y0 = [1.0, 1.0, 1.0]  # 初始状态[x, y, z]
t_span = [0, 20]      # 时间范围
dt = 0.01             # 时间步长

# 求解Lorenz系统
t, states = runge_kutta_4(lorenz_system, t_span, y0, dt, (sigma, rho, beta))
x, y, z = states.T  # 解包状态变量

# ======================
# 4. 混沌调制电流生成
# ======================
# 载波信号参数
f_carrier = 50          # 载波频率 (Hz)
carrier_amp = 1.0       # 载波幅度

# 生成载波信号
carrier = carrier_amp * np.sin(2 * np.pi * f_carrier * t)

# 调制过程：用Lorenz系统的x分量调制载波幅度
modulated_current = x * carrier

# ======================
# 5. 可视化
# ======================
plt.figure(figsize=(8, 7))



# (3) 混沌调制电流

plt.plot(t, modulated_current, 'r', linewidth=0.8)
plt.title("Lorenz-Modulated Current")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")

plt.tight_layout()
plt.show()