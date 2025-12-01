import numpy as np
import matplotlib.pyplot as plt

# 参数设置
t = np.linspace(0, 10, 10000)  # 时间序列
f_carrier = 50                  # 载波频率（Hz）
chaos_amp = 0.5                 # 混沌信号幅度
carrier_amp = 1.0               # 载波幅度

# 1. 生成混沌信号（示例：简化蔡氏电路输出）
def generate_chaos(t):
    # 模拟混沌信号的类噪声特性
    chaos = np.cumsum(np.random.randn(len(t))) * 0.02
    chaos = chaos_amp * (np.sin(2 * np.pi * 0.8 * t) + 0.3 * np.sin(2 * np.pi * 2.3 * t) + chaos)
    return chaos

chaos_signal = generate_chaos(t)

# 2. 生成载波信号（常规正弦波）
carrier = carrier_amp * np.sin(2 * np.pi * f_carrier * t)

# 3. 混沌调制：将混沌信号与载波相乘
modulated_current = chaos_signal * carrier

# 绘制曲线
plt.figure(figsize=(12, 8))

# 混沌信号
plt.subplot(3, 1, 1)
plt.plot(t, chaos_signal, 'b', linewidth=0.8)
plt.title("Chaotic Signal (Modulating Signal)")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")

# 载波信号
plt.subplot(3, 1, 2)
plt.plot(t, carrier, 'g', linewidth=0.8)
plt.title("Carrier Signal (50 Hz Sine Wave)")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")

# 调制后的混沌电流
plt.subplot(3, 1, 3)
plt.plot(t, modulated_current, 'r', linewidth=0.8)
plt.title("Chaos-Modulated Current")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")

plt.tight_layout()
plt.show()