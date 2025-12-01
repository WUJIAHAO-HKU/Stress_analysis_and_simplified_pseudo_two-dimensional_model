import numpy as np
import matplotlib.pyplot as plt

# 创建时间数组
t = np.linspace(0, 10, 1000)

# 创建基础脉冲信号
def generate_pulse_signal(t, positions, amplitudes, base_signal=True):
    pulse_signal = np.zeros_like(t)
    for i in range(len(positions)):
        position = positions[i]
        amplitude = amplitudes[i]
        pulse_signal += np.exp(-((t - position) ** 2) * 10) * amplitude
    if base_signal:
        pulse_signal += 0.3 * np.sin(2 * np.pi * t)  # 添加基础正弦波
    return pulse_signal

# 生成第一条脉冲信号
pulse_positions1 = np.linspace(0, 9, 10)
pulse_amplitudes1 = [1.5, 1.8, 1.65, 1.95, 1.5, 1.8, 1.65, 1.95, 1.5, 1.65]
pulse_signal1 = generate_pulse_signal(t, pulse_positions1, pulse_amplitudes1, base_signal=True)

# 生成第二条脉冲信号（不添加基础正弦波）
pulse_positions2 = np.linspace(0.5, 9.5, 10)
pulse_amplitudes2 = [1.2, 1.5, 1.35, 1.65, 1.2, 1.5, 1.35, 1.65, 1.2, 1.5]
pulse_signal2 = generate_pulse_signal(t, pulse_positions2, pulse_amplitudes2, base_signal=False)

# 创建图形
plt.figure(figsize=(12, 6))

# 绘制第一条曲线
plt.plot(t, pulse_signal1, color='blue', linewidth=1.5, label='Pulse Signal 1')

# 绘制第二条曲线
plt.plot(t, pulse_signal2, color='red', linewidth=1.5, label='Pulse Signal 2')

# 添加图例
plt.legend(fontsize=12)

# 设置标题和坐标轴标签
plt.title('Pulse Waveforms with Spikes', fontsize=16)
plt.xlabel('Time', fontsize=12)
plt.ylabel('Amplitude', fontsize=12)

# 添加网格
plt.grid(True, linestyle='--', alpha=0.5)

# 调整坐标轴范围
plt.ylim(0, 2.5)

# 保存图像
plt.savefig('pulse_waveforms.png', dpi=300)
plt.show()