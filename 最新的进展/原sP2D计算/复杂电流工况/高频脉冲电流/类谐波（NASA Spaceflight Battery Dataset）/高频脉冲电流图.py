import numpy as np
import matplotlib.pyplot as plt

# 创建时间数组
t = np.linspace(0, 10, 1000)

# 创建脉冲波形
pulse = np.zeros_like(t)
pulse[::50] = 1  # 在特定位置设置脉冲

# 创建一个类似于心电图的波形
ecg_signal = np.sin(2 * np.pi * t) * 0.5  # 基础正弦波
ecg_signal += np.sin(4 * np.pi * t) * 0.3  # 添加高频成分
ecg_signal += np.sin(8 * np.pi * t) * 0.2  # 添加更高频成分

# 在波形中加入尖峰（类似于心电图的QRS波群）
for i in range(10):
    peak_position = i * 1.0
    ecg_signal += np.exp(-((t - peak_position) ** 2) * 10) * 1.0

# 添加噪声
noise = np.random.normal(0, 0.1, len(t))
ecg_signal += noise

# 创建图形
plt.figure(figsize=(12, 6))
plt.plot(t, ecg_signal, color='blue', linewidth=1.5)
plt.title('Pulse Waveform', fontsize=16)
plt.xlabel('Time', fontsize=12)
plt.ylabel('Amplitude', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.axis('tight')
plt.tight_layout()

# 保存图像
plt.savefig('pulse_waveform.png', dpi=300)
plt.show()