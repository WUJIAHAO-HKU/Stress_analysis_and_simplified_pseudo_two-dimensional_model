import numpy as np
import matplotlib.pyplot as plt

# 创建时间数组
t = np.linspace(0, 10, 1000)

# 创建多级脉冲信号
signal = np.zeros_like(t)

# 设置多级脉冲的幅度和位置
amplitudes = [0.5, 1.0, 0.7, 1.2, 0.9, 1.1, 0.6, 1.3, 0.8, 1.0]
positions = np.linspace(0, 9, len(amplitudes))

# 生成多级脉冲信号
for i in range(len(amplitudes)):
    start = positions[i]
    end = positions[i] + 0.5  # 每个脉冲持续0.5秒
    idx = (t >= start) & (t < end)
    signal[idx] = amplitudes[i]

    # 添加过渡
    if i < len(amplitudes) - 1:
        next_amp = amplitudes[i + 1]
        transition_idx = (t >= end) & (t < end + 0.1)
        signal[transition_idx] = np.linspace(amplitudes[i], next_amp, np.sum(transition_idx))

# 添加噪声
noise = np.random.normal(0, 0.05, len(t))
signal += noise

# 创建图形
plt.figure(figsize=(12, 6))
plt.plot(t, signal, color='blue', linewidth=1.5)
plt.title('Multi-Level Pulse Waveform', fontsize=16)
plt.xlabel('Time', fontsize=12)
plt.ylabel('Amplitude', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.axis('tight')
plt.tight_layout()

# 保存图像
plt.savefig('multi_level_pulse_waveform.png', dpi=300)
plt.show()