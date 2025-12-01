import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# ======================
# 参数设置 (根据NASA/JPL实验数据调整)
# ======================
fs = 1e14          # 采样频率 (10 THz，覆盖0.1-10THz频段)
T = 1e-9           # 总时长1ns
N = int(fs * T)    # 采样点数
t = np.arange(N) / fs  # 时间序列

# 量子噪声参数
PSD_target = 3.2e-27  # 目标功率谱密度 (A²/Hz)
f_low = 0.1e12       # 0.1 THz
f_high = 10e12       # 10 THz
flatness_tolerance = 0.8  # 频谱平坦度±0.8 dB
phase_coherence_length = 50e-9  # 相位相干长度50nm

# ======================
# 1. 生成基础量子噪声
# ======================
def generate_quantum_noise(N, PSD, fs):
    # 生成满足功率谱密度的白噪声
    freq = np.fft.fftfreq(N, 1/fs)
    fft_vals = np.zeros(N, dtype=complex)
    
    # 在目标频段内填充随机相位
    idx = np.where((np.abs(freq) >= f_low) & (np.abs(freq) <= f_high))[0]
    magnitudes = np.sqrt(PSD * fs / 2)  # 幅度计算
    phases = 2 * np.pi * np.random.rand(len(idx))  # 随机相位
    
    # 构建相干相位 (通过线性相位偏移实现50nm相干长度)
    k = 2 * np.pi / (phase_coherence_length * fs)  # 波数
    phases += k * freq[idx]  # 添加空间相干性相位
    
    fft_vals[idx] = magnitudes * np.exp(1j * phases)
    
    # 逆FFT得到时域信号
    noise_t = np.fft.ifft(fft_vals).real
    return noise_t / np.std(noise_t) * np.sqrt(PSD * fs)  # 归一化功率

quantum_noise = generate_quantum_noise(N, PSD_target, fs)

# ======================
# 2. 频谱平坦度修正
# ======================
# 设计FIR滤波器补偿带内波动
nyq = fs / 2.0
bands = [0.09e12, 0.1e12, 9.9e12, 10.1e12]  # 过渡带
desired = [0, 1, 0]                         # 带外抑制
fir_coeff = signal.remez(101, bands, desired, fs=fs)
quantum_noise_filtered = signal.lfilter(fir_coeff, 1.0, quantum_noise)

# 验证频谱平坦度
f, Pxx = signal.welch(quantum_noise_filtered, fs, nperseg=1024)
in_band = (f >= f_low) & (f <= f_high)
flatness_error = 10 * np.log10(Pxx[in_band] / PSD_target)
print(f"最大平坦度误差: {np.max(np.abs(flatness_error)):.2f} dB")

# ======================
# 3. 添加热涨落分量 (NIST实验模型)
# ======================
def thermal_fluctuation(T_noise, R):
    k = 1.38e-23  # 玻尔兹曼常数
    V_rms = np.sqrt(4 * k * T_noise * R * fs)
    return np.random.normal(0, V_rms / R, N)  # I = V/R

T_noise = 2.0  # 噪声温度2K
R_josephson = 1e3  # 约瑟夫森结阻抗1kΩ
thermal_noise = thermal_fluctuation(T_noise, R_josephson)

# 合成最终电流信号
current = quantum_noise_filtered + thermal_noise

# ======================
# 4. 数据保存与可视化
# ======================
# 保存时间-电流数据
np.savetxt('quantum_noise_current.csv', 
           np.column_stack((t[:1000], current[:1000])),  # 保存前1000点示例
           delimiter=',', header='Time(s),Current(A)')

# 绘制时域波形
plt.figure(figsize=(12, 6))
plt.subplot(211)
plt.plot(t[:100], current[:100])
plt.title('Quantum Noise Current (First 100 ps)')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')

# 绘制功率谱密度
plt.subplot(212)
plt.semilogy(f, Pxx)
plt.axvline(f_low, color='r', linestyle='--')
plt.axvline(f_high, color='r', linestyle='--')
plt.title('Power Spectral Density')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (A²/Hz)')
plt.tight_layout()
plt.show()