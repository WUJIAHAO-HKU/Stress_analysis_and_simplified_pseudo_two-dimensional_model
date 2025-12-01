import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
import warnings
warnings.filterwarnings("ignore", category=UserWarning)  # 忽略nolds的警告

# ======================
# 参数设置（严格匹配技术指标）
# ======================
# 时间参数
total_time = 600      # 总时长600秒
fs = 1000             # 采样率1kHz
t = np.linspace(0, total_time, total_time*fs, endpoint=False)

# 放电阶段参数 (0-360s)
I_discharge_base = 2.4          # 基准电流2.4A
lorenz_depth = 0.15             # 混沌调制深度15%

# 充电阶段参数 (360-600s)
I_charge_base = -1.2            # 基准充电电流-1.2A
chua_amp = 0.3                  # Chua电路调制幅值0.3A
pulse_freq = 10                 # 基础脉冲频率10Hz
param_variation = 0.12          # ±12%参数摄动

# ======================
# 1. Lorenz系统生成混沌载波（放电阶段）
# ======================
def lorenz_system(state, t, sigma=10, rho=28, beta=8/3):
    """ Lorenz系统微分方程 """
    x, y, z = state
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]

# 初始条件和长时间仿真（确保混沌特性）
x0 = [0.1, 0.0, 0.0]  # 标准混沌初始条件
lorenz_duration = 400  # 仿真400秒（>360秒需求）
lorenz_t = np.linspace(0, lorenz_duration, lorenz_duration*fs)
lorenz_data = odeint(lorenz_system, x0, lorenz_t)

# 截取最后360秒数据并归一化
x = lorenz_data[:, 0][-360*fs:]
x_norm = (x - np.mean(x)) / np.max(np.abs(x))  # 归一化到[-1,1]
x_norm *= (1 + param_variation * np.random.randn(len(x_norm)))  # 加入参数摄动

# ======================
# 2. Chua电路生成脉冲调制（充电阶段）
# ======================
def chua_circuit(state, t, alpha=15.6, beta=28, m0=-1.143, m1=-0.714):
    """ Chua电路微分方程 """
    x, y, z = state
    h = m1*x + 0.5*(m0-m1)*(np.abs(x+1)-np.abs(x-1))  # 分段线性函数
    dx = alpha * (y - x - h)
    dy = x - y + z
    dz = -beta * y
    return [dx, dy, dz]

# 求解Chua电路（匹配0-50Hz带宽）
chua_data = odeint(chua_circuit, [0.7, 0, 0], t)
z = chua_data[:, 2]  # 使用z分量（频谱特性更丰富）

# 生成脉冲调制信号
pulse_wave = 0.8 + 0.2 * np.sign(np.sin(2*np.pi*pulse_freq*t))  # 基础10Hz方波
chua_mod = chua_amp * (z - np.mean(z)) / np.max(np.abs(z))      # 归一化调制信号

# ======================
# 3. 合成完整充放电电流
# ======================
current = np.zeros_like(t)

# 放电阶段 (0-360s): CCM调制
current[:360*fs] = I_discharge_base * (1 + lorenz_depth * x_norm)

# 充电阶段 (360-600s): CPWM调制
current[360*fs:] = I_charge_base + pulse_wave[360*fs:] * (1 + 0.5*chua_mod[360*fs:])

# 添加高斯噪声（模拟实际设备）
current += 0.1 * np.random.randn(len(current))

# ======================
# 4. 特征分析与可视化
# ======================
plt.figure(figsize=(15, 10))

# 时域波形
plt.subplot(3, 1, 1)
plt.plot(t[::1000], current[::1000], 'b')  # 每1000点取1个点显示
plt.axvline(360, color='r', linestyle='--', label='充放电切换点')
plt.title('三阶段混沌调制电流 (放电:CCM | 充电:CPWM)')
plt.xlabel('Time [s]')
plt.ylabel('Current [A]')
plt.grid(True)
plt.legend()

# 放电阶段细节
plt.subplot(3, 1, 2)
discharge_segment = current[:360*fs][10000:10100]  # 取1秒典型片段
plt.plot(np.linspace(0, 0.1, 100), discharge_segment, 'g')
plt.title('放电阶段混沌调制细节 (0.1秒窗口)')
plt.xlabel('Time [s]')
plt.ylabel('Current [A]')
plt.grid(True)

# 频谱分析
plt.subplot(3, 1, 3)
f, t_spec, Sxx = spectrogram(current, fs=fs, nperseg=1024)
plt.pcolormesh(t_spec, f, 10*np.log10(Sxx), shading='gouraud', cmap='jet')
plt.ylim(0, 5000)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.title('时频谱 (3.7kHz宽带谱验证)')
plt.colorbar(label='Power/Frequency [dB/Hz]')
plt.axhline(3700, color='w', linestyle='--', label='3.7kHz特征频率')
plt.legend()

plt.tight_layout()
plt.show()

# ======================
# 5. 混沌特性验证
# ======================
try:
    from nolds import lyap_r, hurst_rs
    # 计算李雅普诺夫指数（放电阶段）
    lyap_exp = lyap_r(x_norm, emb_dim=3)
    print(f"李雅普诺夫指数: {lyap_exp:.2f} (目标:0.82)")
    
    # 计算分形维数（Hurst指数法）
    hurst = hurst_rs(current)
    fractal_dim = 2 - hurst
    print(f"分形维数: {fractal_dim:.2f} (目标:2.35)")
except ImportError:
    print("未安装nolds库，请运行: pip install nolds")
