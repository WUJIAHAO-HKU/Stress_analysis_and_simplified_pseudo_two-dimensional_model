import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import chirp

# 设置随机种子
np.random.seed(42)

# LG M50电池参数
rated_capacity = 1.5  # Ah
C_rate = 1.0  # 基准C率
max_C_rate = 4.0  # 最大允许C率 (基于规格书)

# 生成时间列 (0-149.9秒，间隔0.1秒)
time = np.arange(0, 150, 0.1)

def generate_normalized_current(time):
    # 基础信号组成
    base_signal = 0.5 * np.sin(2*np.pi*time/2)  # 低频成分 (5秒周期)
    
    # 添加不同频率的脉冲
    pulse1 = 0.3 * chirp(time, f0=0.5, f1=4, t1=max(time), method='linear')  # 扫频脉冲
    pulse2 = 0.2 * (time % 5 < 2) * np.random.rand(len(time))  # 间歇性随机脉冲
    
    # 组合信号
    combined = base_signal + pulse1 + pulse2
    
    # 归一化到C-rate范围 [-4C, 4C]
    normalized = combined * (2*max_C_rate/np.max(np.abs(combined)))
    
    # 转换为实际电流值 (A)
    current = normalized * rated_capacity
    
    # 添加高频噪声
    current += 0.3 * rated_capacity * np.random.normal(size=len(time))
    
    # 最终电流限制
    return np.clip(current, -max_C_rate*rated_capacity, max_C_rate*rated_capacity)

current = generate_normalized_current(time)

# 创建DataFrame
df = pd.DataFrame({
    'time(s)': np.round(time, 1),
    'current(A)': np.round(current, 3)
})

# 可视化验证
plt.figure(figsize=(12, 5))
plt.plot(df['time(s)'][:500], df['current(A)'][:500])
plt.axhline(y=rated_capacity*max_C_rate, color='r', linestyle='--', label='Max Discharge (4C)')
plt.axhline(y=-rated_capacity*max_C_rate, color='g', linestyle='--', label='Max Charge (4C)')
plt.title('Normalized Pulse Current for LG M50 Battery\n(First 50 seconds)')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# 保存为CSV
df.to_csv('LG_M50_normalized_current.csv', index=False)

# 数据统计信息
print("== 数据统计信息 ==")
print(f"最大放电电流: {df['current(A)'].max():.2f}A ({df['current(A)'].max()/rated_capacity:.1f}C)")
print(f"最大充电电流: {df['current(A)'].min():.2f}A ({df['current(A)'].min()/rated_capacity:.1f}C)")
print(f"平均绝对电流: {np.mean(np.abs(df['current(A)'])):.2f}A")
print(f"数据点数: {len(df)}")
