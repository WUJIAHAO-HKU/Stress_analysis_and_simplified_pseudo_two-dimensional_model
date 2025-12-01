import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import sawtooth

# 定义参数
frequency = 1       # 频率，单位：Hz
amplitude = 1       # 幅值，单位：C
total_time = 3      # 总时间，单位：秒，显示3个周期
time_step = 0.01    # 时间步长，单位：秒

# 生成时间数组
time = np.arange(0, total_time, time_step)

# 计算正弦电流（无直流偏置）
sine_current = amplitude * np.sin(2 * np.pi * frequency * time)

# 计算线性电流（三角波）
triangular_current = amplitude * sawtooth(2 * np.pi * frequency * time, width=0.5)

# 绘制图形
plt.figure(figsize=(10, 6))  # 设置图形大小
plt.plot(time, sine_current, label='Sine waveform', color='#FFB6C1', linewidth=2)
plt.plot(time, triangular_current, label='Linear waveform (triangular wave)', color='teal', linewidth=2)

# 设置图形属性
plt.xlabel('time (s)', fontsize=14)         # x轴标签
plt.ylabel('current(C)', fontsize=14)         # y轴标签
plt.title('Schematic diagram of alternating current', fontsize=16)    # 图标题
plt.xlim(0, total_time)                    # x轴范围
plt.ylim(-1.5, 1.5)                        # y轴范围，略大于幅值
plt.grid(True)                             # 添加网格线
plt.legend(fontsize=12)                    # 添加图例

# 显示图形
plt.show()
