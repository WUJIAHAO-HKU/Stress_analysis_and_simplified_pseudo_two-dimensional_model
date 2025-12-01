import numpy as np
import matplotlib.pyplot as plt

# 定义参数
periods = [10, 15, 20]      # 周期列表，单位：秒
peaks = [0.8, 1.5, 2.4]         # 峰值电流列表，单位：C
total_time = 60             # 总时间，单位：秒
time_step = 0.1             # 时间步长，单位：秒

# 生成时间数组
time = np.arange(0, total_time, time_step)

# 计算三条锯齿电流
currents = [(time % period) / period * peak for period, peak in zip(periods, peaks)]

# 使用科研配色（Matplotlib默认颜色循环的前三种：蓝色、橙色、绿色）
colors = ['C0', 'C1', 'C2']

# 绘制三条锯齿电流
for i, current in enumerate(currents):
    plt.plot(time, current, color=colors[i], linewidth=2, 
             label=f'Period={periods[i]}s, Amplitude={peaks[i]}C')

# 设置图形属性
plt.xlabel('time(s)')              # x轴标签
plt.ylabel('current(c)')              # y轴标签
plt.title('Schematic diagram of sawtooth current')         # 图标题
plt.xlim(0, total_time)            # x轴范围
plt.ylim(0, 2.5)                   # y轴范围，略高于最大峰值2C
plt.grid(True)                     # 添加网格线
plt.legend()                       # 添加图例

plt.gcf().set_size_inches(16, 6) 

# 显示图形
plt.show()