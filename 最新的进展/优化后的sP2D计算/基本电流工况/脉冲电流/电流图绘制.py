import numpy as np
import matplotlib.pyplot as plt

# 定义参数
period1 = 30         # 脉冲1周期，单位：秒
duty_cycle1 = 0.1    # 脉冲1占空比，10%（假设110%为错误）
peak_current1 = 4    # 脉冲1峰值电流，单位：C
period2 = 10         # 脉冲2周期，单位：秒
duty_cycle2 = 0.5    # 脉冲2占空比，50%
peak_current2 = 2    # 脉冲2峰值电流，单位：C
total_time = 200     # 总时间，单位：秒
time_step = 0.1      # 时间步长，单位：秒

# 生成时间数组
time = np.arange(0, total_time, time_step)

# 计算脉冲电流1
position1 = time % period1
on_time1 = duty_cycle1 * period1  # 通电时间：3秒
current1 = np.where(position1 < on_time1, peak_current1, 0)

# 计算脉冲电流2
position2 = time % period2
on_time2 = duty_cycle2 * period2  # 通电时间：5秒
current2 = np.where(position2 < on_time2, peak_current2, 0)

# 绘制两个脉冲电流
plt.step(time, current1, where='post', color='tan', label='Pulse 1: period=30s, duty=10%, peak=4C')
plt.step(time, current2, where='post', color='slateblue', label='Pulse 2: period=10s, duty=50%, peak=2C')

# 设置图形属性
plt.xlabel('Time (s)')               # x轴标签
plt.ylabel('Current (C)')            # y轴标签
plt.title('Pulse Current Profiles')  # 图标题
plt.xlim(0, total_time)             # x轴范围
plt.ylim(0, 5)                      # y轴范围，略高于最大值4C
plt.grid(True)                      # 添加网格线
plt.legend()                        # 添加图例

# 显示图形
plt.show()
