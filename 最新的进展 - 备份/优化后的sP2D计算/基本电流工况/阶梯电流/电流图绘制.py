import matplotlib.pyplot as plt

# 定义时间点和三条电流值
time1 = [0, 100, 200, 300]  # 时间点 (秒)，总时长300秒
time2 = [0, 70, 180, 300]  # 时间点 (秒)，总时长300秒
time3 = [0, 110, 230, 300]  # 时间点 (秒)，总时长300秒

# 第一条：越来越高
current1 = [1, 2, 3, 3.8]    # 0-100秒: 1C, 100-200秒: 2C, 200-300秒: 3C

# 第二条：越来越低
current2 = [3.6, 3.2, 2.5, 1.3]    # 0-100秒: 3C, 100-200秒: 2C, 200-300秒: 1C

# 第三条：忽高忽低
current3 = [0.7, 2.8, 1.8, 2.2]    # 0-100秒: 1C, 100-200秒: 3C, 200-300秒: 1C

# 绘制三条阶梯电流
plt.step(time1, current1, where='post', color='orange', linewidth=2, label='Increasing (1C→3C)')
plt.step(time2, current2, where='post', color='indianred', linewidth=2, label='Decreasing (3C→1C)')
plt.step(time3, current3, where='post', color='mediumslateblue', linewidth=2, label='Alternating (1C↔3C)')

# 设置图形属性
plt.xlabel('time (s)', fontsize=14)         # x轴标签
plt.ylabel('current (C)', fontsize=14)         # y轴标签
plt.title('Schematic diagram of the stepped current', fontsize=16)    # 图标题
plt.xlim(0, 300)                          # x轴范围
plt.ylim(0, 4)                            # y轴范围，略高于最大值3C
plt.grid(True)                            # 添加网格线
plt.legend(fontsize=12)                   # 添加图例

# 显示图形
plt.show()
