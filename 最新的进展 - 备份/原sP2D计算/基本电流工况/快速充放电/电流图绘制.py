import matplotlib.pyplot as plt

# 定义时间点和电流值
time = [0, 0, 360,360,720, 720,1000,1000, 1440]  # 时间点 (秒)
current1 = [0, 5, 5,5,5, -5, -5,-5,-5]    # 第一条电流值 (单位：C)
current2 = [0, 3, 3,-3,-3, -3,-3,-3, -3]    # 第二条电流值 (单位：C)
current3 = [0,2,2, 2, 2, 2, 2,-2,-2]    # 第三条电流值 (单位：C)

# 绘制阶梯图
plt.step(time, current1, where='post', color='g', linewidth=2, label='Current 1 (5C)')
plt.step(time, current2, where='post', color='sandybrown', linewidth=2, label='Current 2 (3C)')
plt.step(time, current3, where='post', color='crimson', linewidth=3, label='Current 3 (2C)')

# 设置图形属性
plt.xlabel('Time (s)', fontsize=14)  # x轴标签，增大字体
plt.ylabel('Current (C)', fontsize=14)  # y轴标签，增大字体
plt.title('Fast Charge-Discharge Current Profile', fontsize=16)  # 图标题，增大字体
plt.xlim(0, 1440)               # x轴范围
plt.ylim(-6, 6)                 # y轴范围

# 增强x轴和y轴显示效果
plt.gca().spines['bottom'].set_linewidth(2)  # 加粗x轴线条
plt.gca().spines['left'].set_linewidth(2)    # 加粗y轴线条
plt.xticks(fontsize=12)                      # 增大x轴刻度字体
plt.yticks(fontsize=12)                      # 增大y轴刻度字体

# 添加2C和-2C的虚线
plt.axhline(y=2, color='r', linestyle='--', linewidth=1, label='2C')
plt.axhline(y=-2, color='r', linestyle='--', linewidth=1, label='-2C')

# 添加网格线和图例
plt.grid(True)
plt.legend()

plt.gcf().set_size_inches(16, 6) 

# 显示图形
plt.show()
