using Plots, JuBat, CSV, DataFrames

# 定义电池参数
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 设置锯齿电流放电条件
max_discharging_current = -4.0  # 最大放电电流（安培）
min_discharging_current = -1.0  # 最小放电电流（安培）
pulse_duration = 0.1  # 脉冲持续时间（秒）
pulse_interval = 1.0  # 脉冲间隔（秒）
num_pulses = 10       # 脉冲数量

# 设置放电电流的锯齿波形（电流逐渐从最大值递减到最小值）
current_steps = [max_discharging_current, min_discharging_current]  # 放电电流的两个阶梯
current_durations = [pulse_duration * num_pulses, pulse_duration * num_pulses]  # 每个电流阶段的持续时间（秒）

# 设置模拟时间
total_time = sum(current_durations)  # 总时间

# 初始化模拟配置
opt.dtType = "fixed"
opt.time = [0, total_time]  # 设置总模拟时间

# 设置电流随时间的变化：锯齿电流放电阶段
opt.Current = (t) -> begin
    # 计算当前时间所在的阶段
    elapsed_time = 0
    for (i, duration) in enumerate(current_durations)
        elapsed_time += duration
        if t <= elapsed_time
            return current_steps[i]  # 返回当前阶段的电流
        end
    end
    return 0.0  # 默认情况下为零电流（理论上不会到达）
end

# 运行模型模拟并比较P2D和SP2D模型
for model in ["sP2D", "P2D"]
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # 绘制模型结果
    plot!(result["time [s]"], result["cell voltage [V]"], label=model, xlabel="时间 (秒)", ylabel="电池电压 (V)", lw=2)
end

# 设置x轴显示范围为total_time
xlims!(0, total_time)

# 保存图形
savefig("sawtooth_discharge_comparison.pdf")
