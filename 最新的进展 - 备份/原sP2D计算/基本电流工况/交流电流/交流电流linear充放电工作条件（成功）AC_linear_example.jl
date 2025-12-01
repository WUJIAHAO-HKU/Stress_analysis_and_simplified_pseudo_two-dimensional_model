using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

# 定义电池参数
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 设置阶梯电流工作条件下的电流参数
current_steps = [4.0, 2.0, 0.0, -2.0, -4.0]  # 阶梯电流（安培）：充电阶段 -> 放电阶段
current_durations = [2000, 2000, 2000, 2000, 2000]  # 每个电流阶段的持续时间（秒）

# 设置模拟时间
total_time = sum(current_durations)  # 总时间

# 初始化模拟配置
opt.dtType = "fixed"
opt.time = [0, total_time]  # 设置总模拟时间

# 设置电流随时间的变化：阶梯电流充电与放电阶段
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
savefig("charge_discharge_step_comparison.pdf")
