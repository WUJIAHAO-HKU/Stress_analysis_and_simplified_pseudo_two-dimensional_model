using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 充放电电流参数
charging_current = 4.0  # 充电电流（安培）
discharging_current = -2.0  # 放电电流（安培）
charge_time = 2000  # 充电时间（秒）
discharge_time = 2000  # 放电时间（秒）
total_time = charge_time + discharge_time  # 总时间

# 脉冲电流参数
pulse_current = 3.0     # 峰值脉冲电流（安培）
low_current = 0.5       # 低电流（安培）
pulse_width = 200       # 脉冲宽度（秒）
pulse_interval = 1000   # 脉冲间隔（秒）

opt.dtType = "fixed"
opt.time = [0, total_time]  # 设置总模拟时间

# 设置电流随时间变化：充电和放电阶段的脉冲电流
opt.Current = (t) -> begin
    if t < charge_time
        # 充电阶段的脉冲电流
        if t % (pulse_width + pulse_interval) < pulse_width
            return charging_current + pulse_current  # 充电阶段加上脉冲电流
        else
            return charging_current + low_current  # 充电阶段加上低电流
        end
    else
        # 放电阶段的脉冲电流
        if t % (pulse_width + pulse_interval) < pulse_width
            return discharging_current + pulse_current  # 放电阶段加上脉冲电流
        else
            return discharging_current + low_current  # 放电阶段加上低电流
        end
    end
end

# 运行模型模拟
for model in ["sP2D", "P2D"]
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # 绘制结果
    plot!(result["time [s]"], result["cell voltage [V]"], label=model, xlabel="time [s]", ylabel="cell voltage [V]", lw=2)
end

# 设置x轴显示范围为total_time
xlims!(0, total_time)

# 保存图形
savefig("pulse_charge_discharge_current_validation.pdf")

# 引用
JuBat.Citation(["ai2024b"])