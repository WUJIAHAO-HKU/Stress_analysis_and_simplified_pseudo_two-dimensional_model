using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 锯齿波流放电参数
charge_amplitude = 4.0    # 放电电流（A）
charge_period = 500    # 放电锯齿波周期（秒）
charge_time = 2000        # 放电总时间（秒）

# 锯齿波充电参数
discharge_amplitude = -2.0 # 充电峰值电流绝对值（A）
discharge_period = 500    # 充电锯齿波周期（秒）
discharge_time = 2000     # 充电总时间（秒）

total_time = charge_time + discharge_time

opt.dtType = "fixed"
opt.time = [0, total_time]

# 定义复合电流函数
opt.Current = (t) -> begin
    if t < charge_time
        # 放电阶段：周期性下降锯齿波
        return charge_amplitude * (1 - (t % charge_period) / charge_period)
    else
        # 充电阶段：周期性下降锯齿波
        return discharge_amplitude * (1 - (t % discharge_period) / discharge_period)
    end
end

# 初始化绘图并设置网格线
plt = plot(grid=true) # <-- 这里初始化时直接启用网格

# 运行模型比较
models = ["sP2D", "P2D"]

for model in models
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    
    # 绘图时保持网格设置
    plot!(plt, result["time [s]"], result["cell voltage [V]"],
          label=model, lw=2, xlabel="Time [s]", ylabel="Cell Voltage [V]", 
          grid=true) # <-- 这里确保每次叠加都保持网格
end

# 图形修饰
xlims!(0, total_time)
title!("Sawtooth Charging → Sawtooth Discharging")
savefig("Sawtooth_charge_sawtooth_discharge.pdf")

# 引用
JuBat.Citation(["ai2024b"])