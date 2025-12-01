using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 交流电流（AC）参数
amplitude = 2.0  # 电流幅度（安培）
frequency = 0.1  # 频率（Hz）
offset = 0.0  # 偏移量（安培），通常设为0

total_time = 5000  # 总模拟时间（秒）

opt.dtType = "fixed"
opt.time = [0, total_time]  # 设置总模拟时间

# 设置为正弦波电流
opt.Current = (t) -> begin
    return amplitude * sin(2 * π * frequency * t) + offset
end

# 运行模型模拟
for model in ["sP2D", "P2D"]
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # 绘制结果
    plot!(result["time [s]"], result["cell voltage [V]"], label=model, xlabel="time [s]", ylabel="cell voltage [V]", lw=2)
end

# 保存图形
savefig("ac_current_validation.pdf")

# 引用
JuBat.Citation(["ai2024b"])