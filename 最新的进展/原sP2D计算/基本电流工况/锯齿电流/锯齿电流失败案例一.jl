using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 锯齿波电流参数
max_current = 2.0  # 最大电流（安培）
min_current = -2.0  # 最小电流（安培）
period = 1000  # 一个周期的时间（秒）
total_time = 5000  # 总模拟时间（秒）

opt.dtType = "fixed"
opt.time = [0, total_time]  # 设置总模拟时间

# 设置电流为锯齿波
opt.Current = (t) -> begin
    return min_current + (max_current - min_current) * (mod(t, period) / period)
end

# 运行模型模拟
for model in ["sP2D", "P2D"]
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # 绘制结果
    plot!(result["time [s]"], result["cell voltage [V]"], label=model, xlabel="时间 [s]", ylabel="电池电压 [V]", lw=2)
end

# 保存图形
savefig("sawtooth_current_validation.pdf")

# 引用
JuBat.Citation(["ai2024b"])
