using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# 定义锯齿波电流参数
amplitude = 3.0  # 锯齿波的幅值（安培）
period = 500     # 锯齿波的周期（秒）
total_time = 3000  # 总模拟时间（秒）

opt.dtType = "fixed"
opt.time = [0, total_time]  # 设置总模拟时间

# 定义锯齿波电流函数
opt.Current = (t) -> amplitude * (t / period - floor(t / period))

# 运行模型模拟
for model in ["sP2D", "P2D"]
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # 绘制结果
    plot!(result["time [s]"], result["cell voltage [V]"], label=model, xlabel="time [s]", ylabel="cell voltage [V]", lw=2)
end

# 保存图形
savefig("sawtooth_current_validation.pdf")

# 引用
JuBat.Citation(["ai2024b"])