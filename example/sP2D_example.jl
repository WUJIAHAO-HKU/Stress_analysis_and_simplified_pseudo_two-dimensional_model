using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
rates = [ 0.5, 1, 2,5]
ts = [3700, 3700, 3600, 3600]

# 定义颜色方案 - 蓝色系列，倍率越小颜色越淡
colors = [:lightsteelblue, :royalblue, :blue, :darkblue]  # sP2D模型
p2d_colors = [:lightsteelblue, :royalblue, :blue, :darkblue]  # P2D模型

# 辅助函数：计算累积容量
function calculate_capacity(time, current)
    capacity = zeros(length(time))
    for i in 2:length(time)
        dt = time[i] - time[i-1]
        capacity[i] = capacity[i-1] + abs(current[i]) * dt / 3600.0  # 转换为Ah
    end
    return capacity
end

for j =1:4
    Crate = rates[j]
    opt.dt = [1, 20]/Crate
    opt.dtType = "fixed"
    i= 5*Crate
    opt.Current = x-> i
    opt.time = [ 0 ts[j]/Crate]
    opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    opt.model = "P2D"
    case2 = JuBat.SetCase(param_dim, opt)
    result2 = JuBat.Solve(case2)

    # 计算容量
    capacity_sp2d = calculate_capacity(result["time [s]"], result["cell current [A]"])
    capacity_p2d = calculate_capacity(result2["time [s]"], result2["cell current [A]"])

    # ploting results with improved styling
    plot!(capacity_p2d, result2["cell voltage [V]"], 
          label="P2D-$(Crate)C", 
          xlabel="Output capacity [Ah]", 
          ylabel="Cell voltage [V]", 
          linecolor=p2d_colors[j],
          linestyle=:solid,
          linewidth=2)
    plot!(capacity_sp2d, result["cell voltage [V]"], 
          label="sP2D-$(Crate)C", 
          linecolor=colors[j],
          linestyle=:dash, 
          linewidth=2.5)
end

savefig("sP2D_validation_LGM50.pdf")
JuBat.Citation(["ai2024b"])