using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

rates = [ 0.1,0.5, 1, 2]
ts = [3700, 3700, 3600, 3600]

# 定义颜色方案 - 蓝色系列，倍率越小颜色越淡
colors = [:lightsteelblue, :royalblue, :blue, :darkblue]  # sP2D模型
p2d_colors = [:lightsteelblue, :royalblue, :blue, :darkblue]  # P2D模型

# 初始化图形
p = plot()

# 辅助函数：计算SP2D和P2D数据的平均值（用于比较）
function calculate_average_data(time_sp2d, data_sp2d, time_p2d, data_p2d)
    # 使用线性插值将P2D的数据值插值到SP2D的时间点上
    data_p2d_interp = zeros(length(time_sp2d))
    
    for i in 1:length(time_sp2d)
        t = time_sp2d[i]
        if t <= time_p2d[1]
            data_p2d_interp[i] = data_p2d[1]
        elseif t >= time_p2d[end]
            data_p2d_interp[i] = data_p2d[end]
        else
            # 线性插值
            for j in 1:(length(time_p2d)-1)
                if time_p2d[j] <= t <= time_p2d[j+1]
                    α = (t - time_p2d[j]) / (time_p2d[j+1] - time_p2d[j])
                    data_p2d_interp[i] = data_p2d[j] * (1 - α) + data_p2d[j+1] * α
                    break
                end
            end
        end
    end
    
    # 计算平均值
    data_avg = (data_sp2d .+ data_p2d_interp) ./ 2
    return data_avg
end

for j =1:4
    Crate = rates[j]
    
    # 针对不同倍率设置不同的时间步长
    if j == 3
        # 2C倍率使用更短的时间步长
        opt.dt = [0.5, 2.0]
    elseif j == 4
        # 5C倍率使用更短的时间步长
        opt.dt = [0.2, 1.0]
    else
        # 0.5C和1C使用原来的时间步长
        opt.dt = [1, 10]/Crate
    end
    
    opt.dtType = "fixed"
    i= 5*Crate
    opt.Current = x-> i
    opt.time = [0, ts[j]/Crate]
    opt.model = "sP2D" # choose model, other options are "SPM" or "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    opt.model = "P2D"
    case2 = JuBat.SetCase(param_dim, opt)
    result2 = JuBat.Solve(case2)

    # 获取负极第一个位置的颗粒表面锂离子浓度
    # 注意：负极表面浓度是一个矩阵，第一个位置对应的是[1,:]
    neg_surf_conc_sp2d = result["negative particle surface lithium concentration [mol/m^3]"][1,:]
    neg_surf_conc_p2d = result2["negative particle surface lithium concentration [mol/m^3]"][1,:]
    
    # 获取时间数据
    time_sp2d = result["time [s]"]
    time_p2d = result2["time [s]"]

    # 计算SP2D和P2D表面浓度的平均值（用于比较）
    neg_surf_conc_avg = calculate_average_data(
        time_sp2d, neg_surf_conc_sp2d,
        time_p2d, neg_surf_conc_p2d
    )

    # 绘制结果，改为时间和浓度的关系
    plot!(p, time_p2d, neg_surf_conc_p2d, 
          label="P2D-$(Crate)C", 
          xlabel="Time [s]", 
          ylabel="Concentration [mol/m^3]", 
          linecolor=p2d_colors[j],
          linestyle=:solid,
          linewidth=2)
    plot!(p, time_sp2d, neg_surf_conc_avg, 
          label="sP2D-$(Crate)C", 
          linecolor=colors[j],
          linestyle=:dash, 
          linewidth=2.5)
end
    savefig(p, "负极表面锂离子浓度比较_LG_M50.pdf")

JuBat.Citation(["ai2024b"])