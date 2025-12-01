using Plots, CSV, DataFrames, Statistics, Interpolations
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("Enertech")
opt = JuBat.Option()

rates = [0.1, 0.5, 1, 2]
ts = [3700, 3700, 3600, 3600]

# 定义颜色方案 - 蓝色系列，倍率越小颜色越淡
colors = [:lightsteelblue, :royalblue, :blue, :darkblue]  # sP2D模型
p2d_colors = [:lightsteelblue, :royalblue, :blue, :darkblue]  # P2D模型

# 初始化图形
p = plot(legendfontsize=16)

# 创建误差汇总数据结构
error_summary = Dict(
    "C_rate" => Float64[],
    "RMSE_voltage" => Float64[],
    "MAE_voltage" => Float64[],
    "max_abs_error" => Float64[],
    "error_at_10pct" => Float64[],
    "error_at_50pct" => Float64[],
    "error_at_90pct" => Float64[],
    "capacity_diff_pct" => Float64[],
    "computation_time_sp2d" => Float64[],
    "computation_time_p2d" => Float64[]
)

# 辅助函数：计算累积容量
function calculate_capacity(time, current)
    capacity = zeros(length(time))
    for i in 2:length(time)
        dt = time[i] - time[i-1]
        capacity[i] = capacity[i-1] + abs(current[i]) * dt / 3600.0  # 转换为Ah
    end
    return capacity
end

# 计算插值后的电压误差
function calculate_voltage_errors(time_sp2d, voltage_sp2d, time_p2d, voltage_p2d, capacity_sp2d, capacity_p2d, max_capacity)
    # 创建基于容量的插值函数
    sp2d_voltage_interp = linear_interpolation(capacity_sp2d, voltage_sp2d, extrapolation_bc=Line())
    p2d_voltage_interp = linear_interpolation(capacity_p2d, voltage_p2d, extrapolation_bc=Line())
    
    # 创建均匀分布的容量点进行比较
    capacity_points = range(0.0, max_capacity, length=1000)
    
    # 计算这些点上的电压值
    sp2d_voltages = sp2d_voltage_interp.(capacity_points)
    p2d_voltages = p2d_voltage_interp.(capacity_points)
    
    # 计算误差
    abs_errors = abs.(sp2d_voltages .- p2d_voltages)
    rel_errors = abs_errors ./ p2d_voltages * 100  # 相对误差百分比
    
    # 计算各种误差指标
    rmse = sqrt(mean(abs_errors.^2))
    mae = mean(abs_errors)
    max_error = maximum(abs_errors)
    
    # 在特定容量点的误差 (10%, 50%, 90%)
    error_10pct = abs_errors[Int(round(length(abs_errors) * 0.1))]
    error_50pct = abs_errors[Int(round(length(abs_errors) * 0.5))]
    error_90pct = abs_errors[Int(round(length(abs_errors) * 0.9))]
    
    return Dict(
        "RMSE" => rmse,
        "MAE" => mae,
        "max_abs_error" => max_error,
        "error_at_10pct" => error_10pct,
        "error_at_50pct" => error_50pct,
        "error_at_90pct" => error_90pct,
        "abs_errors" => abs_errors,
        "rel_errors" => rel_errors,
        "capacity_points" => capacity_points
    )
end

println("="^80)
println("| C-rate | RMSE (V) | MAE (V) | Max Error (V) | Error at 10% | Error at 50% | Error at 90% | Capacity Diff (%) | SP2D Time (s) | P2D Time (s) | Speedup |")
println("|" * "-"^8 * "|" * "-"^10 * "|" * "-"^9 * "|" * "-"^14 * "|" * "-"^13 * "|" * "-"^13 * "|" * "-"^13 * "|" * "-"^18 * "|" * "-"^14 * "|" * "-"^14 * "|" * "-"^9 * "|")

for j = 1:length(rates)
    Crate = rates[j]
    
    # 针对不同倍率设置不同的时间步长
    if j == 3
        # 1C倍率使用更短的时间步长
        opt.dt = [0.5, 2.0]
    elseif j == 4
        # 2C倍率使用更短的时间步长
        opt.dt = [0.2, 1.0]
    else
        # 0.1C和0.5C使用原来的时间步长
        opt.dt = [1, 10]/Crate
    end
    
    opt.dtType = "fixed"
    i = 5*Crate
    opt.Current = x-> i
    opt.time = [0, ts[j]/Crate]
    
    # 计时并运行SP2D模型
    opt.model = "sP2D"
    case = JuBat.SetCase(param_dim, opt)
    sp2d_time = @elapsed result = JuBat.Solve(case)
    
    # 计时并运行P2D模型
    opt.model = "P2D"
    case2 = JuBat.SetCase(param_dim, opt)
    p2d_time = @elapsed result2 = JuBat.Solve(case2)

    # 计算容量
    capacity_sp2d = calculate_capacity(result["time [s]"], result["cell current [A]"])
    capacity_p2d = calculate_capacity(result2["time [s]"], result2["cell current [A]"])
    
    # 计算最大容量和容量差异
    max_capacity_sp2d = maximum(capacity_sp2d)
    max_capacity_p2d = maximum(capacity_p2d)
    max_capacity = min(max_capacity_sp2d, max_capacity_p2d)
    capacity_diff_pct = abs(max_capacity_sp2d - max_capacity_p2d) / max_capacity_p2d * 100
    
    # 计算电压误差
    errors = calculate_voltage_errors(
        result["time [s]"], result["cell voltage [V]"],
        result2["time [s]"], result2["cell voltage [V]"],
        capacity_sp2d, capacity_p2d, max_capacity
    )
    
    # 存储误差数据
    push!(error_summary["C_rate"], Crate)
    push!(error_summary["RMSE_voltage"], errors["RMSE"])
    push!(error_summary["MAE_voltage"], errors["MAE"])
    push!(error_summary["max_abs_error"], errors["max_abs_error"])
    push!(error_summary["error_at_10pct"], errors["error_at_10pct"])
    push!(error_summary["error_at_50pct"], errors["error_at_50pct"])
    push!(error_summary["error_at_90pct"], errors["error_at_90pct"])
    push!(error_summary["capacity_diff_pct"], capacity_diff_pct)
    push!(error_summary["computation_time_sp2d"], sp2d_time)
    push!(error_summary["computation_time_p2d"], p2d_time)
    
    # 打印当前倍率的详细误差统计
    speedup = p2d_time / sp2d_time
    println("| $(Crate)C | $(round(errors["RMSE"], digits=5)) | $(round(errors["MAE"], digits=5)) | $(round(errors["max_abs_error"], digits=5)) | $(round(errors["error_at_10pct"], digits=5)) | $(round(errors["error_at_50pct"], digits=5)) | $(round(errors["error_at_90pct"], digits=5)) | $(round(capacity_diff_pct, digits=3)) | $(round(sp2d_time, digits=3)) | $(round(p2d_time, digits=3)) | $(round(speedup, digits=2))x |")
    
    # ploting results with improved styling
    plot!(p, capacity_p2d, result2["cell voltage [V]"], 
          label="P2D-$(Crate)C", 
          xlabel="Output capacity [Ah]", 
          ylabel="Cell voltage [V]", 
          linecolor=p2d_colors[j],
          linestyle=:solid,
          linewidth=2)
    plot!(p, capacity_sp2d, result["cell voltage [V]"], 
          label="sP2D-$(Crate)C", 
          linecolor=colors[j],
          linestyle=:dash, 
          linewidth=2.5)
end

println("="^80)
println("误差统计汇总:")
println("平均RMSE电压误差: $(mean(error_summary["RMSE_voltage"]))")
println("平均MAE电压误差: $(mean(error_summary["MAE_voltage"]))")
println("最大误差: $(maximum(error_summary["max_abs_error"]))")
println("平均计算加速比: $(mean(error_summary["computation_time_p2d"] ./ error_summary["computation_time_sp2d"]))x")
println()

# 分析误差与倍率的关系
println("误差与倍率关系分析:")
for metric in ["RMSE_voltage", "MAE_voltage", "max_abs_error", "computation_time_sp2d", "computation_time_p2d"]
    # 计算相关系数
    correlation = cor(error_summary["C_rate"], error_summary[metric])
    println("$(metric)与倍率的相关系数: $(round(correlation, digits=4))")
end

# 额外的图形调整
plot!(p,
    legend=:topright,          # 图例位置
    framestyle=:box,           # 边框样式
    grid=false,                # 不显示网格
    size=(800, 600),           # 图形大小
    dpi=300,                   # 分辨率
    margin=5Plots.mm,          # 边距
    legendfonthalign=:left     # 图例对齐方式
)

savefig(p, "sP2D_validation_Enertech.png")

JuBat.Citation(["ai2024b"])

println("\n完整的误差数据已输出。图表已保存为 'sP2D_validation_Enertech.png'")