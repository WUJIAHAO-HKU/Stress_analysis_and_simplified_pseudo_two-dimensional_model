using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("Enertech")
param_dim.cell.v_h = 4.3

# 定义高频脉冲电流参数
total_time = 1000      # Total simulation time in seconds
time_step = 1          # 时间步长

# 创建时间数组
time = collect(Float64, 0:time_step:total_time)

# 创建类似于心电图的电流波形
current = similar(time)
for i in eachindex(time)
    t = time[i]
    # 基础正弦波
    base_wave = sin(2 * π * t) * 0.5 + 0.5
    # 添加高频成分
    high_freq = sin(4 * π * t) * 0.3 - 0.1
    # 添加更高频成分
    higher_freq = sin(8 * π * t) * 0.2 - 0.1
    # 在波形中加入尖峰
    peak = 0.0
    for j in 1:10
        peak_position = j * 1.0
        peak += exp(-((t - peak_position) ^ 2) * 10) * 0.2
    end
    # 添加噪声
    noise = randn() * 0.1
    # 合成电流
    current[i] = base_wave + high_freq + higher_freq + peak + noise
end

# 创建电流插值函数
current_interp = LinearInterpolation(time, current, extrapolation_bc=Flat())

# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = time  # 使用与电流数据相同的时间步长

# 设置电流函数
opt.Current = t -> current_interp(t)

# 创建和运行P2D模型模拟
println("开始模拟高频脉冲电流工况... P2D模型")
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
println("P2D模型模拟完成")

# 测量P2D模型计算时间
p2d_time = @elapsed begin
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
end
println("P2D模型计算时间: ", p2d_time, " 秒")

# 提取P2D模型结果数据
time_result = result["time [s]"]
x_grid_n = case.mesh["negative electrode"].node[:, 1] # 负极网格
x_grid_p = case.mesh["positive electrode"].node[:, 1]  # 正极网格，单位：米 (m)

# 创建归一化的坐标系统
# 根据LG M50电池的实际厚度进行归一化
# 负极厚度: 85.2μm, 隔膜厚度: 12μm, 正极厚度: 75.6μm, 总厚度: 172.8μm
thickness_negative = 85.2e-6  # 负极厚度 (m)
thickness_separator = 12e-6   # 隔膜厚度 (m)  
thickness_positive = 75.6e-6  # 正极厚度 (m)
total_thickness = thickness_negative + thickness_separator + thickness_positive  # 总厚度

# 负极归一化到 [0, 85.2/172.8] ≈ [0, 0.4931]
negative_end = thickness_negative / total_thickness
x_grid_n_normalized = x_grid_n ./ maximum(x_grid_n) .* negative_end

# 正极归一化到 [(85.2+12)/172.8, 1.0] ≈ [0.5625, 1.0]
positive_start = (thickness_negative + thickness_separator) / total_thickness
positive_end = 1.0
x_grid_p_normalized = positive_start .+ (x_grid_p .- minimum(x_grid_p)) ./ (maximum(x_grid_p) - minimum(x_grid_p)) .* (positive_end - positive_start)

# 选择高频脉冲周期中的关键时间点进行空间分布分析
# 选择10个关键时间点，包括尖峰点和谷值点
analysis_times = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
time_indices = []

for t in analysis_times
    idx = argmin(abs.(time_result .- t))
    push!(time_indices, idx)
end

voltage_p2d = result["cell voltage [V]"]
concentration_p2d = result["negative particle surface lithium concentration [mol/m^3]"]
stress_p2d = result["negative particle surface tangential stress[Pa]"]
concentration_p_p2d = result["positive particle surface lithium concentration [mol/m^3]"]
stress_p_p2d = result["positive particle surface tangential stress[Pa]"]

# 无量纲化参数设置
# 浓度无量纲化 (c/c_max)
c_max_neg = param_dim.NE.cs_max  # 负极最大浓度
c_max_pos = param_dim.PE.cs_max  # 正极最大浓度

# 应力无量纲化 (σ/E)
E_neg = param_dim.NE.E  # 负极杨氏模量
E_pos = param_dim.PE.E  # 正极杨氏模量

# 设置sP2D模型
opt.model = "sP2D"  # 使用sP2D模型
println("开始模拟高频脉冲电流工况... sP2D模型")
case_sP2D = JuBat.SetCase(param_dim, opt)
result_sP2D = JuBat.Solve(case_sP2D)
println("sP2D模型模拟完成")

# 测量sP2D模型计算时间
sp2d_time = @elapsed begin
    case_sP2D = JuBat.SetCase(param_dim, opt)
    result_sP2D = JuBat.Solve(case_sP2D)
end

println("P2D模型计算时间: ", p2d_time, " 秒")
println("sP2D模型计算时间: ", sp2d_time, " 秒")

# 提取sP2D模型结果数据
voltage_sP2D = result_sP2D["cell voltage [V]"]
concentration_sP2D = result_sP2D["negative particle surface lithium concentration [mol/m^3]"]
stress_sP2D = result_sP2D["negative particle surface tangential stress[Pa]"]
concentration_p_sP2D = result_sP2D["positive particle surface lithium concentration [mol/m^3]"]
stress_p_sP2D = result_sP2D["positive particle surface tangential stress[Pa]"]

# 提取各时间点的空间分布数据
concentration_p2d_spatial = []
stress_p2d_spatial = []
concentration_p_p2d_spatial = []
stress_p_p2d_spatial = []

concentration_sP2D_spatial = []
stress_sP2D_spatial = []
concentration_p_sP2D_spatial = []
stress_p_sP2D_spatial = []

for (i, t_idx) in enumerate(time_indices)
    # P2D模型数据 - 应用无量纲化
    if ndims(concentration_p2d) > 1
        # 浓度无量纲化 (c/c_max)
        conc_normalized = concentration_p2d[:, t_idx] ./ c_max_neg
        push!(concentration_p2d_spatial, conc_normalized)
    end
    
    if ndims(stress_p2d) > 1
        # 应力无量纲化 (σ/E)
        stress_normalized = stress_p2d[:, t_idx] ./ E_neg
        push!(stress_p2d_spatial, stress_normalized)
    end
    
    if ndims(concentration_p_p2d) > 1
        # 正极浓度无量纲化 (c/c_max)
        conc_p_normalized = concentration_p_p2d[:, t_idx] ./ c_max_pos
        push!(concentration_p_p2d_spatial, conc_p_normalized)
    end
    
    if ndims(stress_p_p2d) > 1
        # 正极应力无量纲化 (σ/E)
        stress_p_normalized = stress_p_p2d[:, t_idx] ./ E_pos
        push!(stress_p_p2d_spatial, stress_p_normalized)
    end
    
    # sP2D模型数据 - 应用无量纲化
    if ndims(concentration_sP2D) > 1
        # 浓度无量纲化 (c/c_max)
        conc_normalized = concentration_sP2D[:, t_idx] ./ c_max_neg
        push!(concentration_sP2D_spatial, conc_normalized)
    end
    
    if ndims(stress_sP2D) > 1
        # 应力无量纲化 (σ/E)
        stress_normalized = stress_sP2D[:, t_idx] ./ E_neg
        push!(stress_sP2D_spatial, stress_normalized)
    end
    
    if ndims(concentration_p_sP2D) > 1
        # 正极浓度无量纲化 (c/c_max)
        conc_p_normalized = concentration_p_sP2D[:, t_idx] ./ c_max_pos
        push!(concentration_p_sP2D_spatial, conc_p_normalized)
    end
    
    if ndims(stress_p_sP2D) > 1
        # 正极应力无量纲化 (σ/E)
        stress_p_normalized = stress_p_sP2D[:, t_idx] ./ E_pos
        push!(stress_p_sP2D_spatial, stress_p_normalized)
    end
end

# 定义科学论文风格的配色方案
# 使用专业期刊常用的配色方案，确保色盲友好和打印友好
time_colors = [
    :steelblue, :darkorange, :forestgreen, :crimson, :mediumpurple,  # 前5个时间点
    :saddlebrown, :hotpink, :dimgray, :darkolivegreen, :darkcyan     # 后5个时间点
]

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time_result, voltage_p2d, 
    label="P2D Model", 
    xlabel="Time [s]", 
    ylabel="Voltage [V]", 
    title="Cell Voltage Comparison",
    lw=2.5, color=:steelblue,
    grid=true,
    gridalpha=0.3,
    framestyle=:box,
    size=(800, 400),dpi = 600)

plot!(time_result, voltage_sP2D, 
    label="sP2D Model", 
    linestyle=:dash, color=:crimson, lw=2.5)

# 添加分析时间点的垂直线
for (i, t) in enumerate(analysis_times)
    vline!([t], color=:gray, alpha=0.4, linestyle=:dot, linewidth=1, label="")
end

# 创建浓度空间分布图
p2 = plot(xlabel="Normalized Position", 
         ylabel="c/c_max (dimensionless)", 
         title="Electrode Concentration Distribution - Normalized",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的浓度分布
for i in 1:length(time_indices)
    line_style_p2d = i <= 5 ? :solid : :dash      # P2D模型线型
    line_style_sp2d = i <= 5 ? :dot : :dashdot    # sP2D模型线型
    line_width = i <= 5 ? 2.0 : 1.5
    
    # 正极浓度 - P2D (实线)
    plot!(x_grid_p_normalized, concentration_p_p2d_spatial[i], 
        label=(i==1 ? "Positive (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=0.8)
    
    # 正极浓度 - sP2D (虚线)
    plot!(x_grid_p_normalized, concentration_p_sP2D_spatial[i], 
        label=(i==1 ? "Positive (sP2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_sp2d, alpha=0.8)
    
    # 负极浓度 - P2D (实线)
    plot!(x_grid_n_normalized, concentration_p2d_spatial[i], 
        label=(i==1 ? "Negative (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=0.8)
    
    # 负极浓度 - sP2D (虚线)
    plot!(x_grid_n_normalized, concentration_sP2D_spatial[i], 
        label=(i==1 ? "Negative (sP2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_sp2d, alpha=0.8)
end

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域填充
separator_mid = (negative_end + positive_start) / 2
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="Separator")

# 区域标识
neg_center = negative_end / 2
pos_center = (positive_start + positive_end) / 2
annotate!(neg_center, maximum([maximum(conc) for conc in concentration_p2d_spatial])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center, maximum([maximum(conc) for conc in concentration_p_p2d_spatial])*0.95, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(conc) for conc in concentration_p2d_spatial])*0.85, 
         text("Separator", 8, :center, :italic, :gray))

# 创建应力空间分布图
p3 = plot(xlabel="Normalized Position", 
         ylabel="σ/E (dimensionless)", 
         title="Electrode Stress Distribution - Normalized",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的应力分布
for i in 1:length(time_indices)
    line_style_p2d = i <= 5 ? :solid : :dash      # P2D模型线型
    line_style_sp2d = i <= 5 ? :dot : :dashdot    # sP2D模型线型
    line_width = i <= 5 ? 2.0 : 1.5
    
    # 正极应力 - P2D (实线)
    plot!(x_grid_p_normalized, stress_p_p2d_spatial[i], 
        label=(i==1 ? "Positive (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=0.8)
    
    # 正极应力 - sP2D (虚线)
    plot!(x_grid_p_normalized, stress_p_sP2D_spatial[i], 
        label=(i==1 ? "Positive (sP2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_sp2d, alpha=0.8)
    
    # 负极应力 - P2D (实线)
    plot!(x_grid_n_normalized, stress_p2d_spatial[i], 
        label=(i==1 ? "Negative (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=0.8)
    
    # 负极应力 - sP2D (虚线)
    plot!(x_grid_n_normalized, stress_sP2D_spatial[i], 
        label=(i==1 ? "Negative (sP2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_sp2d, alpha=0.8)
end

# 添加分隔线和区域标识  
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识
annotate!(neg_center, maximum([maximum(stress) for stress in stress_p2d_spatial])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center, maximum([maximum(stress) for stress in stress_p_p2d_spatial])*0.95, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(stress) for stress in stress_p2d_spatial])*0.85, 
         text("Separator", 8, :center, :italic, :gray))

# 创建电流图
p4 = plot(time, current, 
    label="Current Profile", 
    xlabel="Time [s]", 
    ylabel="Current [A]", 
    title="High Frequency Pulse Current Profile",
    lw=2.5, color=:darkorange,
    grid=true,
    gridalpha=0.3,
    framestyle=:box,
    size=(800, 400),dpi = 600)

# 添加分析时间点的垂直线
for (i, t) in enumerate(analysis_times)
    vline!([t], color=:gray, alpha=0.4, linestyle=:dot, linewidth=1, label="")
end

# 组合图表（原始版本）
plot_combined = plot(p1, p2, p3, p4, 
                    layout=(4, 1), 
                    size=(1000, 1200))

# 保存图表
savefig(p1, "pulse_voltage_spatial_comparison.pdf")
savefig(p2, "pulse_concentration_spatial_comparison.png")
savefig(p3, "pulse_stress_spatial_comparison.png")
savefig(p4, "pulse_current_profile_spatial.pdf")
savefig(plot_combined, "pulse_spatial_analysis_comparison_super.pdf")

# 计算空间分布误差
concentration_spatial_errors_n = []  # 负极浓度误差
stress_spatial_errors_n = []         # 负极应力误差
concentration_spatial_errors_p = []  # 正极浓度误差
stress_spatial_errors_p = []         # 正极应力误差

for i in 1:length(time_indices)
    # 负极浓度误差
    conc_error_n = abs.(concentration_sP2D_spatial[i] .- concentration_p2d_spatial[i])
    push!(concentration_spatial_errors_n, conc_error_n)
    
    # 负极应力误差
    stress_error_n = abs.(stress_sP2D_spatial[i] .- stress_p2d_spatial[i])
    push!(stress_spatial_errors_n, stress_error_n)
    
    # 正极浓度误差
    conc_error_p = abs.(concentration_p_sP2D_spatial[i] .- concentration_p_p2d_spatial[i])
    push!(concentration_spatial_errors_p, conc_error_p)
    
    # 正极应力误差
    stress_error_p = abs.(stress_p_sP2D_spatial[i] .- stress_p_p2d_spatial[i])
    push!(stress_spatial_errors_p, stress_error_p)
end

# 创建空间误差分布图 - 浓度误差
p5 = plot(xlabel="Normalized Position", 
         ylabel="Concentration Error (c/c_max)", 
         title="Spatial Concentration Error Distribution - Normalized",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用不同颜色和线型
for i in 1:length(time_indices)
    line_style = i <= 5 ? :solid : :dash  # 前5个用实线，后5个用虚线
    line_width = i <= 5 ? 2.0 : 1.5      # 前5个线宽稍粗
    
    plot!(x_grid_n_normalized, concentration_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=line_width, 
        color=time_colors[i],
        linestyle=line_style,
        alpha=0.8)
end

# 绘制正极误差 - 使用相同颜色但不同线型
for i in 1:length(time_indices)
    line_style = i <= 5 ? :dot : :dashdot  # 前5个用点线，后5个用点划线
    line_width = i <= 5 ? 2.0 : 1.5
    
    plot!(x_grid_p_normalized, concentration_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=line_width, 
        color=time_colors[i],
        linestyle=line_style,
        alpha=0.8)
end

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识
annotate!(neg_center, maximum([maximum(err) for err in concentration_spatial_errors_n])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center, maximum([maximum(err) for err in concentration_spatial_errors_p])*0.95, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in concentration_spatial_errors_n])*0.85, 
         text("Separator", 8, :center, :italic, :gray))

# 创建空间误差分布图 - 应力误差
p6 = plot(xlabel="Normalized Position", 
         ylabel="Stress Error (σ/E)", 
         title="Spatial Stress Error Distribution - Normalized",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用不同颜色和线型
for i in 1:length(time_indices)
    line_style = i <= 5 ? :solid : :dash  # 前5个用实线，后5个用虚线
    line_width = i <= 5 ? 2.0 : 1.5      # 前5个线宽稍粗
    
    plot!(x_grid_n_normalized, stress_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=line_width, 
        color=time_colors[i],
        linestyle=line_style,
        alpha=0.8)
end

# 绘制正极误差 - 使用相同颜色但不同线型
for i in 1:length(time_indices)
    line_style = i <= 5 ? :dot : :dashdot  # 前5个用点线，后5个用点划线
    line_width = i <= 5 ? 2.0 : 1.5
    
    plot!(x_grid_p_normalized, stress_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=line_width, 
        color=time_colors[i],
        linestyle=line_style,
        alpha=0.8)
end

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识
annotate!(neg_center, maximum([maximum(err) for err in stress_spatial_errors_n])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center, maximum([maximum(err) for err in stress_spatial_errors_p])*0.95, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in stress_spatial_errors_n])*0.85, 
         text("Separator", 8, :center, :italic, :gray))

savefig(p5, "pulse_concentration_spatial_error.png")
savefig(p6, "pulse_stress_spatial_error.png")

# 创建合并的浓度分布图（双y轴：原始数据+误差）
p7 = plot(xlabel="Normalized Position", 
         ylabel="Concentration (c/c_max)", 
         title="Electrode Concentration: Distribution & Error Analysis (Dual Y-axis)",
         legend=:outertopright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(1200, 700),dpi = 600)

# 左y轴：绘制原始浓度分布
for i in 1:length(time_indices)
    line_style_p2d = i <= 5 ? :solid : :solid      # P2D模型线型
    line_style_sp2d = i <= 5 ? :dot : :dashdot    # sP2D模型线型
    line_width = i <= 5 ? 2.5 : 2.0
    alpha_val = 0.9
    
    # 正极浓度 - P2D (较粗的实线)
    plot!(x_grid_p_normalized, concentration_p_p2d_spatial[i], 
        label=(i==1 ? "Positive (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=alpha_val)
    
    # 正极浓度 - sP2D (点线)
    plot!(x_grid_p_normalized, concentration_p_sP2D_spatial[i], 
        label=(i==1 ? "Positive (sP2D)" : ""), 
        lw=line_width-0.5, color=time_colors[i], linestyle=line_style_sp2d, alpha=alpha_val)
    
    # 负极浓度 - P2D (较粗的实线)
    plot!(x_grid_n_normalized, concentration_p2d_spatial[i], 
        label=(i==1 ? "Negative (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=alpha_val)
    
    # 负极浓度 - sP2D (点线)
    plot!(x_grid_n_normalized, concentration_sP2D_spatial[i], 
        label=(i==1 ? "Negative (sP2D)" : ""), 
        lw=line_width-0.5, color=time_colors[i], linestyle=line_style_sp2d, alpha=alpha_val)
end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="Separator")

# 区域标识
max_conc_val = maximum([maximum(conc) for conc in concentration_p2d_spatial])
annotate!(neg_center, max_conc_val*0.95, 
         text("Negative Electrode", 10, :center, :bold))
annotate!(pos_center, max_conc_val*0.95, 
         text("Positive Electrode", 10, :center, :bold))
annotate!(separator_mid, max_conc_val*0.85, 
         text("Separator", 9, :center, :italic, :gray))

# 创建右y轴用于显示误差
p7_twin = twinx(p7)
plot!(p7_twin, xlabel="", ylabel="Concentration Error (c/c_max)", 
      legend=:outerbottomright, framestyle=:box)

# 右y轴：绘制误差分布
for i in 1:length(time_indices)
    error_line_width = i <= 5 ? 2.0 : 1.5
    error_alpha = 0.8
    
    # 负极误差
    plot!(p7_twin, x_grid_n_normalized, concentration_spatial_errors_n[i], 
        label=(i==1 ? "Error (Negative)" : ""), 
        lw=error_line_width, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=error_alpha)
    
    # 正极误差
    plot!(p7_twin, x_grid_p_normalized, concentration_spatial_errors_p[i], 
        label=(i==1 ? "Error (Positive)" : ""), 
        lw=error_line_width, 
        color=time_colors[i],
        linestyle=:dashdot,
        alpha=error_alpha)
end

# 创建合并的应力分布图（双y轴：原始数据+误差）
p8 = plot(xlabel="Normalized Position", 
         ylabel="Stress (σ/E)", 
         title="Electrode Stress: Distribution & Error Analysis (Dual Y-axis)",
         legend=:outertopright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(1200, 700),dpi = 600)

# 左y轴：绘制原始应力分布
for i in 1:length(time_indices)
    line_style_p2d = i <= 5 ? :solid : :solid      # P2D模型线型
    line_style_sp2d = i <= 5 ? :dot : :dashdot    # sP2D模型线型
    line_width = i <= 5 ? 2.5 : 2.0
    alpha_val = 0.9
    
    # 正极应力 - P2D (较粗的实线)
    plot!(x_grid_p_normalized, stress_p_p2d_spatial[i], 
        label=(i==1 ? "Positive (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=alpha_val)
    
    # 正极应力 - sP2D (点线)
    plot!(x_grid_p_normalized, stress_p_sP2D_spatial[i], 
        label=(i==1 ? "Positive (sP2D)" : ""), 
        lw=line_width-0.5, color=time_colors[i], linestyle=line_style_sp2d, alpha=alpha_val)
    
    # 负极应力 - P2D (较粗的实线)
    plot!(x_grid_n_normalized, stress_p2d_spatial[i], 
        label=(i==1 ? "Negative (P2D)" : ""), 
        lw=line_width, color=time_colors[i], linestyle=line_style_p2d, alpha=alpha_val)
    
    # 负极应力 - sP2D (点线)
    plot!(x_grid_n_normalized, stress_sP2D_spatial[i], 
        label=(i==1 ? "Negative (sP2D)" : ""), 
        lw=line_width-0.5, color=time_colors[i], linestyle=line_style_sp2d, alpha=alpha_val)
end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="Separator")

# 区域标识
max_stress_val = maximum([maximum(stress) for stress in stress_p2d_spatial])
annotate!(neg_center, max_stress_val*0.95, 
         text("Negative Electrode", 10, :center, :bold))
annotate!(pos_center, max_stress_val*0.95, 
         text("Positive Electrode", 10, :center, :bold))
annotate!(separator_mid, max_stress_val*0.85, 
         text("Separator", 9, :center, :italic, :gray))

# 创建右y轴用于显示误差
p8_twin = twinx(p8)
plot!(p8_twin, xlabel="", ylabel="Stress Error (σ/E)", 
      legend=:outerbottomright, framestyle=:box)

# 右y轴：绘制误差分布
for i in 1:length(time_indices)
    error_line_width = i <= 5 ? 2.0 : 1.5
    error_alpha = 0.8
    
    # 负极应力误差
    plot!(p8_twin, x_grid_n_normalized, stress_spatial_errors_n[i], 
        label=(i==1 ? "Error (Negative)" : ""), 
        lw=error_line_width, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=error_alpha)
    
    # 正极应力误差
    plot!(p8_twin, x_grid_p_normalized, stress_spatial_errors_p[i], 
        label=(i==1 ? "Error (Positive)" : ""), 
        lw=error_line_width, 
        color=time_colors[i],
        linestyle=:dashdot,
        alpha=error_alpha)
end

# 保存双y轴合并图
savefig(p7, "pulse_concentration_spatial_dual_axis.png")
savefig(p8, "pulse_stress_spatial_dual_axis.png")

# 创建完整的双y轴合并布局
plot_combined_dual_axis = plot(p1, p7, p8, p4, 
                              layout=(4, 1), 
                              size=(1200, 1600))

savefig(plot_combined_dual_axis, "pulse_spatial_analysis_dual_axis.pdf")

# 计算和打印统计信息
max_spatial_conc_errors_n = [maximum(err) for err in concentration_spatial_errors_n]
max_spatial_stress_errors_n = [maximum(err) for err in stress_spatial_errors_n]
mean_spatial_conc_errors_n = [mean(err) for err in concentration_spatial_errors_n]
mean_spatial_stress_errors_n = [mean(err) for err in stress_spatial_errors_n]

max_spatial_conc_errors_p = [maximum(err) for err in concentration_spatial_errors_p]
max_spatial_stress_errors_p = [maximum(err) for err in stress_spatial_errors_p]
mean_spatial_conc_errors_p = [mean(err) for err in concentration_spatial_errors_p]
mean_spatial_stress_errors_p = [mean(err) for err in stress_spatial_errors_p]

# 合并正负极误差数据用于总体统计
all_conc_errors = vcat(concentration_spatial_errors_n, concentration_spatial_errors_p)
all_stress_errors = vcat(stress_spatial_errors_n, stress_spatial_errors_p)

println("=== 高频脉冲电流工况空间分布分析结果 ===")
println("分析时间点: ", analysis_times, " 秒")
println("负极最大空间浓度误差: ", max_spatial_conc_errors_n, " (c/c_max)")
println("负极平均空间浓度误差: ", mean_spatial_conc_errors_n, " (c/c_max)")
println("负极最大空间应力误差: ", max_spatial_stress_errors_n, " (σ/E)")
println("负极平均空间应力误差: ", mean_spatial_stress_errors_n, " (σ/E)")
println("正极最大空间浓度误差: ", max_spatial_conc_errors_p, " (c/c_max)")
println("正极平均空间浓度误差: ", mean_spatial_conc_errors_p, " (c/c_max)")
println("正极最大空间应力误差: ", max_spatial_stress_errors_p, " (σ/E)")
println("正极平均空间应力误差: ", mean_spatial_stress_errors_p, " (σ/E)")

# 总体最大误差
overall_max_conc_error = maximum([maximum(err) for err in all_conc_errors])
overall_max_stress_error = maximum([maximum(err) for err in all_stress_errors])
overall_mean_conc_error = mean([mean(err) for err in all_conc_errors])
overall_mean_stress_error = mean([mean(err) for err in all_stress_errors])

println("\n=== 总体统计 ===")
println("总体最大浓度误差: $overall_max_conc_error (c/c_max)")
println("总体平均浓度误差: $overall_mean_conc_error (c/c_max)")
println("总体最大应力误差: $overall_max_stress_error (σ/E)")
println("总体平均应力误差: $overall_mean_stress_error (σ/E)")

# 分析高频特性对空间分布的影响
println("\n=== 高频特性分析 ===")
high_freq_indices = [1, 3, 5, 7, 9]  # 选择奇数时间点分析高频特性
low_freq_indices = [2, 4, 6, 8, 10]  # 选择偶数时间点作为对比

high_freq_conc_errors_n = [concentration_spatial_errors_n[i] for i in high_freq_indices]
low_freq_conc_errors_n = [concentration_spatial_errors_n[i] for i in low_freq_indices]
high_freq_conc_errors_p = [concentration_spatial_errors_p[i] for i in high_freq_indices]
low_freq_conc_errors_p = [concentration_spatial_errors_p[i] for i in low_freq_indices]

high_freq_max_error_n = maximum([maximum(err) for err in high_freq_conc_errors_n])
low_freq_max_error_n = maximum([maximum(err) for err in low_freq_conc_errors_n])
high_freq_max_error_p = maximum([maximum(err) for err in high_freq_conc_errors_p])
low_freq_max_error_p = maximum([maximum(err) for err in low_freq_conc_errors_p])

println("负极高频段最大浓度误差: $high_freq_max_error_n (c/c_max)")
println("负极低频段最大浓度误差: $low_freq_max_error_n (c/c_max)")
println("负极高频影响系数: $(high_freq_max_error_n / low_freq_max_error_n)")
println("正极高频段最大浓度误差: $high_freq_max_error_p (c/c_max)")
println("正极低频段最大浓度误差: $low_freq_max_error_p (c/c_max)")
println("正极高频影响系数: $(high_freq_max_error_p / low_freq_max_error_p)")

# 保存空间分析数据到CSV
# 为负极创建DataFrame
negative_df = DataFrame()
negative_df[!, "Position"] = x_grid_n
negative_df[!, "Normalized_Position"] = x_grid_n_normalized

for i in 1:length(time_indices)
    t = analysis_times[i]
    negative_df[!, "Conc_P2D_t$(t)s_(c_cmax)"] = concentration_p2d_spatial[i]
    negative_df[!, "Conc_sP2D_t$(t)s_(c_cmax)"] = concentration_sP2D_spatial[i]
    negative_df[!, "Conc_Error_t$(t)s_(c_cmax)"] = concentration_spatial_errors_n[i]
    negative_df[!, "Stress_P2D_t$(t)s_(sigma_E)"] = stress_p2d_spatial[i]
    negative_df[!, "Stress_sP2D_t$(t)s_(sigma_E)"] = stress_sP2D_spatial[i]
    negative_df[!, "Stress_Error_t$(t)s_(sigma_E)"] = stress_spatial_errors_n[i]
    # 为每个位置添加对应的电流值（重复标量值）
    negative_df[!, "Current_t$(t)s"] = fill(current_interp(t), length(x_grid_n))
end

# 为正极创建DataFrame
positive_df = DataFrame()
positive_df[!, "Position"] = x_grid_p
positive_df[!, "Normalized_Position"] = x_grid_p_normalized

for i in 1:length(time_indices)
    t = analysis_times[i]
    positive_df[!, "Conc_P2D_t$(t)s_(c_cmax)"] = concentration_p_p2d_spatial[i]
    positive_df[!, "Conc_sP2D_t$(t)s_(c_cmax)"] = concentration_p_sP2D_spatial[i]
    positive_df[!, "Conc_Error_t$(t)s_(c_cmax)"] = concentration_spatial_errors_p[i]
    positive_df[!, "Stress_P2D_t$(t)s_(sigma_E)"] = stress_p_p2d_spatial[i]
    positive_df[!, "Stress_sP2D_t$(t)s_(sigma_E)"] = stress_p_sP2D_spatial[i]
    positive_df[!, "Stress_Error_t$(t)s_(sigma_E)"] = stress_spatial_errors_p[i]
    # 为每个位置添加对应的电流值（重复标量值）
    positive_df[!, "Current_t$(t)s"] = fill(current_interp(t), length(x_grid_p))
end

# 保存两个DataFrame
CSV.write("pulse_spatial_analysis_negative_electrode.csv", negative_df)
CSV.write("pulse_spatial_analysis_positive_electrode.csv", positive_df)

println("高频脉冲电流空间分布分析完成，结果已保存")
println("\n=== 图表文件说明 ===")
println("生成的图表文件：")
println("  原始分布图:")
println("    - pulse_concentration_spatial_comparison.png (浓度分布)")
println("    - pulse_stress_spatial_comparison.png (应力分布)")
println("  误差分布图:")
println("    - pulse_concentration_spatial_error.png (浓度误差)")
println("    - pulse_stress_spatial_error.png (应力误差)")
println("  双y轴分析图 (解决量级差异问题):")
println("    - pulse_concentration_spatial_dual_axis.png (浓度分布+误差，双y轴)")
println("    - pulse_stress_spatial_dual_axis.png (应力分布+误差，双y轴)")
println("  综合布局:")
println("    - pulse_spatial_analysis_comparison_super.pdf (4张分离图)")
println("    - pulse_spatial_analysis_dual_axis.pdf (双y轴合并分析图)")
println()
println("双y轴图表特色:")
println("  1. 左y轴：显示原始数据（浓度/应力分布）")
println("  2. 右y轴：显示误差数据（解决量级差异问题）")
println("  3. 实线和点线显示P2D与sP2D模型结果")
println("  4. 虚线和点划线显示误差分布")
println("  5. 准确反映电池物理结构（负极-隔膜-正极）")
println("  6. 解决原始数据与误差量级差异导致的显示问题")

# 无量纲化说明
println("\n=== 无量纲化处理说明 ===")
println("空间坐标归一化:")
println("  负极范围: [0, $(round(negative_end, digits=4))] (对应实际厚度: 85.2μm)")
println("  隔膜范围: [$(round(negative_end, digits=4)), $(round(positive_start, digits=4))] (对应实际厚度: 12μm)")
println("  正极范围: [$(round(positive_start, digits=4)), 1.0] (对应实际厚度: 75.6μm)")
println("  总厚度: $(thickness_negative*1e6 + thickness_separator*1e6 + thickness_positive*1e6)μm")
println()
println("浓度无量纲化: c/c_max (荷电状态SOC，0-1范围)")
println("  负极最大浓度: $c_max_neg mol/m³")
println("  正极最大浓度: $c_max_pos mol/m³")
println("应力无量纲化: σ/E (相对应力水平)")
println("  负极杨氏模量: $E_neg Pa")
println("  正极杨氏模量: $E_pos Pa")
println("无量纲化优势:")
println("  1. 提高数值稳定性")
println("  2. 便于不同工况比较")
println("  3. 物理意义更清晰")
println("  4. 符合工程实践")
println("  5. 准确反映电池结构")

# 引用信息
JuBat.Citation()
