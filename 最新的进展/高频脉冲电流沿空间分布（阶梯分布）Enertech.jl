using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("Enertech")
param_dim.cell.v_h = 4.3

# 定义高频脉冲电流参数
total_time = 1000      # Total simulation time in seconds
time_step = 0.01          # 时间步长

# 创建时间数组
time = collect(Float64, 0:time_step:total_time)

# 多级脉冲参数（基于Tesla电池测试工况设计）
# Tesla电池测试特点：
# - 超级充电站V3：250kW功率，对应约2.5C充电倍率
# - 高速行驶工况：持续1-2C放电，间歇性4-6C加速
# - 热管理系统：需要在高频脉冲下维持电池温度在15-45°C
# - 电池管理系统：需要实时监控电压、温度、SOC变化
pulse_interval = 5   # 每12秒一个脉冲（模拟充电站间歇性高功率输出）
pulse_width = 0.5      # 脉冲持续时间5ms（模拟IGBT开关频率）
I_base = 5          # 基准电流(1C=5A，对应LG M50电池)

current = zeros(length(time))
for i in eachindex(time)
    t = time[i]
    
    # 判断是否在脉冲区间内
    if mod(t, pulse_interval) < pulse_width
        # 三级阶梯划分
        phase = mod(t, pulse_width) / pulse_width  # 脉冲内相对位置[0,1)
        
        if phase < 0.3
            # 第一级：3C上升沿
            current[i] = 3.0 * I_base * (1 + 0.2*sin(2π*1000*t))  # 叠加1kHz调制
        elseif phase < 0.6
            # 第二级：4C平台
            current[i] = 4.0 * I_base * (1 + 0.1*sin(2π*5000*t))  # 叠加5kHz调制
        else
            # 第三级：5C下降沿
            current[i] = 5.0 * I_base * exp(-(phase-0.6)^2/0.01)   # 高斯衰减
        end
        
    else
        # 脉冲间隔期：维持1.5C背景电流
        current[i] = 1.5 * I_base
    end
    
    # 添加宽频噪声（保留原有特征）
    current[i] += 0.05 * I_base * randn()
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

# 提取电化学相关数据
voltage_p2d = result["cell voltage [V]"]
concentration_p2d = result["negative particle surface lithium concentration [mol/m^3]"]
stress_p2d = result["negative particle surface tangential stress[Pa]"]
concentration_p_p2d = result["positive particle surface lithium concentration [mol/m^3]"]
stress_p_p2d = result["positive particle surface tangential stress[Pa]"]

# 提取反应电流、过电压和电势数据
reaction_current_p2d = result["negative electrode interfacial current density [A/m^2]"]
reaction_current_p_p2d = result["positive electrode interfacial current density [A/m^2]"]
overpotential_p2d = result["negative electrode overpotential [V]"]
overpotential_p_p2d = result["positive electrode overpotential [V]"]
potential_p2d = result["negative electrode open circuit potential [V]"]
potential_p_p2d = result["positive electrode open circuit potential [V]"]

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
# 选择3个关键时间点，时间越后颜色越深
analysis_times = [40.3, 500.0, 1000.2]
time_indices = []

for t in analysis_times
    idx = argmin(abs.(time_result .- t))
    push!(time_indices, idx)
end

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

# 提取sP2D模型的反应电流、过电压和电势数据
reaction_current_sP2D = result_sP2D["negative electrode interfacial current density [A/m^2]"]
reaction_current_p_sP2D = result_sP2D["positive electrode interfacial current density [A/m^2]"]
overpotential_sP2D = result_sP2D["negative electrode overpotential [V]"]
overpotential_p_sP2D = result_sP2D["positive electrode overpotential [V]"]
potential_sP2D = result_sP2D["negative electrode open circuit potential [V]"]
potential_p_sP2D = result_sP2D["positive electrode open circuit potential [V]"]

# 提取各时间点的空间分布数据
concentration_p2d_spatial = []
stress_p2d_spatial = []
concentration_p_p2d_spatial = []
stress_p_p2d_spatial = []

concentration_sP2D_spatial = []
stress_sP2D_spatial = []
concentration_p_sP2D_spatial = []
stress_p_sP2D_spatial = []

# 提取电化学参数的空间分布数据
reaction_current_p2d_spatial = []
reaction_current_p_p2d_spatial = []
overpotential_p2d_spatial = []
overpotential_p_p2d_spatial = []
potential_p2d_spatial = []
potential_p_p2d_spatial = []

reaction_current_sP2D_spatial = []
reaction_current_p_sP2D_spatial = []
overpotential_sP2D_spatial = []
overpotential_p_sP2D_spatial = []
potential_sP2D_spatial = []
potential_p_sP2D_spatial = []

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
    
    # 提取P2D模型电化学参数的空间分布
    if ndims(reaction_current_p2d) > 1
        push!(reaction_current_p2d_spatial, reaction_current_p2d[:, t_idx])
    end
    
    if ndims(reaction_current_p_p2d) > 1
        push!(reaction_current_p_p2d_spatial, reaction_current_p_p2d[:, t_idx])
    end
    
    if ndims(overpotential_p2d) > 1
        push!(overpotential_p2d_spatial, overpotential_p2d[:, t_idx])
    end
    
    if ndims(overpotential_p_p2d) > 1
        push!(overpotential_p_p2d_spatial, overpotential_p_p2d[:, t_idx])
    end
    
    if ndims(potential_p2d) > 1
        push!(potential_p2d_spatial, potential_p2d[:, t_idx])
    end
    
    if ndims(potential_p_p2d) > 1
        push!(potential_p_p2d_spatial, potential_p_p2d[:, t_idx])
    end
    
    # 提取sP2D模型电化学参数的空间分布
    if ndims(reaction_current_sP2D) > 1
        push!(reaction_current_sP2D_spatial, reaction_current_sP2D[:, t_idx])
    end
    
    if ndims(reaction_current_p_sP2D) > 1
        push!(reaction_current_p_sP2D_spatial, reaction_current_p_sP2D[:, t_idx])
    end
    
    if ndims(overpotential_sP2D) > 1
        push!(overpotential_sP2D_spatial, overpotential_sP2D[:, t_idx])
    end
    
    if ndims(overpotential_p_sP2D) > 1
        push!(overpotential_p_sP2D_spatial, overpotential_p_sP2D[:, t_idx])
    end
    
    if ndims(potential_sP2D) > 1
        push!(potential_sP2D_spatial, potential_sP2D[:, t_idx])
    end
    
    if ndims(potential_p_sP2D) > 1
        push!(potential_p_sP2D_spatial, potential_p_sP2D[:, t_idx])
    end
end

# 调试：检查数组长度
println("调试信息 - 数组长度检查:")
println("x_grid_n_normalized 长度: ", length(x_grid_n_normalized))
println("x_grid_p_normalized 长度: ", length(x_grid_p_normalized))
println("concentration_p2d_spatial[1] 长度: ", length(concentration_p2d_spatial[1]))
println("concentration_p_p2d_spatial[1] 长度: ", length(concentration_p_p2d_spatial[1]))
println("stress_p2d_spatial[1] 长度: ", length(stress_p2d_spatial[1]))
println("stress_p_p2d_spatial[1] 长度: ", length(stress_p_p2d_spatial[1]))

# 确保所有数组长度一致
if length(concentration_p2d_spatial[1]) != length(x_grid_n_normalized)
    println("警告：负极浓度数组长度不匹配！")
    println("x_grid_n_normalized: ", length(x_grid_n_normalized))
    println("concentration_p2d_spatial[1]: ", length(concentration_p2d_spatial[1]))
end

if length(concentration_p_p2d_spatial[1]) != length(x_grid_p_normalized)
    println("警告：正极浓度数组长度不匹配！")
    println("x_grid_p_normalized: ", length(x_grid_p_normalized))
    println("concentration_p_p2d_spatial[1]: ", length(concentration_p_p2d_spatial[1]))
end

# 定义科学论文风格的配色方案
# 使用蓝色系，时间越后颜色越深
time_colors = [
    :royalblue, :blue4, :purple  # 从浅到深的蓝色
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
         title="Electrode Concentration Distribution(Enertech)",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的浓度分布
for i in 1:length(time_indices)
    # P2D模型用实线
    # 正极浓度 - P2D (实线)
    plot!(x_grid_p_normalized, concentration_p_p2d_spatial[i], 
        label="P2D (Solid Line) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # 负极浓度 - P2D (实线)
    plot!(x_grid_n_normalized, concentration_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # sP2D模型用五角星标记 + 虚线连接
    # 正极浓度 - sP2D (五角星 + 虚线)
    plot!(x_grid_p_normalized, concentration_p_sP2D_spatial[i], 
        label="sP2D (Star+Dash) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)
    
    # 负极浓度 - sP2D (五角星 + 虚线)
    plot!(x_grid_n_normalized, concentration_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

end

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域填充
separator_mid = (negative_end + positive_start) / 2
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="Separator")

# 区域标识 - 浓度图专用
neg_center_conc = negative_end / 2
pos_center_conc = (positive_start + positive_end) / 2
annotate!(neg_center_conc, maximum([maximum(conc) for conc in concentration_p2d_spatial])*1.1, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_conc, maximum([maximum(conc) for conc in concentration_p_p2d_spatial])*0.3, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(conc) for conc in concentration_p2d_spatial])*0.80, 
         text("Separator", 9, :center, :italic, :gray))

# 创建应力空间分布图
p3 = plot(xlabel="Normalized Position", 
         ylabel="σ/E (dimensionless)", 
         title="Electrode Stress Distribution(Enertech)",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的应力分布
for i in 1:length(time_indices)
    # P2D模型用实线
    # 正极应力 - P2D (实线)
    plot!(x_grid_p_normalized, stress_p_p2d_spatial[i], 
        label="P2D (Solid Line) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # 负极应力 - P2D (实线)
    plot!(x_grid_n_normalized, stress_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # sP2D模型用五角星标记 + 虚线连接
    # 正极应力 - sP2D (五角星 + 虚线)
    plot!(x_grid_p_normalized, stress_p_sP2D_spatial[i], 
        label="sP2D (Star+Dash) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

    
    # 负极应力 - sP2D (五角星 + 虚线)
    plot!(x_grid_n_normalized, stress_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

end

# 添加分隔线和区域标识  
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 应力图专用
neg_center_stress = negative_end / 2
pos_center_stress = (positive_start + positive_end) / 2
annotate!(neg_center_stress, maximum([maximum(stress) for stress in stress_p2d_spatial])*1.05, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_stress, maximum([maximum(stress) for stress in stress_p_p2d_spatial])*2.0, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(stress) for stress in stress_p2d_spatial])*0.65, 
         text("Separator", 9, :center, :italic, :gray))

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

# 创建反应电流空间分布图
p7 = plot(xlabel="Normalized Position", 
         ylabel="Reaction Current Density [A/m²]", 
         title="Reaction Current Distribution(Enertech)",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的反应电流分布
for i in 1:length(time_indices)
    # P2D模型用实线
    # 正极反应电流 - P2D (实线)
    plot!(x_grid_p_normalized, reaction_current_p_p2d_spatial[i], 
        label="P2D (Solid Line) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # 负极反应电流 - P2D (实线)
    plot!(x_grid_n_normalized, reaction_current_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # sP2D模型用五角星标记 + 虚线连接
    # 正极反应电流 - sP2D (五角星 + 虚线)
    plot!(x_grid_p_normalized, reaction_current_p_sP2D_spatial[i], 
        label="sP2D (Star+Dash) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

    
    # 负极反应电流 - sP2D (五角星 + 虚线)
    plot!(x_grid_n_normalized, reaction_current_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="Separator")

# 区域标识 - 反应电流图专用
neg_center_reaction = negative_end / 2
pos_center_reaction = (positive_start + positive_end) / 2
annotate!(neg_center_reaction, maximum([maximum(curr) for curr in reaction_current_p2d_spatial])*0.95, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_reaction, maximum([maximum(curr) for curr in reaction_current_p_p2d_spatial])*0.1, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(curr) for curr in reaction_current_p2d_spatial])*0.2, 
         text("Separator", 9, :center, :italic, :gray))

# 创建过电压空间分布图
p8 = plot(xlabel="Normalized Position", 
         ylabel="Overpotential [V]", 
         title="Overpotential Distribution(Enertech)",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的过电压分布
for i in 1:length(time_indices)
    # P2D模型用实线
    # 正极过电压 - P2D (实线)
    plot!(x_grid_p_normalized, overpotential_p_p2d_spatial[i], 
        label="P2D (Solid Line) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # 负极过电压 - P2D (实线)
    plot!(x_grid_n_normalized, overpotential_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # sP2D模型用五角星标记 + 虚线连接
    # 正极过电压 - sP2D (五角星 + 虚线)
    plot!(x_grid_p_normalized, overpotential_p_sP2D_spatial[i], 
        label="sP2D (Star+Dash) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

    
    # 负极过电压 - sP2D (五角星 + 虚线)
    plot!(x_grid_n_normalized, overpotential_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8, markersize=6.5, marker=:star5, markerstrokewidth=0)

end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 过电压图专用
neg_center_over = negative_end / 2
pos_center_over = (positive_start + positive_end) / 2
annotate!(neg_center_over, maximum([maximum(over) for over in overpotential_p2d_spatial])*0.2, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_over, maximum([maximum(over) for over in overpotential_p_p2d_spatial])*0.05, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(over) for over in overpotential_p2d_spatial])*0.4, 
         text("Separator", 9, :center, :italic, :gray))

# 创建电势空间分布图
p9 = plot(xlabel="Normalized Position", 
         ylabel="Electrode Potential [V]", 
         title="Electrode Potential Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的电势分布
for i in 1:length(time_indices)
    # P2D模型用实线
    # 正极电势 - P2D (实线)
    plot!(x_grid_p_normalized, potential_p_p2d_spatial[i], 
        label="P2D (Solid Line) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # 负极电势 - P2D (实线)
    plot!(x_grid_n_normalized, potential_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:solid, alpha=0.8)
    
    # sP2D模型用五角星标记 + 虚线连接
    # 正极电势 - sP2D (五角星 + 虚线)
    plot!(x_grid_p_normalized, potential_p_sP2D_spatial[i], 
        label="sP2D (Star+Dash) - $(analysis_times[i])s", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8)
    scatter!(x_grid_p_normalized, potential_p_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors[i], marker=:star5, alpha=0.8)
    
    # 负极电势 - sP2D (五角星 + 虚线)
    plot!(x_grid_n_normalized, potential_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors[i], linestyle=:dash, alpha=0.8)
    scatter!(x_grid_n_normalized, potential_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors[i], marker=:star5, alpha=0.8)
end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 电势图专用
neg_center_pot = negative_end / 2
pos_center_pot = (positive_start + positive_end) / 2
annotate!(neg_center_pot, maximum([maximum(pot) for pot in potential_p2d_spatial])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_pot, maximum([maximum(pot) for pot in potential_p_p2d_spatial])*1.2, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(pot) for pot in potential_p2d_spatial])*0.85, 
         text("Separator", 8, :center, :italic, :gray))

# 组合图表
plot_combined = plot(p1, p2, p3, p4, p7, p8, p9, 
                    layout=(7, 1), 
                    size=(1000, 1800))

# 保存图表
savefig(p1, "pulse_voltage_spatial_comparison.png")
savefig(p2, "pulse_concentration_spatial_comparison.png")
savefig(p3, "pulse_stress_spatial_comparison.png")
savefig(p4, "pulse_current_profile_spatial.png")
savefig(p7, "pulse_reaction_current_spatial_comparison.png")
savefig(p8, "pulse_overpotential_spatial_comparison.png")
savefig(p9, "pulse_potential_spatial_comparison.png")
savefig(plot_combined, "pulse_spatial_analysis_comparison_super.png")

# 计算空间分布误差
concentration_spatial_errors_n = []  # 负极浓度误差
stress_spatial_errors_n = []         # 负极应力误差
concentration_spatial_errors_p = []  # 正极浓度误差
stress_spatial_errors_p = []         # 正极应力误差

# 计算电化学参数的空间分布误差
reaction_current_spatial_errors_n = []  # 负极反应电流误差
reaction_current_spatial_errors_p = []  # 正极反应电流误差
overpotential_spatial_errors_n = []     # 负极过电压误差
overpotential_spatial_errors_p = []     # 正极过电压误差
potential_spatial_errors_n = []         # 负极电势误差
potential_spatial_errors_p = []         # 正极电势误差

for i in 1:length(time_indices)
    # 负极浓度相对误差 (%)
    conc_error_n = abs.(concentration_sP2D_spatial[i] .- concentration_p2d_spatial[i]) ./ abs.(concentration_p2d_spatial[i]) .* 100
    push!(concentration_spatial_errors_n, conc_error_n)
    
    # 负极应力相对误差 (%)
    stress_error_n = abs.(stress_sP2D_spatial[i] .- stress_p2d_spatial[i]) ./ abs.(stress_p2d_spatial[i]) .* 100
    push!(stress_spatial_errors_n, stress_error_n)
    
    # 正极浓度相对误差 (%)
    conc_error_p = abs.(concentration_p_sP2D_spatial[i] .- concentration_p_p2d_spatial[i]) ./ abs.(concentration_p_p2d_spatial[i]) .* 100
    push!(concentration_spatial_errors_p, conc_error_p)
    
    # 正极应力相对误差 (%)
    stress_error_p = abs.(stress_p_sP2D_spatial[i] .- stress_p_p2d_spatial[i]) ./ abs.(stress_p_p2d_spatial[i]) .* 100
    push!(stress_spatial_errors_p, stress_error_p)
    
    # 负极反应电流相对误差 (%)
    reaction_error_n = abs.(reaction_current_sP2D_spatial[i] .- reaction_current_p2d_spatial[i]) ./ abs.(reaction_current_p2d_spatial[i]) .* 100
    push!(reaction_current_spatial_errors_n, reaction_error_n)
    
    # 正极反应电流相对误差 (%)
    reaction_error_p = abs.(reaction_current_p_sP2D_spatial[i] .- reaction_current_p_p2d_spatial[i]) ./ abs.(reaction_current_p_p2d_spatial[i]) .* 100
    push!(reaction_current_spatial_errors_p, reaction_error_p)
    
    # 负极过电压相对误差 (%)
    over_error_n = abs.(overpotential_sP2D_spatial[i] .- overpotential_p2d_spatial[i]) ./ abs.(overpotential_p2d_spatial[i]) .* 100
    push!(overpotential_spatial_errors_n, over_error_n)
    
    # 正极过电压相对误差 (%)
    over_error_p = abs.(overpotential_p_sP2D_spatial[i] .- overpotential_p_p2d_spatial[i]) ./ abs.(overpotential_p_p2d_spatial[i]) .* 100
    push!(overpotential_spatial_errors_p, over_error_p)
    
    # 负极电势相对误差 (%)
    pot_error_n = abs.(potential_sP2D_spatial[i] .- potential_p2d_spatial[i]) ./ abs.(potential_p2d_spatial[i]) .* 100
    push!(potential_spatial_errors_n, pot_error_n)
    
    # 正极电势相对误差 (%)
    pot_error_p = abs.(potential_p_sP2D_spatial[i] .- potential_p_p2d_spatial[i]) ./ abs.(potential_p_p2d_spatial[i]) .* 100
    push!(potential_spatial_errors_p, pot_error_p)
end

# 创建空间误差分布图 - 浓度误差
p5 = plot(xlabel="Normalized Position", 
         ylabel="Concentration Relative Error (%)", 
         title="Spatial Concentration Relative Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用实线
for i in 1:length(time_indices)
    plot!(x_grid_n_normalized, concentration_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, concentration_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 浓度误差图专用
neg_center_conc_err = negative_end / 2
pos_center_conc_err = (positive_start + positive_end) / 2
annotate!(neg_center_conc_err, maximum([maximum(err) for err in concentration_spatial_errors_n])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_conc_err, maximum([maximum(err) for err in concentration_spatial_errors_p])*4.0, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in concentration_spatial_errors_n])*0.7, 
         text("Separator", 8, :center, :italic, :gray))

# 创建空间误差分布图 - 应力误差
p6 = plot(xlabel="Normalized Position", 
         ylabel="Stress Relative Error (%)", 
         title="Spatial Stress Relative Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用实线
for i in 1:length(time_indices)
    plot!(x_grid_n_normalized, stress_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, stress_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 应力误差图专用
neg_center_stress_err = negative_end / 2
pos_center_stress_err = (positive_start + positive_end) / 2
annotate!(neg_center_stress_err, maximum([maximum(err) for err in stress_spatial_errors_n])*0.8, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_stress_err, maximum([maximum(err) for err in stress_spatial_errors_p])*2.0, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in stress_spatial_errors_n])*0.5, 
         text("Separator", 8, :center, :italic, :gray))

savefig(p5, "pulse_concentration_spatial_relative_error.png")
savefig(p6, "pulse_stress_spatial_relative_error.png")

# 创建电化学参数误差分布图
# 反应电流误差分布图
p10 = plot(xlabel="Normalized Position", 
         ylabel="Reaction Current Relative Error (%)", 
         title="Spatial Reaction Current Relative Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用实线
for i in 1:length(time_indices)
    plot!(x_grid_n_normalized, reaction_current_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, reaction_current_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 反应电流误差图专用
neg_center_reaction_err = negative_end / 2
pos_center_reaction_err = (positive_start + positive_end) / 2
annotate!(neg_center_reaction_err, maximum([maximum(err) for err in reaction_current_spatial_errors_n])*0.5, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_reaction_err, maximum([maximum(err) for err in reaction_current_spatial_errors_p])*2.0, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in reaction_current_spatial_errors_n])*0.8, 
         text("Separator", 8, :center, :italic, :gray))

# 过电压误差分布图
p11 = plot(xlabel="Normalized Position", 
         ylabel="Overpotential Relative Error (%)", 
         title="Spatial Overpotential Relative Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用实线
for i in 1:length(time_indices)
    plot!(x_grid_n_normalized, overpotential_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, overpotential_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 过电压误差图专用
neg_center_over_err = negative_end / 2
pos_center_over_err = (positive_start + positive_end) / 2
annotate!(neg_center_over_err, maximum([maximum(err) for err in overpotential_spatial_errors_n])*0.5, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_over_err, maximum([maximum(err) for err in overpotential_spatial_errors_p])*3.0, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in overpotential_spatial_errors_n])*0.4, 
         text("Separator", 8, :center, :italic, :gray))

# 电势误差分布图
p12 = plot(xlabel="Normalized Position", 
         ylabel="Potential Relative Error (%)", 
         title="Spatial Potential Relative Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制负极误差 - 使用实线
for i in 1:length(time_indices)
    plot!(x_grid_n_normalized, potential_spatial_errors_n[i], 
        label="Negative t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, potential_spatial_errors_p[i], 
        label="Positive t=$(analysis_times[i])s", 
        lw=2.0, 
        color=time_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 电势误差图专用
neg_center_pot_err = negative_end / 2
pos_center_pot_err = (positive_start + positive_end) / 2
annotate!(neg_center_pot_err, maximum([maximum(err) for err in potential_spatial_errors_n])*0.95, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_pot_err, maximum([maximum(err) for err in potential_spatial_errors_p])*1.2, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in potential_spatial_errors_n])*0.85, 
         text("Separator", 8, :center, :italic, :gray))

savefig(p10, "pulse_reaction_current_spatial_relative_error.png")
savefig(p11, "pulse_overpotential_spatial_relative_error.png")
savefig(p12, "pulse_potential_spatial_relative_error.png")

# 计算和打印统计信息
max_spatial_conc_errors_n = [maximum(err) for err in concentration_spatial_errors_n]
max_spatial_stress_errors_n = [maximum(err) for err in stress_spatial_errors_n]
mean_spatial_conc_errors_n = [mean(err) for err in concentration_spatial_errors_n]
mean_spatial_stress_errors_n = [mean(err) for err in stress_spatial_errors_n]

max_spatial_conc_errors_p = [maximum(err) for err in concentration_spatial_errors_p]
max_spatial_stress_errors_p = [maximum(err) for err in stress_spatial_errors_p]
mean_spatial_conc_errors_p = [mean(err) for err in concentration_spatial_errors_p]
mean_spatial_stress_errors_p = [mean(err) for err in stress_spatial_errors_p]

# 电化学参数误差统计
max_spatial_reaction_errors_n = [maximum(err) for err in reaction_current_spatial_errors_n]
max_spatial_over_errors_n = [maximum(err) for err in overpotential_spatial_errors_n]
max_spatial_pot_errors_n = [maximum(err) for err in potential_spatial_errors_n]
mean_spatial_reaction_errors_n = [mean(err) for err in reaction_current_spatial_errors_n]
mean_spatial_over_errors_n = [mean(err) for err in overpotential_spatial_errors_n]
mean_spatial_pot_errors_n = [mean(err) for err in potential_spatial_errors_n]

max_spatial_reaction_errors_p = [maximum(err) for err in reaction_current_spatial_errors_p]
max_spatial_over_errors_p = [maximum(err) for err in overpotential_spatial_errors_p]
max_spatial_pot_errors_p = [maximum(err) for err in potential_spatial_errors_p]
mean_spatial_reaction_errors_p = [mean(err) for err in reaction_current_spatial_errors_p]
mean_spatial_over_errors_p = [mean(err) for err in overpotential_spatial_errors_p]
mean_spatial_pot_errors_p = [mean(err) for err in potential_spatial_errors_p]

# 合并正负极误差数据用于总体统计
all_conc_errors = vcat(concentration_spatial_errors_n, concentration_spatial_errors_p)
all_stress_errors = vcat(stress_spatial_errors_n, stress_spatial_errors_p)
all_reaction_errors = vcat(reaction_current_spatial_errors_n, reaction_current_spatial_errors_p)
all_over_errors = vcat(overpotential_spatial_errors_n, overpotential_spatial_errors_p)
all_pot_errors = vcat(potential_spatial_errors_n, potential_spatial_errors_p)

# ===== sP2D简化机制误差分析 =====
println("\n" * "="^60)
println("sP2D简化机制误差分析")
println("="^60)

# 1. 过电压误差与应力误差的相关性分析
println("\n1. 过电压误差与应力误差的相关性分析")
println("-"^50)

for i in 1:length(time_indices)
    t = analysis_times[i]
    
    # 负极相关性
    if length(overpotential_spatial_errors_n[i]) > 1 && length(stress_spatial_errors_n[i]) > 1
        cor_over_stress_n = cor(overpotential_spatial_errors_n[i], stress_spatial_errors_n[i])
        println("时间点 $(t)s - 负极:")
        println("  过电压误差与应力误差相关系数: $(round(cor_over_stress_n, digits=4))")
        
        # 分析误差分布特征
        max_over_err_n = maximum(overpotential_spatial_errors_n[i])
        max_stress_err_n = maximum(stress_spatial_errors_n[i])
        println("  最大过电压相对误差: $(round(max_over_err_n, digits=2))%")
        println("  最大应力相对误差: $(round(max_stress_err_n, digits=2))%")
        println("  误差比值 (应力/过电压): $(round(max_stress_err_n/max_over_err_n, digits=2))")
    end
    
    # 正极相关性
    if length(overpotential_spatial_errors_p[i]) > 1 && length(stress_spatial_errors_p[i]) > 1
        cor_over_stress_p = cor(overpotential_spatial_errors_p[i], stress_spatial_errors_p[i])
        println("时间点 $(t)s - 正极:")
        println("  过电压误差与应力误差相关系数: $(round(cor_over_stress_p, digits=4))")
        
        max_over_err_p = maximum(overpotential_spatial_errors_p[i])
        max_stress_err_p = maximum(stress_spatial_errors_p[i])
        println("  最大过电压相对误差: $(round(max_over_err_p, digits=2))%")
        println("  最大应力相对误差: $(round(max_stress_err_p, digits=2))%")
        println("  误差比值 (应力/过电压): $(round(max_stress_err_p/max_over_err_p, digits=2))")
    end
    println()
end

# 2. 反应电流误差与浓度误差的相关性分析
println("\n2. 反应电流误差与浓度误差的相关性分析")
println("-"^50)

for i in 1:length(time_indices)
    t = analysis_times[i]
    
    # 负极相关性
    if length(reaction_current_spatial_errors_n[i]) > 1 && length(concentration_spatial_errors_n[i]) > 1
        cor_reaction_conc_n = cor(reaction_current_spatial_errors_n[i], concentration_spatial_errors_n[i])
        println("时间点 $(t)s - 负极:")
        println("  反应电流误差与浓度误差相关系数: $(round(cor_reaction_conc_n, digits=4))")
        
        max_reaction_err_n = maximum(reaction_current_spatial_errors_n[i])
        max_conc_err_n = maximum(concentration_spatial_errors_n[i])
        println("  最大反应电流相对误差: $(round(max_reaction_err_n, digits=2))%")
        println("  最大浓度相对误差: $(round(max_conc_err_n, digits=2))%")
    end
    
    # 正极相关性
    if length(reaction_current_spatial_errors_p[i]) > 1 && length(concentration_spatial_errors_p[i]) > 1
        cor_reaction_conc_p = cor(reaction_current_spatial_errors_p[i], concentration_spatial_errors_p[i])
        println("时间点 $(t)s - 正极:")
        println("  反应电流误差与浓度误差相关系数: $(round(cor_reaction_conc_p, digits=4))")
        
        max_reaction_err_p = maximum(reaction_current_spatial_errors_p[i])
        max_conc_err_p = maximum(concentration_spatial_errors_p[i])
        println("  最大反应电流相对误差: $(round(max_reaction_err_p, digits=2))%")
        println("  最大浓度相对误差: $(round(max_conc_err_p, digits=2))%")
    end
    println()
end

# 3. sP2D简化机制误差解释
println("\n3. sP2D简化机制误差解释")
println("-"^50)

println("过电压误差产生原因:")
println("  A. 电势分布简化:")
println("     - 边界层效应: sP2D用指数函数拟合电极表面快速电势变化")
println("     - 内部区域: 用多项式拟合，忽略局部细节变化")
println("     - 区域边界: 不同拟合函数在边界处可能不连续")
println()
println("  B. 浓度-电势耦合简化:")
println("     - P2D: ∇²φe = -∇·(κ∇φe) + ∇·(κ∇ln(ce)) (完全耦合)")
println("     - sP2D: 解耦处理，忽略浓度梯度对电势的影响")
println("     - 迁移-扩散耦合被简化处理")
println()

println("反应电流误差产生原因:")
println("  A. 过电压误差传递:")
println("     - 反应电流: j = j₀[exp(αₐFη/RT) - exp(-αcFη/RT)]")
println("     - 过电压小误差被指数函数放大")
println("     - η_error → j_error (误差传递链)")
println()
println("  B. 交换电流密度简化:")
println("     - P2D: j₀ = k(ce)^0.5(cs)^0.5(1-cs)^0.5 (精确计算)")
println("     - sP2D: j₀ ≈ k(ce_avg)^0.5(cs_avg)^0.5(1-cs_avg)^0.5 (平均值)")
println("     - 忽略局部浓度梯度和反应活性空间分布")
println()

# 4. 高频脉冲工况下的特殊影响
println("\n4. 高频脉冲工况下的特殊影响")
println("-"^50)

println("高频脉冲对sP2D误差的影响:")
println("  A. 快速变化挑战:")
println("     - 脉冲电流快速变化 → 电势响应滞后")
println("     - sP2D简化模型响应速度可能跟不上P2D")
println("     - 边界层厚度动态变化，固定拟合参数不够精确")
println()
println("  B. 空间-时间耦合:")
println("     - 高频脉冲产生复杂的空间-时间耦合效应")
println("     - sP2D的区域划分可能在高频下不够精细")
println("     - 不同时间点的误差模式可能不同")
println()

# 5. 误差随时间的演化分析
println("\n5. 误差随时间的演化分析")
println("-"^50)

overpotential_error_trend = [mean(err) for err in overpotential_spatial_errors_n]
reaction_current_error_trend = [mean(err) for err in reaction_current_spatial_errors_n]
stress_error_trend = [mean(err) for err in stress_spatial_errors_n]
concentration_error_trend = [mean(err) for err in concentration_spatial_errors_n]

println("平均误差随时间变化趋势:")
for i in 1:length(analysis_times)
    t = analysis_times[i]
    println("时间点 $(t)s:")
    println("  过电压相对误差: $(round(overpotential_error_trend[i], digits=2))%")
    println("  反应电流相对误差: $(round(reaction_current_error_trend[i], digits=2))%")
    println("  应力相对误差: $(round(stress_error_trend[i], digits=2))%")
    println("  浓度相对误差: $(round(concentration_error_trend[i], digits=2))%")
end

# 计算误差增长率
if length(overpotential_error_trend) > 1
    over_error_growth = (overpotential_error_trend[end] - overpotential_error_trend[1]) / overpotential_error_trend[1] * 100
    reaction_error_growth = (reaction_current_error_trend[end] - reaction_current_error_trend[1]) / reaction_current_error_trend[1] * 100
    stress_error_growth = (stress_error_trend[end] - stress_error_trend[1]) / stress_error_trend[1] * 100
    conc_error_growth = (concentration_error_trend[end] - concentration_error_trend[1]) / concentration_error_trend[1] * 100
    
    println("\n相对误差增长率 (从$(analysis_times[1])s到$(analysis_times[end])s):")
    println("  过电压相对误差增长率: $(round(over_error_growth, digits=1))%")
    println("  反应电流相对误差增长率: $(round(reaction_error_growth, digits=1))%")
    println("  应力相对误差增长率: $(round(stress_error_growth, digits=1))%")
    println("  浓度相对误差增长率: $(round(conc_error_growth, digits=1))%")
end

println("\n" * "="^60)

# ===== 原始统计信息输出 =====
println("\n=== 高频脉冲电流工况空间分布分析结果 ===")
println("分析时间点: ", analysis_times, " 秒")
println("负极最大空间浓度误差: ", max_spatial_conc_errors_n, " (c/c_max)")
println("负极平均空间浓度误差: ", mean_spatial_conc_errors_n, " (c/c_max)")
println("负极最大空间应力误差: ", max_spatial_stress_errors_n, " (σ/E)")
println("负极平均空间应力误差: ", mean_spatial_stress_errors_n, " (σ/E)")
println("正极最大空间浓度误差: ", max_spatial_conc_errors_p, " (c/c_max)")
println("正极平均空间浓度误差: ", mean_spatial_conc_errors_p, " (c/c_max)")
println("正极最大空间应力误差: ", max_spatial_stress_errors_p, " (σ/E)")
println("正极平均空间应力误差: ", mean_spatial_stress_errors_p, " (σ/E)")

println("\n=== 电化学参数误差统计 ===")
println("负极最大反应电流误差: ", max_spatial_reaction_errors_n, " [A/m²]")
println("负极平均反应电流误差: ", mean_spatial_reaction_errors_n, " [A/m²]")
println("负极最大过电压误差: ", max_spatial_over_errors_n, " [V]")
println("负极平均过电压误差: ", mean_spatial_over_errors_n, " [V]")
println("负极最大电势误差: ", max_spatial_pot_errors_n, " [V]")
println("负极平均电势误差: ", mean_spatial_pot_errors_n, " [V]")

println("正极最大反应电流误差: ", max_spatial_reaction_errors_p, " [A/m²]")
println("正极平均反应电流误差: ", mean_spatial_reaction_errors_p, " [A/m²]")
println("正极最大过电压误差: ", max_spatial_over_errors_p, " [V]")
println("正极平均过电压误差: ", mean_spatial_over_errors_p, " [V]")
println("正极最大电势误差: ", max_spatial_pot_errors_p, " [V]")
println("正极平均电势误差: ", mean_spatial_pot_errors_p, " [V]")

# 总体最大误差
overall_max_conc_error = maximum([maximum(err) for err in all_conc_errors])
overall_max_stress_error = maximum([maximum(err) for err in all_stress_errors])
overall_max_reaction_error = maximum([maximum(err) for err in all_reaction_errors])
overall_max_over_error = maximum([maximum(err) for err in all_over_errors])
overall_max_pot_error = maximum([maximum(err) for err in all_pot_errors])

overall_mean_conc_error = mean([mean(err) for err in all_conc_errors])
overall_mean_stress_error = mean([mean(err) for err in all_stress_errors])
overall_mean_reaction_error = mean([mean(err) for err in all_reaction_errors])
overall_mean_over_error = mean([mean(err) for err in all_over_errors])
overall_mean_pot_error = mean([mean(err) for err in all_pot_errors])

println("\n=== 总体统计 ===")
println("总体最大浓度误差: $overall_max_conc_error (c/c_max)")
println("总体平均浓度误差: $overall_mean_conc_error (c/c_max)")
println("总体最大应力误差: $overall_max_stress_error (σ/E)")
println("总体平均应力误差: $overall_mean_stress_error (σ/E)")
println("总体最大反应电流误差: $overall_max_reaction_error [A/m²]")
println("总体平均反应电流误差: $overall_mean_reaction_error [A/m²]")
println("总体最大过电压误差: $overall_max_over_error [V]")
println("总体平均过电压误差: $overall_mean_over_error [V]")
println("总体最大电势误差: $overall_max_pot_error [V]")
println("总体平均电势误差: $overall_mean_pot_error [V]")

# 分析高频特性对空间分布的影响
println("\n=== 高频特性分析 ===")
high_freq_indices = [1, 3]  # 选择前、后时间点分析高频特性
low_freq_indices = [2]  # 选择中间时间点作为对比

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
    
    # 添加电化学参数
    negative_df[!, "ReactionCurrent_P2D_t$(t)s_[A_m2]"] = reaction_current_p2d_spatial[i]
    negative_df[!, "ReactionCurrent_sP2D_t$(t)s_[A_m2]"] = reaction_current_sP2D_spatial[i]
    negative_df[!, "ReactionCurrent_Error_t$(t)s_[A_m2]"] = reaction_current_spatial_errors_n[i]
    negative_df[!, "Overpotential_P2D_t$(t)s_[V]"] = overpotential_p2d_spatial[i]
    negative_df[!, "Overpotential_sP2D_t$(t)s_[V]"] = overpotential_sP2D_spatial[i]
    negative_df[!, "Overpotential_Error_t$(t)s_[V]"] = overpotential_spatial_errors_n[i]
    negative_df[!, "Potential_P2D_t$(t)s_[V]"] = potential_p2d_spatial[i]
    negative_df[!, "Potential_sP2D_t$(t)s_[V]"] = potential_sP2D_spatial[i]
    negative_df[!, "Potential_Error_t$(t)s_[V]"] = potential_spatial_errors_n[i]
    
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
    
    # 添加电化学参数
    positive_df[!, "ReactionCurrent_P2D_t$(t)s_[A_m2]"] = reaction_current_p_p2d_spatial[i]
    positive_df[!, "ReactionCurrent_sP2D_t$(t)s_[A_m2]"] = reaction_current_p_sP2D_spatial[i]
    positive_df[!, "ReactionCurrent_Error_t$(t)s_[A_m2]"] = reaction_current_spatial_errors_p[i]
    positive_df[!, "Overpotential_P2D_t$(t)s_[V]"] = overpotential_p_p2d_spatial[i]
    positive_df[!, "Overpotential_sP2D_t$(t)s_[V]"] = overpotential_p_sP2D_spatial[i]
    positive_df[!, "Overpotential_Error_t$(t)s_[V]"] = overpotential_spatial_errors_p[i]
    positive_df[!, "Potential_P2D_t$(t)s_[V]"] = potential_p_p2d_spatial[i]
    positive_df[!, "Potential_sP2D_t$(t)s_[V]"] = potential_p_sP2D_spatial[i]
    positive_df[!, "Potential_Error_t$(t)s_[V]"] = potential_spatial_errors_p[i]
    
    # 为每个位置添加对应的电流值（重复标量值）
    positive_df[!, "Current_t$(t)s"] = fill(current_interp(t), length(x_grid_p))
end

# 保存两个DataFrame
CSV.write("pulse_spatial_analysis_negative_electrode.csv", negative_df)
CSV.write("pulse_spatial_analysis_positive_electrode.csv", positive_df)

println("高频脉冲电流空间分布分析完成，结果已保存")

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

# 修复RMSE计算部分的代码

# ===== RMSE计算与输出 =====
println("\n" * "="^60)
println("RMSE (均方根误差) 详细分析")
println("="^60)

# 定义RMSE计算函数
function calculate_rmse(predicted, actual)
    squared_errors = (predicted .- actual).^2
    mean_squared_error = mean(squared_errors)
    return sqrt(mean_squared_error)
end

# 存储各参数的RMSE结果
rmse_results = Dict()

# 计算各时间点的负极浓度RMSE
rmse_conc_n = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(concentration_sP2D_spatial[i], concentration_p2d_spatial[i])
    push!(rmse_conc_n, rmse)
end
rmse_results["负极浓度"] = rmse_conc_n

# 计算各时间点的正极浓度RMSE
rmse_conc_p = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(concentration_p_sP2D_spatial[i], concentration_p_p2d_spatial[i])
    push!(rmse_conc_p, rmse)
end
rmse_results["正极浓度"] = rmse_conc_p

# 计算各时间点的负极应力RMSE
rmse_stress_n = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(stress_sP2D_spatial[i], stress_p2d_spatial[i])
    push!(rmse_stress_n, rmse)
end
rmse_results["负极应力"] = rmse_stress_n

# 计算各时间点的正极应力RMSE
rmse_stress_p = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(stress_p_sP2D_spatial[i], stress_p_p2d_spatial[i])
    push!(rmse_stress_p, rmse)
end
rmse_results["正极应力"] = rmse_stress_p

# 计算各时间点的负极反应电流RMSE
rmse_reaction_n = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(reaction_current_sP2D_spatial[i], reaction_current_p2d_spatial[i])
    push!(rmse_reaction_n, rmse)
end
rmse_results["负极反应电流"] = rmse_reaction_n

# 计算各时间点的正极反应电流RMSE
rmse_reaction_p = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(reaction_current_p_sP2D_spatial[i], reaction_current_p_p2d_spatial[i])
    push!(rmse_reaction_p, rmse)
end
rmse_results["正极反应电流"] = rmse_reaction_p

# 计算各时间点的负极过电压RMSE
rmse_over_n = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(overpotential_sP2D_spatial[i], overpotential_p2d_spatial[i])
    push!(rmse_over_n, rmse)
end
rmse_results["负极过电压"] = rmse_over_n

# 计算各时间点的正极过电压RMSE
rmse_over_p = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(overpotential_p_sP2D_spatial[i], overpotential_p_p2d_spatial[i])
    push!(rmse_over_p, rmse)
end
rmse_results["正极过电压"] = rmse_over_p

# 计算各时间点的负极电势RMSE
rmse_pot_n = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(potential_sP2D_spatial[i], potential_p2d_spatial[i])
    push!(rmse_pot_n, rmse)
end
rmse_results["负极电势"] = rmse_pot_n

# 计算各时间点的正极电势RMSE
rmse_pot_p = []
for i in 1:length(time_indices)
    rmse = calculate_rmse(potential_p_sP2D_spatial[i], potential_p_p2d_spatial[i])
    push!(rmse_pot_p, rmse)
end
rmse_results["正极电势"] = rmse_pot_p

# 计算电压RMSE
rmse_voltage = calculate_rmse(voltage_sP2D, voltage_p2d)
rmse_results["电池电压"] = rmse_voltage

# 归一化RMSE
# 定义最大值用于归一化
max_values = Dict(
    "负极浓度" => c_max_neg,
    "正极浓度" => c_max_pos,
    "负极应力" => E_neg,
    "正极应力" => E_pos,
    "负极反应电流" => maximum([maximum(abs.(curr)) for curr in reaction_current_p2d_spatial]),
    "正极反应电流" => maximum([maximum(abs.(curr)) for curr in reaction_current_p_p2d_spatial]),
    "负极过电压" => maximum([maximum(abs.(over)) for over in overpotential_p2d_spatial]),
    "正极过电压" => maximum([maximum(abs.(over)) for over in overpotential_p_p2d_spatial]),
    "负极电势" => maximum([maximum(abs.(pot)) for pot in potential_p2d_spatial]),
    "正极电势" => maximum([maximum(abs.(pot)) for pot in potential_p_p2d_spatial]),
    "电池电压" => maximum(abs.(voltage_p2d))
)

# 打印详细的RMSE数据表格
println("\n电化学参数RMSE分析（按时间点）")
println("-"^80)
println("| 参数名称       | 单位       | ", join(["t=$(round(t, digits=1))s" |> (s -> rpad(s, 10)) for t in analysis_times], " | "), " | 平均RMSE   | 归一化RMSE |")
println("|", "-"^14, "|", "-"^11, "|", join(["-"^12 for _ in 1:length(analysis_times)], "|"), "|", "-"^11, "|", "-"^11, "|")

# 输出按参数归类的RMSE表格
for param in ["负极浓度", "正极浓度", "负极应力", "正极应力", 
             "负极反应电流", "正极反应电流", "负极过电压", "正极过电压", 
             "负极电势", "正极电势"]
    
    # 确定参数单位
    unit = if contains(param, "浓度")
        "mol/m³"
    elseif contains(param, "应力")
        "Pa"
    elseif contains(param, "反应电流")
        "A/m²"
    elseif contains(param, "过电压") || contains(param, "电势")
        "V"
    else
        ""
    end
    
    values = rmse_results[param]
    avg_rmse = mean(values)
    norm_rmse = avg_rmse / max_values[param]
    
    # 优雅地打印RMSE表格行
    if isa(values, Vector)
        print("| ", rpad(param, 14), "| ", rpad(unit, 10), "| ")
        for val in values
            print(rpad(string(round(val, digits=6)), 10), " | ")
        end
        println(rpad(string(round(avg_rmse, digits=6)), 10), " | ", rpad(string(round(norm_rmse, digits=6)), 10), " |")
    end
end

# 输出电压RMSE
println("| ", rpad("电池电压", 14), "| ", rpad("V", 10), "| ", "-"^(12*length(analysis_times)-1), " | ", 
        rpad(string(round(rmse_voltage, digits=6)), 10), " | ", 
        rpad(string(round(rmse_voltage/max_values["电池电压"], digits=6)), 10), " |")

println("-"^80)

# 计算总体RMSE指标
overall_rmse = Dict()
for param_type in ["浓度", "应力", "反应电流", "过电压", "电势"]
    neg_key = "负极" * param_type
    pos_key = "正极" * param_type
    
    # 计算负极和正极的平均RMSE
    neg_avg = mean(rmse_results[neg_key])
    pos_avg = mean(rmse_results[pos_key])
    
    # 归一化
    neg_norm = neg_avg / max_values[neg_key]
    pos_norm = pos_avg / max_values[pos_key]
    
    # 总体归一化RMSE (平均)
    overall_norm = (neg_norm + pos_norm) / 2
    
    overall_rmse[param_type] = overall_norm
end

# 输出总体归一化RMSE
println("\n总体归一化RMSE (平均值)")
println("-"^40)
for (param, value) in overall_rmse
    println("| ", rpad(param, 10), " | ", round(value, digits=16))
end
println("-"^40)

# 计算RMSE时间趋势
println("\nRMSE时间变化趋势分析")
println("-"^60)
for param in ["浓度", "应力", "反应电流", "过电压", "电势"]
    neg_key = "负极" * param
    pos_key = "正极" * param
    
    neg_values = rmse_results[neg_key]
    pos_values = rmse_results[pos_key]
    
    # 计算变化率 (最后点/首点 - 1) * 100%
    if length(neg_values) > 1
        neg_change = (neg_values[end] - neg_values[1]) / neg_values[1] * 100
        pos_change = (pos_values[end] - pos_values[1]) / pos_values[1] * 100
        
        println("$param RMSE变化率:")
        println("  负极: $(round(neg_change, digits=2))%")
        println("  正极: $(round(pos_change, digits=2))%")
    end
end

# 计算参数间RMSE相关性
println("\nRMSE参数间相关性分析")
println("-"^60)

# 特别分析以下关键关系对
correlation_pairs = [
    ("过电压", "反应电流"),
    ("反应电流", "浓度"),
    ("浓度", "应力"),
    ("过电压", "应力")
]

for (param1, param2) in correlation_pairs
    neg_key1 = "负极" * param1
    neg_key2 = "负极" * param2
    pos_key1 = "正极" * param1
    pos_key2 = "正极" * param2
    
    # 计算负极相关性
    neg_cor = cor(rmse_results[neg_key1], rmse_results[neg_key2])
    
    # 计算正极相关性
    pos_cor = cor(rmse_results[pos_key1], rmse_results[pos_key2])
    
    println("$param1-$param2 RMSE相关性:")
    println("  负极: $(round(neg_cor, digits=4))")
    println("  正极: $(round(pos_cor, digits=4))")
end

function correlation_std_error(r, n)
    # 相关系数的标准差近似公式
    se_r = sqrt((1 - r^2)/(n-2))
    return se_r
end

println("\nRMSE参数间相关性分析 (含标准差)")
println("-"^70)

for (param1, param2) in correlation_pairs
    neg_key1 = "负极" * param1
    neg_key2 = "负极" * param2
    pos_key1 = "正极" * param1
    pos_key2 = "正极" * param2
    
    # 计算负极相关性及标准差
    neg_cor = cor(rmse_results[neg_key1], rmse_results[neg_key2])
    neg_std = correlation_std_error(neg_cor, length(rmse_results[neg_key1]))
    
    # 计算正极相关性及标准差
    pos_cor = cor(rmse_results[pos_key1], rmse_results[pos_key2])
    pos_std = correlation_std_error(pos_cor, length(rmse_results[pos_key1]))
    
    println("$param1-$param2 RMSE相关性:")
    println("  负极: $(round(neg_cor, digits=4)) ± $(round(neg_std, digits=4))")
    println("  正极: $(round(pos_cor, digits=4)) ± $(round(pos_std, digits=4))")
end

println("\n" * "="^60)

# 引用信息
JuBat.Citation()
