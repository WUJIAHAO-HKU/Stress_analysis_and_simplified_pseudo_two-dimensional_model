using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
include("../src/JuBat.jl") 

# 设置数据路径
path = "D:/竞赛和课程文件/课程文件/毕业设计/SRC/data/drive_cycles/"

# 读取UDDS数据
udds_data = CSV.read(path * "UDDS2.csv", DataFrame, header = 1)
time_udds = udds_data[:, "Time"]
current_udds = udds_data[:, " Current"]

# 创建电流插值函数，用于获取任意时间点的电流值
current_interp = LinearInterpolation(time_udds, current_udds, extrapolation_bc=Flat())

time = collect(Float64, 0:0.01:maximum(time_udds))

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.5

# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = time  # 使用与电流数据相同的时间步长

# 设置电流函数
opt.Current = t -> current_interp(t)

# 创建和运行P2D模型模拟
println("开始模拟UDDS电流工况... P2D模型")
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

# 选择UDDS周期中的十个平均分隔点进行空间分布分析
# 将总时间均匀分为十份
max_time = maximum(time_result)
analysis_times = collect(range(max_time * 0.1, max_time, length=10))
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
println("开始模拟UDDS电流工况... sP2D模型")
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

# 定义科学论文风格的配色方案
# 修改颜色方案定义
# 将颜色定义从电极类型改为模型类型
p2d_color = :blueviolet  # P2D模型颜色 - 蓝紫色
sp2d_color = :orangered  # SP2D模型颜色 - 橙红色

# 透明度递减，越晚的时间点越透明
# 计算10个时间点的透明度，从0.95到0.35，时间越靠后越透明
time_alphas = collect(range(0.30, 0.95, length=10))

# 重新定义time_colors变量，使其包含P2D和SP2D模型的颜色及透明度
time_colors = Dict(
    "p2d_color" => p2d_color,
    "sp2d_color" => sp2d_color,
    "alphas" => time_alphas
)

# 为误差图保留原来的彩色方案
error_colors = [:dodgerblue, :royalblue, :steelblue, :midnightblue, :darkturquoise, 
            :crimson, :firebrick, :tomato, :coral, :brown]

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
         legendrows=1,
         legendcolumns=2,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的浓度分布
for i in 1:length(time_indices)
    current_alpha = time_colors["alphas"][i]  # 获取当前时间点的透明度
    
    # P2D模型用实线，颜色为蓝紫色
    # 正极浓度 - P2D
    plot!(x_grid_p_normalized, concentration_p_p2d_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # 负极浓度 - P2D
    plot!(x_grid_n_normalized, concentration_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # sP2D模型用五角星标记 + 虚线连接，颜色为橙红色
    # 正极浓度 - sP2D
    plot!(x_grid_p_normalized, concentration_p_sP2D_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_p_normalized, concentration_p_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
    
    # 负极浓度 - sP2D
    plot!(x_grid_n_normalized, concentration_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_n_normalized, concentration_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
end

# 添加模型类型图例
plot!([], [], label="P2D (Solid Line)", lw=2.0, color=time_colors["p2d_color"], linestyle=:solid)
# 对于sP2D，绘制线条和五角星在同一个图例中
plot!([], [], label="sP2D (Dot Line + Star)", lw=2.0, color=time_colors["sp2d_color"], 
      linestyle=:dash, marker=:star5, markersize=4, markerstrokewidth=0)

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域填充
separator_mid = (negative_end + positive_start) / 2
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 浓度图专用
neg_center_conc = negative_end / 2
pos_center_conc = (positive_start + positive_end) / 2
annotate!(neg_center_conc, maximum([maximum(conc) for conc in concentration_p2d_spatial])*0.70, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_conc, maximum([maximum(conc) for conc in concentration_p_p2d_spatial])*1.25, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(conc) for conc in concentration_p2d_spatial])*0.55, 
         text("Separator", 9, :center, :italic, :gray))

# 创建应力空间分布图
p3 = plot(xlabel="Normalized Position", 
         ylabel="σ/E (dimensionless)", 
         title="Electrode Stress Distribution - Normalized",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的应力分布
for i in 1:length(time_indices)
    current_alpha = time_colors["alphas"][i]  # 获取当前时间点的透明度
    
    # P2D模型用实线，颜色为蓝紫色
    # 正极应力 - P2D
    plot!(x_grid_p_normalized, stress_p_p2d_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # 负极应力 - P2D
    plot!(x_grid_n_normalized, stress_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # sP2D模型用五角星标记 + 虚线连接，颜色为橙红色
    # 正极应力 - sP2D
    plot!(x_grid_p_normalized, stress_p_sP2D_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_p_normalized, stress_p_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
    
    # 负极应力 - sP2D
    plot!(x_grid_n_normalized, stress_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_n_normalized, stress_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
end

# 添加模型类型图例
plot!([], [], label="P2D (Solid Line)", lw=2.0, color=time_colors["p2d_color"], linestyle=:solid)
# 对于sP2D，绘制线条和五角星在同一个图例中
plot!([], [], label="sP2D (Dash Line + Star)", lw=2.0, color=time_colors["sp2d_color"], 
      linestyle=:dash, marker=:star5, markersize=4, markerstrokewidth=0)

# 添加分隔线和区域标识  
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 应力图专用
neg_center_stress = negative_end / 2
pos_center_stress = (positive_start + positive_end) / 2
annotate!(neg_center_stress, maximum([maximum(stress) for stress in stress_p2d_spatial])-0.00006, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_stress, maximum([maximum(stress) for stress in stress_p_p2d_spatial])+ 0.00025, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(stress) for stress in stress_p2d_spatial])*0.54, 
         text("Separator", 9, :center, :italic, :gray))

# 创建电流图
p4 = plot(time_udds, current_udds, 
    label="Current Profile", 
    xlabel="Time [s]", 
    ylabel="Current [A]", 
    title="UDDS Current Profile",
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
         title="Electrode Reaction Current Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的反应电流分布
for i in 1:length(time_indices)
    current_alpha = time_colors["alphas"][i]  # 获取当前时间点的透明度
    
    # P2D模型用实线，颜色为蓝紫色
    # 正极反应电流 - P2D
    plot!(x_grid_p_normalized, reaction_current_p_p2d_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # 负极反应电流 - P2D
    plot!(x_grid_n_normalized, reaction_current_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # sP2D模型用五角星标记 + 虚线连接，颜色为橙红色
    # 正极反应电流 - sP2D
    plot!(x_grid_p_normalized, reaction_current_p_sP2D_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_p_normalized, reaction_current_p_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
    
    # 负极反应电流 - sP2D
    plot!(x_grid_n_normalized, reaction_current_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_n_normalized, reaction_current_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
end

# 添加模型类型图例
plot!([], [], label="P2D (Solid Line)", lw=2.0, color=time_colors["p2d_color"], linestyle=:solid)
# 对于sP2D，绘制线条和五角星在同一个图例中
plot!([], [], label="sP2D (Dash Line + Star)", lw=2.0, color=time_colors["sp2d_color"], 
      linestyle=:dash, marker=:star5, markersize=4, markerstrokewidth=0)

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="Separator")

# 区域标识 - 反应电流图专用
neg_center_reaction = negative_end / 2
pos_center_reaction = (positive_start + positive_end) / 2
annotate!(neg_center_reaction, maximum([maximum(curr) for curr in reaction_current_p2d_spatial])*0.80, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_reaction, maximum([maximum(curr) for curr in reaction_current_p_p2d_spatial]) + 0.7, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(curr) for curr in reaction_current_p2d_spatial])- 3.3, 
         text("Separator", 9, :center, :italic, :gray))

# 创建过电压空间分布图
p8 = plot(xlabel="Normalized Position", 
         ylabel="Overpotential [V]", 
         title="Electrode Overpotential Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=false,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500),dpi = 600)

# 绘制各时间点的过电压分布
for i in 1:length(time_indices)
    current_alpha = time_colors["alphas"][i]  # 获取当前时间点的透明度
    
    # P2D模型用实线，颜色为蓝紫色
    # 正极过电压 - P2D
    plot!(x_grid_p_normalized, overpotential_p_p2d_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # 负极过电压 - P2D
    plot!(x_grid_n_normalized, overpotential_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # sP2D模型用五角星标记 + 点线连接，颜色为橙红色
    # 正极过电压 - sP2D
    plot!(x_grid_p_normalized, overpotential_p_sP2D_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_p_normalized, overpotential_p_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
    
    # 负极过电压 - sP2D
    plot!(x_grid_n_normalized, overpotential_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dot, alpha=current_alpha)
    scatter!(x_grid_n_normalized, overpotential_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
end

# 添加模型类型图例
plot!([], [], label="P2D (Solid Line)", lw=2.0, color=time_colors["p2d_color"], linestyle=:solid)
# 对于sP2D，绘制线条和五角星在同一个图例中
plot!([], [], label="sP2D (Dash Line + Star)", lw=2.0, color=time_colors["sp2d_color"], 
      linestyle=:dash, marker=:star5, markersize=4, markerstrokewidth=0)

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 过电压图专用
neg_center_over = negative_end / 2
pos_center_over = (positive_start + positive_end) / 2
annotate!(neg_center_over, maximum([maximum(over) for over in overpotential_p2d_spatial])-0.220, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_over, maximum([maximum(over) for over in overpotential_p_p2d_spatial])+0.07, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(over) for over in overpotential_p2d_spatial])- 0.20, 
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
    current_alpha = time_colors["alphas"][i]  # 获取当前时间点的透明度
    
    # P2D模型用实线
    # 正极电势 - P2D (实线) - 使用橙红色
    plot!(x_grid_p_normalized, potential_p_p2d_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # 负极电势 - P2D (实线) - 使用蓝紫色
    plot!(x_grid_n_normalized, potential_p2d_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["p2d_color"], linestyle=:solid, alpha=current_alpha)
    
    # sP2D模型用五角星标记 + 虚线连接
    # 正极电势 - sP2D (五角星 + 虚线) - 使用橙红色
    plot!(x_grid_p_normalized, potential_p_sP2D_spatial[i], 
        label="", # 不在这里添加图例
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dash, alpha=current_alpha)
    scatter!(x_grid_p_normalized, potential_p_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
    
    # 负极电势 - sP2D (五角星 + 虚线) - 使用蓝紫色
    plot!(x_grid_n_normalized, potential_sP2D_spatial[i], 
        label="", 
        lw=2.0, color=time_colors["sp2d_color"], linestyle=:dash, alpha=current_alpha)
    scatter!(x_grid_n_normalized, potential_sP2D_spatial[i], 
        label="", 
        markersize=4, color=time_colors["sp2d_color"], marker=:star5, alpha=current_alpha)
end

# 添加时间点图例 - 为所有时间点添加透明度图例
for i in 1:length(time_indices)
    current_alpha = time_colors["alphas"][i]
    # 在图表的右侧添加时间点图例
    plot!([], [], 
        label="t = $(round(analysis_times[i], digits=1))s", 
        lw=2.0, color=:black, linestyle=:solid, alpha=current_alpha,
        legend=:outerright)
end

# 添加模型和电极类型图例
plot!([], [], label="正极 P2D (Solid Line)", lw=2.0, color=time_colors["p2d_color"], linestyle=:solid)
plot!([], [], label="负极 P2D (Solid Line)", lw=2.0, color=time_colors["p2d_color"], linestyle=:solid)
plot!([], [], label="正极 sP2D (Dash Line + Star)", lw=2.0, color=time_colors["sp2d_color"], linestyle=:dash)
scatter!([], [], label="", markersize=4, color=time_colors["sp2d_color"], marker=:star5)
plot!([], [], label="负极 sP2D (Dash Line + Star)", lw=2.0, color=time_colors["sp2d_color"], linestyle=:dash)
scatter!([], [], label="", markersize=4, color=time_colors["sp2d_color"], marker=:star5)

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
savefig(p1, "udds_voltage_spatial_comparison.pdf")
savefig(p2, "udds_concentration_spatial_comparison.png")
savefig(p3, "udds_stress_spatial_comparison.png")
savefig(p4, "udds_current_profile_spatial.pdf")
savefig(p7, "udds_reaction_current_spatial_comparison.png")
savefig(p8, "udds_overpotential_spatial_comparison.png")
savefig(p9, "udds_potential_spatial_comparison.png")
savefig(plot_combined, "udds_spatial_analysis_comparison_super.pdf")

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
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, concentration_spatial_errors_p[i], 
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加时间点图例 - 为所有时间点添加颜色图例
for i in 1:length(time_indices)
    # 在图表的右侧添加时间点图例
    plot!([], [], 
        label="t = $(round(analysis_times[i], digits=1))s", 
        lw=2.0, color=error_colors[i], linestyle=:solid, alpha=0.8,
        legend=:outerright)
end

# 添加电极类型图例
plot!([], [], label="Negative Electrode (Solid Line)", lw=2.0, color=:black, linestyle=:solid)
plot!([], [], label="Positive Electrode (Dash Line)", lw=2.0, color=:black, linestyle=:dash)

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 浓度误差图专用
neg_center_conc_err = negative_end / 2
pos_center_conc_err = (positive_start + positive_end) / 2
annotate!(neg_center_conc_err, maximum([maximum(err) for err in concentration_spatial_errors_n])*1.3, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_conc_err, maximum([maximum(err) for err in concentration_spatial_errors_p])*0.90, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in concentration_spatial_errors_n])*3.5, 
         text("Separator", 9, :center, :italic, :gray))

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
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, stress_spatial_errors_p[i], 
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加时间点图例 - 为所有时间点添加颜色图例
for i in 1:length(time_indices)
    # 在图表的右侧添加时间点图例
    plot!([], [], 
        label="t = $(round(analysis_times[i], digits=1))s", 
        lw=2.0, color=error_colors[i], linestyle=:solid, alpha=0.8,
        legend=:outerright)
end

# 添加电极类型图例
plot!([], [], label="Negative Electrode (Solid Line)", lw=2.0, color=:black, linestyle=:solid)
plot!([], [], label="Positive Electrode (Dash Line)", lw=2.0, color=:black, linestyle=:dash)

# 添加分隔线和区域标识
# 添加隔膜区域的边界线
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
# 添加隔膜区域标识
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 应力误差图专用
neg_center_stress_err = negative_end / 2
pos_center_stress_err = (positive_start + positive_end) / 2
annotate!(neg_center_stress_err, maximum([maximum(err) for err in stress_spatial_errors_n])*4.0, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_stress_err, maximum([maximum(err) for err in stress_spatial_errors_p])*0.7, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in stress_spatial_errors_n])*3.0, 
         text("Separator", 8, :center, :italic, :gray))

savefig(p5, "udds_concentration_spatial_relative_error.png")
savefig(p6, "udds_stress_spatial_relative_error.png")

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
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, reaction_current_spatial_errors_p[i], 
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加时间点图例 - 为所有时间点添加颜色图例
for i in 1:length(time_indices)
    # 在图表的右侧添加时间点图例
    plot!([], [], 
        label="t = $(round(analysis_times[i], digits=1))s", 
        lw=2.0, color=error_colors[i], linestyle=:solid, alpha=0.8,
        legend=:outerright)
end

# 添加电极类型图例
plot!([], [], label="Negative Electrode (Solid Line)", lw=2.0, color=:black, linestyle=:solid)
plot!([], [], label="Positive Electrode (Dash Line)", lw=2.0, color=:black, linestyle=:dash)

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 反应电流误差图专用
neg_center_reaction_err = negative_end / 2
pos_center_reaction_err = (positive_start + positive_end) / 2
annotate!(neg_center_reaction_err, maximum([maximum(err) for err in reaction_current_spatial_errors_n])+400, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_reaction_err, maximum([maximum(err) for err in reaction_current_spatial_errors_p])*0.85, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in reaction_current_spatial_errors_n])+200, 
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
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, overpotential_spatial_errors_p[i], 
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加时间点图例 - 为所有时间点添加颜色图例
for i in 1:length(time_indices)
    # 在图表的右侧添加时间点图例
    plot!([], [], 
        label="t = $(round(analysis_times[i], digits=1))s", 
        lw=2.0, color=error_colors[i], linestyle=:solid, alpha=0.8,
        legend=:outerright)
end

# 添加电极类型图例
plot!([], [], label="Negative Electrode (Solid Line)", lw=2.0, color=:black, linestyle=:solid)
plot!([], [], label="Positive Electrode (Dash Line)", lw=2.0, color=:black, linestyle=:dash)

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识 - 过电压误差图专用
neg_center_over_err = negative_end / 2
pos_center_over_err = (positive_start + positive_end) / 2
annotate!(neg_center_over_err, maximum([maximum(err) for err in overpotential_spatial_errors_n])*10.0, 
         text("Negative Electrode", 9, :center, :bold))
annotate!(pos_center_over_err, maximum([maximum(err) for err in overpotential_spatial_errors_p])*0.7, 
         text("Positive Electrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(err) for err in overpotential_spatial_errors_n])*20.0, 
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

# 只为第1、第5和第10个时间点显示图例
key_indices = [1, 5, 10]

# 绘制负极误差 - 使用实线
for i in 1:length(time_indices)
    plot!(x_grid_n_normalized, potential_spatial_errors_n[i], 
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:solid,
        alpha=0.8)
end

# 绘制正极误差 - 使用虚线
for i in 1:length(time_indices)
    plot!(x_grid_p_normalized, potential_spatial_errors_p[i], 
        label="", # 不在这里添加图例
        lw=2.0, 
        color=error_colors[i],
        linestyle=:dash,
        alpha=0.8)
end

# 添加时间点图例 - 为所有时间点添加颜色图例
for i in 1:length(time_indices)
    # 在图表的右侧添加时间点图例
    plot!([], [], 
        label="t = $(round(analysis_times[i], digits=1))s", 
        lw=2.0, color=error_colors[i], linestyle=:solid, alpha=0.8,
        legend=:outerright)
end

# 添加电极类型图例
plot!([], [], label="Negative Electrode (Solid Line)", lw=2.0, color=:black, linestyle=:solid)
plot!([], [], label="Positive Electrode (Dash Line)", lw=2.0, color=:black, linestyle=:dash)

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

savefig(p10, "udds_reaction_current_spatial_relative_error.png")
savefig(p11, "udds_overpotential_spatial_relative_error.png")
savefig(p12, "udds_potential_spatial_relative_error.png")

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

# 4. UDDS工况下的特殊影响
println("\n4. UDDS工况下的特殊影响")
println("-"^50)

println("UDDS对sP2D误差的影响:")
println("  A. 快速变化挑战:")
println("     - 脉冲电流快速变化 → 电势响应滞后")
println("     - sP2D简化模型响应速度可能跟不上P2D")
println("     - 边界层厚度动态变化，固定拟合参数不够精确")
println()
println("  B. 空间-时间耦合:")
println("     - UDDS产生复杂的空间-时间耦合效应")
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
println("\n=== UDDS电流工况空间分布分析结果 ===")
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
CSV.write("udds_spatial_analysis_negative_electrode.csv", negative_df)
CSV.write("udds_spatial_analysis_positive_electrode.csv", positive_df)

println("UDDS电流空间分布分析完成，结果已保存")

println("\n=== 分析时间点的电流值 ===")
for t in analysis_times
    current_value = current_interp(t)
    println("时间 $(t)s 的电流值: $(round(current_value, digits=4)) A")
end

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

# ===== 计算各空间点在所有时间线下的应力RMSE =====
println("\n" * "="^60)
println("空间点应力RMSE分析 - 寻找危险点")
println("="^60)

# 计算空间点应力RMSE和RRMSE的函数
function calculate_spatial_stress_rmse()
    # 计算负极各空间点在所有时间点下的应力RMSE和RRMSE
    n_positions = length(x_grid_n)
    negative_spatial_rmse = zeros(n_positions)
    negative_spatial_rrmse = zeros(n_positions)  # 相对均方根误差
    
    # 遍历负极每个空间点
    for pos_idx in 1:n_positions
        stress_errors = Float64[]
        stress_relative_errors = Float64[]
        
        # 遍历所有时间点
        for time_idx in 1:length(time_indices)
            # 提取该空间点在该时间点的P2D和sP2D应力值
            p2d_stress = stress_p2d_spatial[time_idx][pos_idx]
            sp2d_stress = stress_sP2D_spatial[time_idx][pos_idx]
            
            # 计算该时间点的应力误差(绝对误差的平方)
            stress_error = (sp2d_stress - p2d_stress)^2
            push!(stress_errors, stress_error)
            
            # 计算该时间点的应力相对误差(相对误差的平方)
            # 避免除以零的问题
            if abs(p2d_stress) > 1e-10
                stress_relative_error = ((sp2d_stress - p2d_stress) / p2d_stress)^2
            else
                stress_relative_error = 0.0
            end
            push!(stress_relative_errors, stress_relative_error)
        end
        
        # 计算该空间点在所有时间点下的RMSE
        if !isempty(stress_errors)
            negative_spatial_rmse[pos_idx] = sqrt(mean(stress_errors))
        end
        
        # 计算该空间点在所有时间点下的RRMSE (相对均方根误差)
        if !isempty(stress_relative_errors)
            negative_spatial_rrmse[pos_idx] = sqrt(mean(stress_relative_errors))
        end
    end
    
    # 计算正极各空间点在所有时间点下的应力RMSE和RRMSE
    p_positions = length(x_grid_p)
    positive_spatial_rmse = zeros(p_positions)
    positive_spatial_rrmse = zeros(p_positions)  # 相对均方根误差
    
    # 遍历正极每个空间点
    for pos_idx in 1:p_positions
        stress_errors = Float64[]
        stress_relative_errors = Float64[]
        
        # 遍历所有时间点
        for time_idx in 1:length(time_indices)
            # 提取该空间点在该时间点的P2D和sP2D应力值
            p2d_stress = stress_p_p2d_spatial[time_idx][pos_idx]
            sp2d_stress = stress_p_sP2D_spatial[time_idx][pos_idx]
            
            # 计算该时间点的应力误差
            stress_error = (sp2d_stress - p2d_stress)^2
            push!(stress_errors, stress_error)
            
            # 计算该时间点的应力相对误差
            if abs(p2d_stress) > 1e-10
                stress_relative_error = ((sp2d_stress - p2d_stress) / p2d_stress)^2
            else
                stress_relative_error = 0.0
            end
            push!(stress_relative_errors, stress_relative_error)
        end
        
        # 计算该空间点在所有时间点下的RMSE
        if !isempty(stress_errors)
            positive_spatial_rmse[pos_idx] = sqrt(mean(stress_errors))
        end
        
        # 计算该空间点在所有时间点下的RRMSE
        if !isempty(stress_relative_errors)
            positive_spatial_rrmse[pos_idx] = sqrt(mean(stress_relative_errors))
        end
    end
    
    return negative_spatial_rmse, positive_spatial_rmse, negative_spatial_rrmse, positive_spatial_rrmse
end

# 计算空间点应力RMSE和RRMSE
negative_spatial_stress_rmse, positive_spatial_stress_rmse, negative_spatial_stress_rrmse, positive_spatial_stress_rrmse = calculate_spatial_stress_rmse()

# === 绝对均方根误差(RMSE)分析 ===
# 找出负极应力RMSE最大的点
max_n_rmse, max_n_idx = findmax(negative_spatial_stress_rmse)
max_n_position = x_grid_n[max_n_idx]
max_n_position_normalized = x_grid_n_normalized[max_n_idx]

# 找出正极应力RMSE最大的点
max_p_rmse, max_p_idx = findmax(positive_spatial_stress_rmse)
max_p_position = x_grid_p[max_p_idx]
max_p_position_normalized = x_grid_p_normalized[max_p_idx]

# 计算最大应力RMSE点的相对位置
n_relative_position = max_n_position_normalized / negative_end
p_relative_position = (max_p_position_normalized - positive_start) / (1.0 - positive_start)

# === 相对均方根误差(RRMSE)分析 ===
# 找出负极应力RRMSE最大的点
max_n_rrmse, max_n_rrmse_idx = findmax(negative_spatial_stress_rrmse)
max_n_rrmse_position = x_grid_n[max_n_rrmse_idx]
max_n_rrmse_position_normalized = x_grid_n_normalized[max_n_rrmse_idx]

# 找出正极应力RRMSE最大的点
max_p_rrmse, max_p_rrmse_idx = findmax(positive_spatial_stress_rrmse)
max_p_rrmse_position = x_grid_p[max_p_rrmse_idx]
max_p_rrmse_position_normalized = x_grid_p_normalized[max_p_rrmse_idx]

# 计算最大应力RRMSE点的相对位置
n_rrmse_relative_position = max_n_rrmse_position_normalized / negative_end
p_rrmse_relative_position = (max_p_rrmse_position_normalized - positive_start) / (1.0 - positive_start)

# 输出分析结果
println("\n=== 绝对均方根误差(RMSE)分析 ===")
println("\n负极应力RMSE分析:")
println("  最大应力RMSE: $(round(max_n_rmse, digits=6)) Pa")
println("  最大应力RMSE位置: $(round(max_n_position*1e6, digits=2)) μm (绝对位置)")
println("  归一化位置: $(round(max_n_position_normalized, digits=4)) (归一化位置)")
println("  相对电极位置: $(round(n_relative_position*100, digits=1))% (从集流体算起)")

println("\n正极应力RMSE分析:")
println("  最大应力RMSE: $(round(max_p_rmse, digits=6)) Pa")
println("  最大应力RMSE位置: $(round(max_p_position*1e6, digits=2)) μm (绝对位置)")
println("  归一化位置: $(round(max_p_position_normalized, digits=4)) (归一化位置)")
println("  相对电极位置: $(round(p_relative_position*100, digits=1))% (从集流体算起)")

# 比较负极和正极哪个是绝对误差的整体危险点
if max_n_rmse > max_p_rmse
    println("\nRMSE分析整体危险点位于负极，RMSE比正极高 $(round(max_n_rmse/max_p_rmse, digits=2)) 倍")
    overall_dangerous_position_rmse = max_n_position
    overall_dangerous_normalized_rmse = max_n_position_normalized
    println("RMSE整体危险点绝对位置: $(round(overall_dangerous_position_rmse*1e6, digits=2)) μm")
    println("RMSE整体危险点归一化位置: $(round(overall_dangerous_normalized_rmse, digits=4))")
else
    println("\nRMSE分析整体危险点位于正极，RMSE比负极高 $(round(max_p_rmse/max_n_rmse, digits=2)) 倍")
    overall_dangerous_position_rmse = max_p_position
    overall_dangerous_normalized_rmse = max_p_position_normalized
    println("RMSE整体危险点绝对位置: $(round(overall_dangerous_position_rmse*1e6, digits=2)) μm")
    println("RMSE整体危险点归一化位置: $(round(overall_dangerous_normalized_rmse, digits=4))")
end

println("\n=== 相对均方根误差(RRMSE)分析 ===")
println("\n负极应力RRMSE分析:")
println("  最大应力RRMSE: $(round(max_n_rrmse, digits=6)) (无量纲)")
println("  最大应力RRMSE位置: $(round(max_n_rrmse_position*1e6, digits=2)) μm (绝对位置)")
println("  归一化位置: $(round(max_n_rrmse_position_normalized, digits=4)) (归一化位置)")
println("  相对电极位置: $(round(n_rrmse_relative_position*100, digits=1))% (从集流体算起)")

println("\n正极应力RRMSE分析:")
println("  最大应力RRMSE: $(round(max_p_rrmse, digits=6)) (无量纲)")
println("  最大应力RRMSE位置: $(round(max_p_rrmse_position*1e6, digits=2)) μm (绝对位置)")
println("  归一化位置: $(round(max_p_rrmse_position_normalized, digits=4)) (归一化位置)")
println("  相对电极位置: $(round(p_rrmse_relative_position*100, digits=1))% (从集流体算起)")

# 比较负极和正极哪个是相对误差的整体危险点
if max_n_rrmse > max_p_rrmse
    println("\nRRMSE分析整体危险点位于负极，RRMSE比正极高 $(round(max_n_rrmse/max_p_rrmse, digits=2)) 倍")
    overall_dangerous_position_rrmse = max_n_rrmse_position
    overall_dangerous_normalized_rrmse = max_n_rrmse_position_normalized
    println("RRMSE整体危险点绝对位置: $(round(overall_dangerous_position_rrmse*1e6, digits=2)) μm")
    println("RRMSE整体危险点归一化位置: $(round(overall_dangerous_normalized_rrmse, digits=4))")
else
    println("\nRRMSE分析整体危险点位于正极，RRMSE比负极高 $(round(max_p_rrmse/max_n_rrmse, digits=2)) 倍")
    overall_dangerous_position_rrmse = max_p_rrmse_position
    overall_dangerous_normalized_rrmse = max_p_rrmse_position_normalized
    println("RRMSE整体危险点绝对位置: $(round(overall_dangerous_position_rrmse*1e6, digits=2)) μm")
    println("RRMSE整体危险点归一化位置: $(round(overall_dangerous_normalized_rrmse, digits=4))")
end

# 确定最终采用的危险点指标（这里可以根据需要选择RMSE或RRMSE作为最终判断标准）
println("\n=== 危险点综合判断 ===")
println("基于RMSE和RRMSE的综合分析，推荐采用RRMSE作为危险点判断标准，因为它能消除量纲影响，更好地反映相对误差水平。")

if max_n_rrmse > max_p_rrmse
    overall_dangerous_position = max_n_rrmse_position
    overall_dangerous_normalized = max_n_rrmse_position_normalized
    println("最终确定的整体危险点位于负极")
else
    overall_dangerous_position = max_p_rrmse_position
    overall_dangerous_normalized = max_p_rrmse_position_normalized
    println("最终确定的整体危险点位于正极")
end
println("整体危险点绝对位置: $(round(overall_dangerous_position*1e6, digits=2)) μm")
println("整体危险点归一化位置: $(round(overall_dangerous_normalized, digits=4))")

# 绘制空间点应力RMSE分布图
p_stress_rmse = plot(xlabel="Normalized Position", 
                    ylabel="Stress RMSE [Pa]", 
                    title="Spatial Stress RMSE Distribution",
                    legend=:outerright,
                    xlims=(0,1),
                    grid=false,
                    gridalpha=0.3,
                    framestyle=:box,
                    size=(800, 500),dpi = 600)

# 绘制负极RMSE分布
plot!(x_grid_n_normalized, negative_spatial_stress_rmse, 
    label="Negative Electrode", 
    lw=2.5, color=:royalblue, linestyle=:solid)

# 绘制正极RMSE分布
plot!(x_grid_p_normalized, positive_spatial_stress_rmse, 
    label="Positive Electrode", 
    lw=2.5, color=:firebrick, linestyle=:solid)

# 标记负极最大RMSE点
scatter!([max_n_position_normalized], [max_n_rmse], 
    label="Negative Max RMSE", 
    markersize=8, color=:blue, marker=:star8)

# 标记正极最大RMSE点
scatter!([max_p_position_normalized], [max_p_rmse], 
    label="Positive Max RMSE", 
    markersize=8, color=:red, marker=:star8)

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识
neg_center_stress_rmse = negative_end / 2
pos_center_stress_rmse = (positive_start + positive_end) / 2
annotate!(neg_center_stress_rmse, maximum(negative_spatial_stress_rmse)*0.8, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_stress_rmse, maximum(positive_spatial_stress_rmse)*0.8, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(negative_spatial_stress_rmse), maximum(positive_spatial_stress_rmse)])*0.5, 
         text("Separator", 9, :center, :italic, :gray))

# 为最大RMSE点添加文本标记
if max_n_rmse > max_p_rmse
    annotate!(max_n_position_normalized, max_n_rmse*1.15, 
         text("RMSE Dangerous Point", 8, :center, :bold, :blue))
else
    annotate!(max_p_position_normalized, max_p_rmse*1.15, 
         text("RMSE Dangerous Point", 8, :center, :bold, :red))
end

savefig(p_stress_rmse, "spatial_stress_rmse_distribution.png")

# 绘制空间点应力相对RMSE分布图
p_stress_rrmse = plot(xlabel="Normalized Position", 
                    ylabel="Stress RRMSE (dimensionless)", 
                    title="Spatial Stress Relative RMSE Distribution",
                    legend=:outerright,
                    xlims=(0,1),
                    grid=false,
                    gridalpha=0.3,
                    framestyle=:box,
                    size=(800, 500),dpi = 600)

# 绘制负极RRMSE分布
plot!(x_grid_n_normalized, negative_spatial_stress_rrmse, 
    label="Negative Electrode", 
    lw=2.5, color=:royalblue, linestyle=:solid)

# 绘制正极RRMSE分布
plot!(x_grid_p_normalized, positive_spatial_stress_rrmse, 
    label="Positive Electrode", 
    lw=2.5, color=:firebrick, linestyle=:solid)

# 标记负极最大RRMSE点
scatter!([max_n_rrmse_position_normalized], [max_n_rrmse], 
    label="Negative Max RRMSE", 
    markersize=8, color=:blue, marker=:star8)

# 标记正极最大RRMSE点
scatter!([max_p_rrmse_position_normalized], [max_p_rrmse], 
    label="Positive Max RRMSE", 
    markersize=8, color=:red, marker=:star8)

# 添加分隔线和区域标识
vline!([negative_end], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([positive_start], color=:gray, linestyle=:dash, alpha=0.5, linewidth=1.5, label="")
vline!([separator_mid], color=:lightgray, linestyle=:dot, alpha=0.3, linewidth=3, label="")

# 区域标识
neg_center_stress_rrmse = negative_end / 2
pos_center_stress_rrmse = (positive_start + positive_end) / 2
annotate!(neg_center_stress_rrmse, maximum(negative_spatial_stress_rrmse)*0.8, 
         text("Negative\nElectrode", 9, :center, :bold))
annotate!(pos_center_stress_rrmse, maximum(positive_spatial_stress_rrmse)*0.8, 
         text("Positive\nElectrode", 9, :center, :bold))
annotate!(separator_mid, maximum([maximum(negative_spatial_stress_rrmse), maximum(positive_spatial_stress_rrmse)])*0.5, 
         text("Separator", 9, :center, :italic, :gray))

# 为最大RRMSE点添加文本标记
if max_n_rrmse > max_p_rrmse
    annotate!(max_n_rrmse_position_normalized, max_n_rrmse*1.15, 
         text("RRMSE Dangerous Point", 8, :center, :bold, :blue))
else
    annotate!(max_p_rrmse_position_normalized, max_p_rrmse*1.15, 
         text("RRMSE Dangerous Point", 8, :center, :bold, :red))
end

savefig(p_stress_rrmse, "spatial_stress_rrmse_distribution.png")

# 创建组合图：RMSE和RRMSE对比
p_stress_comparison = plot(
    p_stress_rmse, p_stress_rrmse,
    layout = (2, 1),
    size = (800, 1000),
    dpi = 600
)
savefig(p_stress_comparison, "spatial_stress_rmse_rrmse_comparison.png")

# 进一步分析最危险点
println("\n" * "="^60)
println("危险点详细分析")
println("="^60)

# === 分析RMSE危险点 ===
println("\n==== 基于RMSE的危险点分析 ====")

# 确定RMSE最危险点所在的电极和索引
dangerous_electrode_rmse = max_n_rmse > max_p_rmse ? "negative" : "positive"
dangerous_idx_rmse = max_n_rmse > max_p_rmse ? max_n_idx : max_p_idx

# 提取RMSE危险点在各时间点的应力值
dangerous_p2d_stress_rmse = []
dangerous_sp2d_stress_rmse = []
dangerous_stress_error_rmse = []
dangerous_stress_rel_error_rmse = []

if dangerous_electrode_rmse == "negative"
    for time_idx in 1:length(time_indices)
        p2d_val = stress_p2d_spatial[time_idx][dangerous_idx_rmse]
        sp2d_val = stress_sP2D_spatial[time_idx][dangerous_idx_rmse]
        push!(dangerous_p2d_stress_rmse, p2d_val)
        push!(dangerous_sp2d_stress_rmse, sp2d_val)
        error = abs(sp2d_val - p2d_val)
        rel_error = abs(p2d_val) > 1e-10 ? abs((sp2d_val - p2d_val) / p2d_val) : 0.0
        push!(dangerous_stress_error_rmse, error)
        push!(dangerous_stress_rel_error_rmse, rel_error)
    end
else
    for time_idx in 1:length(time_indices)
        p2d_val = stress_p_p2d_spatial[time_idx][dangerous_idx_rmse]
        sp2d_val = stress_p_sP2D_spatial[time_idx][dangerous_idx_rmse]
        push!(dangerous_p2d_stress_rmse, p2d_val)
        push!(dangerous_sp2d_stress_rmse, sp2d_val)
        error = abs(sp2d_val - p2d_val)
        rel_error = abs(p2d_val) > 1e-10 ? abs((sp2d_val - p2d_val) / p2d_val) : 0.0
        push!(dangerous_stress_error_rmse, error)
        push!(dangerous_stress_rel_error_rmse, rel_error)
    end
end

# 找出RMSE危险点应力误差最大的时间点
max_error_value_rmse, max_error_idx_rmse = findmax(dangerous_stress_error_rmse)
critical_time_rmse = analysis_times[max_error_idx_rmse]
critical_current_rmse = current_interp(critical_time_rmse)

println("\nRMSE危险点随时间的应力变化:")
println("  危险点位于: $(dangerous_electrode_rmse == "negative" ? "负极" : "正极")")
println("  危险时刻: $(round(critical_time_rmse, digits=2)) s")
println("  危险时刻电流: $(round(critical_current_rmse, digits=4)) A")
println("  最大绝对应力误差: $(round(max_error_value_rmse, digits=6)) Pa")
println("  对应相对误差: $(round(dangerous_stress_rel_error_rmse[max_error_idx_rmse]*100, digits=2))%")

# === 分析RRMSE危险点 ===
println("\n==== 基于RRMSE的危险点分析 ====")

# 确定RRMSE最危险点所在的电极和索引
dangerous_electrode_rrmse = max_n_rrmse > max_p_rrmse ? "negative" : "positive"
dangerous_idx_rrmse = max_n_rrmse > max_p_rrmse ? max_n_rrmse_idx : max_p_rrmse_idx

# 提取RRMSE危险点在各时间点的应力值
dangerous_p2d_stress_rrmse = []
dangerous_sp2d_stress_rrmse = []
dangerous_stress_error_rrmse = []
dangerous_stress_rel_error_rrmse = []

if dangerous_electrode_rrmse == "negative"
    for time_idx in 1:length(time_indices)
        p2d_val = stress_p2d_spatial[time_idx][dangerous_idx_rrmse]
        sp2d_val = stress_sP2D_spatial[time_idx][dangerous_idx_rrmse]
        push!(dangerous_p2d_stress_rrmse, p2d_val)
        push!(dangerous_sp2d_stress_rrmse, sp2d_val)
        error = abs(sp2d_val - p2d_val)
        rel_error = abs(p2d_val) > 1e-10 ? abs((sp2d_val - p2d_val) / p2d_val) : 0.0
        push!(dangerous_stress_error_rrmse, error)
        push!(dangerous_stress_rel_error_rrmse, rel_error)
    end
else
    for time_idx in 1:length(time_indices)
        p2d_val = stress_p_p2d_spatial[time_idx][dangerous_idx_rrmse]
        sp2d_val = stress_p_sP2D_spatial[time_idx][dangerous_idx_rrmse]
        push!(dangerous_p2d_stress_rrmse, p2d_val)
        push!(dangerous_sp2d_stress_rrmse, sp2d_val)
        error = abs(sp2d_val - p2d_val)
        rel_error = abs(p2d_val) > 1e-10 ? abs((sp2d_val - p2d_val) / p2d_val) : 0.0
        push!(dangerous_stress_error_rrmse, error)
        push!(dangerous_stress_rel_error_rrmse, rel_error)
    end
end

# 找出RRMSE危险点相对应力误差最大的时间点
max_rel_error_value_rrmse, max_rel_error_idx_rrmse = findmax(dangerous_stress_rel_error_rrmse)
critical_time_rrmse = analysis_times[max_rel_error_idx_rrmse]
critical_current_rrmse = current_interp(critical_time_rrmse)

println("\nRRMSE危险点随时间的应力变化:")
println("  危险点位于: $(dangerous_electrode_rrmse == "negative" ? "负极" : "正极")")
println("  危险时刻: $(round(critical_time_rrmse, digits=2)) s")
println("  危险时刻电流: $(round(critical_current_rrmse, digits=4)) A")
println("  最大相对应力误差: $(round(max_rel_error_value_rrmse*100, digits=2))%")
println("  对应绝对误差: $(round(dangerous_stress_error_rrmse[max_rel_error_idx_rrmse], digits=6)) Pa")

# 选择用于进一步可视化的危险点（默认使用RRMSE）
dangerous_electrode = dangerous_electrode_rrmse
dangerous_idx = dangerous_idx_rrmse
dangerous_p2d_stress = dangerous_p2d_stress_rrmse
dangerous_sp2d_stress = dangerous_sp2d_stress_rrmse
dangerous_stress_error = dangerous_stress_error_rrmse
max_error_idx = max_rel_error_idx_rrmse
critical_time = critical_time_rrmse
critical_current = critical_current_rrmse

# 绘制危险点应力随时间的变化
p_dangerous = plot(xlabel="Time [s]", 
                   ylabel="Stress [Pa]", 
                   title="Dangerous Point Stress vs Time",
                   legend=:outerright,
                   grid=false,
                   gridalpha=0.3,
                   framestyle=:box,
                   size=(800, 500),dpi = 600)

# 绘制P2D和sP2D在危险点的应力
plot!(analysis_times, dangerous_p2d_stress, 
    label="P2D Model", 
    lw=2.5, color=time_colors["p2d_color"], linestyle=:solid)

plot!(analysis_times, dangerous_sp2d_stress, 
    label="sP2D Model", 
    lw=2.5, color=time_colors["sp2d_color"], linestyle=:solid)

# 标记最大误差点
scatter!([critical_time], [dangerous_p2d_stress[max_error_idx]], 
    label="Critical Point (P2D)", 
    markersize=8, color=time_colors["p2d_color"], marker=:star8)

scatter!([critical_time], [dangerous_sp2d_stress[max_error_idx]], 
    label="Critical Point (sP2D)", 
    markersize=8, color=time_colors["sp2d_color"], marker=:star8)

# 添加对应的电流变化到次坐标轴
p_dangerous_current = twinx(p_dangerous)
plot!(p_dangerous_current, analysis_times, [current_interp(t) for t in analysis_times], 
    label="Current", 
    lw=2.0, color=:darkgreen, linestyle=:dash, 
    ylabel="Current [A]")

# 标记临界电流点
scatter!(p_dangerous_current, [critical_time], [critical_current], 
    label="Critical Current", 
    markersize=6, color=:darkgreen, marker=:circle)

savefig(p_dangerous, "dangerous_point_stress_time_analysis.png")

# 分析危险点特征
println("\n危险点特征分析:")
println("  1. 位置特征: 位于$(dangerous_electrode == "negative" ? "负极" : "正极")内部")
println("  2. 应力特征: 最大RMSE为$(round(max_n_rmse > max_p_rmse ? max_n_rmse : max_p_rmse, digits=4)) Pa")
println("  3. 时间特征: 在$(round(critical_time, digits=2))s时达到最大误差")
println("  4. 电流特征: 临界电流为$(round(critical_current, digits=4)) A")

println("危险点形成机理:")
if dangerous_electrode == "negative"
    println("  1. 负极石墨材料在循环过程中容易形成锂浓度梯度")
    println("  2. 在高电流密度区域，锂离子浓度变化剧烈导致应力集中")
    println("  3. sP2D模型在此区域的简化处理导致较大预测误差")
else
    println("  1. 正极材料在高充放电率下应力梯度大")
    println("  2. 电极内部与界面处应力传递复杂，sP2D简化模型难以精确捕捉")
    println("  3. 电位分布与机械应力耦合效应在高电流工况下更明显")
end

println("\n实际应用建议:")
println("  1. 在电池设计时应特别关注已识别的危险点区域")
println("  2. 优化电极厚度分布和微结构以减轻应力集中")
println("  3. 动态工况下的电池管理系统应考虑危险点应力集中效应")
println("  4. 使用精细模型在关键区域进行更精确的预测")

println("\n" * "="^60)

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

# 打印时间尺度分析所需参数
println("--- 时间尺度分析所需参数 (LGM50) ---")
# 负极 (Negative Electrode)
println("负极 (Graphite):")
τ_diff_n = param_dim.NE.rs^2 / param_dim.NE.Ds
println("  - 粒子半径 (rs): ", param_dim.NE.rs, " m")
println("  - 固相扩散系数 (Ds): ", param_dim.NE.Ds, " m^2/s")
println("  - 固相电导率 (sigma): ", param_dim.NE.sig, " S/m")
println("  - 电极厚度 (thickness): ", param_dim.NE.thickness, " m")
println("  - 固相扩散时间尺度 (τ_diff = rs^2 / Ds): ", τ_diff_n, " s")

# 正极 (Positive Electrode)
println("\n正极 (NMC):")
τ_diff_p = param_dim.PE.rs^2 / param_dim.PE.Ds
println("  - 粒子半径 (rs): ", param_dim.PE.rs, " m")
println("  - 固相扩散系数 (Ds): ", param_dim.PE.Ds, " m^2/s")
println("  - 固相电导率 (sigma): ", param_dim.PE.sig, " S/m")
println("  - 电极厚度 (thickness): ", param_dim.PE.thickness, " m")
println("  - 固相扩散时间尺度 (τ_diff = rs^2 / Ds): ", τ_diff_p, " s")
println("------------------------------------------\n")