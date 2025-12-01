using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
if !isdefined(Main, :JuBat)
    include("../src/JuBat.jl")
end

# 设置数据路径
path = "D:/竞赛和课程文件/课程文件/毕业设计/SRC/data/drive_cycles/"

# 读取US06数据
us06_data = CSV.read(path * "US061.csv", DataFrame, header = 1)
time_us06 = us06_data[:, "time"]
current_us06 = us06_data[:, "current"]

# 创建电流插值函数，用于获取任意时间点的电流值
current_interp = LinearInterpolation(time_us06, current_us06, extrapolation_bc=Flat())

# 设置电池参数
param_dim = JuBat.ChooseCell("Enertech")
param_dim.cell.v_h = 4.3

# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = collect(Float64, 0:1:maximum(time_us06))  # 以1秒的间隔进行模拟

# 设置电流函数为US06表格中的数据
opt.Current = t -> current_interp(t)

# 创建和运行P2D模型模拟
println("开始模拟US06工况... P2D模型")
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
# 负极坐标归一化到 [0, 0.5]
x_grid_n_normalized = x_grid_n ./ maximum(x_grid_n) .* 0.5
# 正极坐标归一化到 [0.5, 1.0]
x_grid_p_normalized = 0.5 .+ (x_grid_p .- minimum(x_grid_p)) ./ (maximum(x_grid_p) - minimum(x_grid_p)) .* 0.5

# 选择US06周期中的关键时间点进行空间分布分析
# US06周期大约600秒，选择10个关键时间点
max_time = maximum(time_us06)
analysis_times = Int.(round.(collect(range(max_time/10, max_time, length=10))))
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

# 设置sP2D模型
opt.model = "sP2D"  # 使用sP2D模型
println("开始模拟US06工况... sP2D模型")
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
    # P2D模型数据
    if ndims(concentration_p2d) > 1
        push!(concentration_p2d_spatial, concentration_p2d[:, t_idx])
    end
    
    if ndims(stress_p2d) > 1
        push!(stress_p2d_spatial, stress_p2d[:, t_idx])
    end
    
    if ndims(concentration_p_p2d) > 1
        push!(concentration_p_p2d_spatial, concentration_p_p2d[:, t_idx])
    end
    
    if ndims(stress_p_p2d) > 1
        push!(stress_p_p2d_spatial, stress_p_p2d[:, t_idx])
    end
    
    # sP2D模型数据
    if ndims(concentration_sP2D) > 1
        push!(concentration_sP2D_spatial, concentration_sP2D[:, t_idx])
    end
    
    if ndims(stress_sP2D) > 1
        push!(stress_sP2D_spatial, stress_sP2D[:, t_idx])
    end
    
    if ndims(concentration_p_sP2D) > 1
        push!(concentration_p_sP2D_spatial, concentration_p_sP2D[:, t_idx])
    end
    
    if ndims(stress_p_sP2D) > 1
        push!(stress_p_sP2D_spatial, stress_p_sP2D[:, t_idx])
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
    size=(800, 400))

plot!(time_result, voltage_sP2D, 
    label="sP2D Model", 
    linestyle=:dash, color=:crimson, lw=2.5)

# 添加分析时间点的垂直线
for (i, t) in enumerate(analysis_times)
    vline!([t], color=:gray, alpha=0.4, linestyle=:dot, linewidth=1, label="")
end

# 创建浓度空间分布图
p2 = plot(xlabel="Normalized Position", 
         ylabel="Concentration [mol/m³]", 
         title="Electrode Concentration Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500))

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
vline!([0.5], color=:black, linestyle=:dash, alpha=0.7, linewidth=2, label="")
annotate!(0.25, maximum([maximum(conc) for conc in concentration_p2d_spatial])*0.95, 
         text("Negative Electrode", 10, :center, :bold))
annotate!(0.75, maximum([maximum(conc) for conc in concentration_p_p2d_spatial])*0.95, 
         text("Positive Electrode", 10, :center, :bold))

# 创建应力空间分布图
p3 = plot(xlabel="Normalized Position", 
         ylabel="Stress [Pa]", 
         title="Electrode Stress Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500))

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
vline!([0.5], color=:black, linestyle=:dash, alpha=0.7, linewidth=2, label="")
annotate!(0.25, maximum([maximum(stress) for stress in stress_p2d_spatial])*0.95, 
         text("Negative Electrode", 10, :center, :bold))
annotate!(0.75, maximum([maximum(stress) for stress in stress_p_p2d_spatial])*0.95, 
         text("Positive Electrode", 10, :center, :bold))

# 创建电流图
p4 = plot(time_result, [current_interp(t) for t in time_result], 
    label="Current Profile", 
    xlabel="Time [s]", 
    ylabel="Current [A]", 
    title="US06 Current Profile",
    lw=2.5, color=:darkorange,
    grid=true,
    gridalpha=0.3,
    framestyle=:box,
    size=(800, 400))

# 添加分析时间点的垂直线
for (i, t) in enumerate(analysis_times)
    vline!([t], color=:gray, alpha=0.4, linestyle=:dot, linewidth=1, label="")
end

# 组合图表
plot_combined = plot(p1, p2, p3, p4, 
                    layout=(4, 1), 
                    size=(1000, 1200))

# 保存图表
savefig(p1, "us06_voltage_comparison.pdf")
savefig(p2, "us06_concentration_spatial_comparison.png")
savefig(p3, "us06_stress_spatial_comparison.png")
savefig(p4, "us06_current_profile.pdf")
savefig(plot_combined, "us06_spatial_analysis_comparison_super.pdf")

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
         ylabel="Concentration Error [mol/m³]", 
         title="Spatial Concentration Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500))

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
vline!([0.5], color=:black, linestyle=:dash, alpha=0.7, linewidth=2, label="")
annotate!(0.25, maximum([maximum(err) for err in concentration_spatial_errors_n])*0.95, 
         text("Negative Electrode", 10, :center, :bold))
annotate!(0.75, maximum([maximum(err) for err in concentration_spatial_errors_p])*0.95, 
         text("Positive Electrode", 10, :center, :bold))

# 创建空间误差分布图 - 应力误差
p6 = plot(xlabel="Normalized Position", 
         ylabel="Stress Error [Pa]", 
         title="Spatial Stress Error Distribution",
         legend=:outerright,
         xlims=(0,1),
         grid=true,
         gridalpha=0.3,
         framestyle=:box,
         size=(800, 500))

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
vline!([0.5], color=:black, linestyle=:dash, alpha=0.7, linewidth=2, label="")
annotate!(0.25, maximum([maximum(err) for err in stress_spatial_errors_n])*0.95, 
         text("Negative Electrode", 10, :center, :bold))
annotate!(0.75, maximum([maximum(err) for err in stress_spatial_errors_p])*0.95, 
         text("Positive Electrode", 10, :center, :bold))

savefig(p5, "us06_concentration_spatial_error.png")
savefig(p6, "us06_stress_spatial_error.png")

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

println("=== US06工况空间分布分析结果 ===")
println("分析时间点: ", analysis_times, " 秒")
println("负极最大空间浓度误差: ", max_spatial_conc_errors_n, " mol/m³")
println("负极平均空间浓度误差: ", mean_spatial_conc_errors_n, " mol/m³")
println("负极最大空间应力误差: ", max_spatial_stress_errors_n, " Pa")
println("负极平均空间应力误差: ", mean_spatial_stress_errors_n, " Pa")
println("正极最大空间浓度误差: ", max_spatial_conc_errors_p, " mol/m³")
println("正极平均空间浓度误差: ", mean_spatial_conc_errors_p, " mol/m³")
println("正极最大空间应力误差: ", max_spatial_stress_errors_p, " Pa")
println("正极平均空间应力误差: ", mean_spatial_stress_errors_p, " Pa")

# 总体最大误差
overall_max_conc_error = maximum([maximum(err) for err in all_conc_errors])
overall_max_stress_error = maximum([maximum(err) for err in all_stress_errors])
overall_mean_conc_error = mean([mean(err) for err in all_conc_errors])
overall_mean_stress_error = mean([mean(err) for err in all_stress_errors])

println("\n=== 总体统计 ===")
println("总体最大浓度误差: $overall_max_conc_error mol/m³")
println("总体平均浓度误差: $overall_mean_conc_error mol/m³")
println("总体最大应力误差: $overall_max_stress_error Pa")
println("总体平均应力误差: $overall_mean_stress_error Pa")

# 分析US06工况特性对空间分布的影响
println("\n=== US06工况特性分析 ===")
early_indices = [1, 2, 3, 4, 5]  # 选择前5个时间点分析早期特性
late_indices = [6, 7, 8, 9, 10]  # 选择后5个时间点作为对比

early_conc_errors_n = [concentration_spatial_errors_n[i] for i in early_indices]
late_conc_errors_n = [concentration_spatial_errors_n[i] for i in late_indices]
early_conc_errors_p = [concentration_spatial_errors_p[i] for i in early_indices]
late_conc_errors_p = [concentration_spatial_errors_p[i] for i in late_indices]

early_max_error_n = maximum([maximum(err) for err in early_conc_errors_n])
late_max_error_n = maximum([maximum(err) for err in late_conc_errors_n])
early_max_error_p = maximum([maximum(err) for err in early_conc_errors_p])
late_max_error_p = maximum([maximum(err) for err in late_conc_errors_p])

println("负极早期段最大浓度误差: $early_max_error_n mol/m³")
println("负极后期段最大浓度误差: $late_max_error_n mol/m³")
println("负极误差变化系数: $(early_max_error_n / late_max_error_n)")
println("正极早期段最大浓度误差: $early_max_error_p mol/m³")
println("正极后期段最大浓度误差: $late_max_error_p mol/m³")
println("正极误差变化系数: $(early_max_error_p / late_max_error_p)")

# 保存空间分析数据到CSV
# 为负极创建DataFrame
negative_df = DataFrame()
negative_df[!, "Position"] = x_grid_n
negative_df[!, "Normalized_Position"] = x_grid_n_normalized

for i in 1:length(time_indices)
    t = analysis_times[i]
    negative_df[!, "Conc_P2D_t$(t)s"] = concentration_p2d_spatial[i]
    negative_df[!, "Conc_sP2D_t$(t)s"] = concentration_sP2D_spatial[i]
    negative_df[!, "Conc_Error_t$(t)s"] = concentration_spatial_errors_n[i]
    negative_df[!, "Stress_P2D_t$(t)s"] = stress_p2d_spatial[i]
    negative_df[!, "Stress_sP2D_t$(t)s"] = stress_sP2D_spatial[i]
    negative_df[!, "Stress_Error_t$(t)s"] = stress_spatial_errors_n[i]
    # 为每个位置添加对应的电流值（重复标量值）
    negative_df[!, "Current_t$(t)s"] = fill(current_interp(t), length(x_grid_n))
end

# 为正极创建DataFrame
positive_df = DataFrame()
positive_df[!, "Position"] = x_grid_p
positive_df[!, "Normalized_Position"] = x_grid_p_normalized

for i in 1:length(time_indices)
    t = analysis_times[i]
    positive_df[!, "Conc_P2D_t$(t)s"] = concentration_p_p2d_spatial[i]
    positive_df[!, "Conc_sP2D_t$(t)s"] = concentration_p_sP2D_spatial[i]
    positive_df[!, "Conc_Error_t$(t)s"] = concentration_spatial_errors_p[i]
    positive_df[!, "Stress_P2D_t$(t)s"] = stress_p_p2d_spatial[i]
    positive_df[!, "Stress_sP2D_t$(t)s"] = stress_p_sP2D_spatial[i]
    positive_df[!, "Stress_Error_t$(t)s"] = stress_spatial_errors_p[i]
    # 为每个位置添加对应的电流值（重复标量值）
    positive_df[!, "Current_t$(t)s"] = fill(current_interp(t), length(x_grid_p))
end

# 保存两个DataFrame
CSV.write("us06_spatial_analysis_negative_electrode.csv", negative_df)
CSV.write("us06_spatial_analysis_positive_electrode.csv", positive_df)

println("US06工况空间分布分析完成，结果已保存")

# 引用信息
JuBat.Citation()
