using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
using Plots.PlotMeasures
include("../src/JuBat.jl") 

# 设置数据路径
path = "D:/竞赛和课程文件/课程文件/毕业设计/SRC/data/drive_cycles/"

# 读取UDDS数据
udds_data = CSV.read(path * "UDDS1.csv", DataFrame, header = 1)
time_udds = udds_data[:, "Time"]
current_udds = udds_data[:, " Current"]

time = collect(Float64, 0:0.01:maximum(time_udds))

# 创建电流插值函数，用于获取任意时间点的电流值
current_interp = LinearInterpolation(time_udds, current_udds, extrapolation_bc=Flat())

# 设置电池参数
param_dim = JuBat.ChooseCell("Enertech")
param_dim.cell.v_h = 4.3

# UDDS工况测试背景说明
# 基于Tesla电池实际使用工况选择UDDS测试：
# 1. UDDS (Urban Dynamometer Driving Schedule) 是EPA标准城市驾驶循环
# 2. Tesla Model S/X在城市工况下的典型使用场景，包含频繁启停
# 3. 根据SAE J2464和GB/T 31467.3标准，电池测试需要覆盖实际驾驶工况
# 4. UDDS循环包含加速、减速、巡航、怠速等多种工况，模拟真实城市驾驶
# 5. 测试时间约1370秒，能够充分验证电池在动态负载下的性能表现

# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = time  # 以1秒的间隔进行模拟

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
voltage_p2d = result["cell voltage [V]"]
concentration_p2d_raw = result["negative particle surface lithium concentration [mol/m^3]"]
stress_neg_p2d_raw = result["negative particle surface tangential stress[Pa]"]
stress_pos_p2d_raw = result["positive particle surface tangential stress[Pa]"]

# 处理P2D模型结果的维度
if ndims(concentration_p2d_raw) > 1
    # 计算整个负极的平均浓度（沿空间维度平均）
    concentration_p2d_raw_avg = mean(concentration_p2d_raw, dims=1)[1, :]
    # 保留第一个位置索引的浓度用于对比
    concentration_p2d_first = concentration_p2d_raw[1, :]
else
    concentration_p2d_raw_avg = concentration_p2d_raw
    concentration_p2d_first = concentration_p2d_raw
end

# 无量纲化浓度 (c/c_max)
# 使用负极最大浓度进行归一化，c/c_max表示荷电状态(SOC)
# c_max_neg = 28700 mol/m³ (Enertech负极最大浓度)
c_max_neg = param_dim.NE.cs_max  # 负极最大浓度
concentration_p2d = concentration_p2d_raw_avg ./ c_max_neg
concentration_p2d_first = concentration_p2d_first ./ c_max_neg

# 处理负极应力数据（取最后一个位置索引）
if ndims(stress_neg_p2d_raw) > 1
    stress_neg_p2d_raw_last = stress_neg_p2d_raw[end, :]  # 负极最后一个位置
else
    stress_neg_p2d_raw_last = stress_neg_p2d_raw
end

# 处理正极应力数据（取第一个位置索引）
if ndims(stress_pos_p2d_raw) > 1
    stress_pos_p2d_raw_first = stress_pos_p2d_raw[1, :]  # 正极第一个位置
else
    stress_pos_p2d_raw_first = stress_pos_p2d_raw
end

# 无量纲化应力 (σ/E)
# 使用负极和正极的杨氏模量进行归一化，σ/E表示相对应力水平
# E_neg = 15 GPa (Enertech负极杨氏模量)
# E_pos = 375 GPa (Enertech正极杨氏模量)
E_neg = param_dim.NE.E  # 负极杨氏模量
E_pos = param_dim.PE.E  # 正极杨氏模量

stress_neg_p2d = stress_neg_p2d_raw_last ./ E_neg
stress_pos_p2d = stress_pos_p2d_raw_first ./ E_pos

# 调试：检查无量纲化是否有效
println("P2D负极应力无量纲化检查:")
println("  原始应力范围: $(minimum(stress_neg_p2d_raw_last)) 到 $(maximum(stress_neg_p2d_raw_last)) Pa")
println("  杨氏模量: $E_neg Pa")
println("  无量纲化后范围: $(minimum(stress_neg_p2d)) 到 $(maximum(stress_neg_p2d))")
println("P2D正极应力无量纲化检查:")
println("  原始应力范围: $(minimum(stress_pos_p2d_raw_first)) 到 $(maximum(stress_pos_p2d_raw_first)) Pa")
println("  杨氏模量: $E_pos Pa")
println("  无量纲化后范围: $(minimum(stress_pos_p2d)) 到 $(maximum(stress_pos_p2d))")

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
println("sP2D模型计算时间: ", sp2d_time, " 秒")


# 提取sP2D模型结果数据
voltage_sP2D = result_sP2D["cell voltage [V]"]
concentration_sP2D_raw = result_sP2D["negative particle surface lithium concentration [mol/m^3]"]
stress_neg_sP2D_raw = result_sP2D["negative particle surface tangential stress[Pa]"]
stress_pos_sP2D_raw = result_sP2D["positive particle surface tangential stress[Pa]"]

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D_raw) > 1
    # 计算整个负极的平均浓度（沿空间维度平均）
    concentration_sP2D_raw_avg = mean(concentration_sP2D_raw, dims=1)[1, :]
    # 保留第一个位置索引的浓度用于对比
    concentration_sP2D_first = concentration_sP2D_raw[1, :]
else
    concentration_sP2D_raw_avg = concentration_sP2D_raw
    concentration_sP2D_first = concentration_sP2D_raw
end

# 无量纲化浓度 (c/c_max)
concentration_sP2D = concentration_sP2D_raw_avg ./ c_max_neg
concentration_sP2D_first = concentration_sP2D_first ./ c_max_neg

# 处理负极应力数据（取最后一个位置索引）
if ndims(stress_neg_sP2D_raw) > 1
    stress_neg_sP2D_raw_last = stress_neg_sP2D_raw[end, :]  # 负极最后一个位置
else
    stress_neg_sP2D_raw_last = stress_neg_sP2D_raw
end

# 处理正极应力数据（取第一个位置索引）
if ndims(stress_pos_sP2D_raw) > 1
    stress_pos_sP2D_raw_first = stress_pos_sP2D_raw[1, :]  # 正极第一个位置
else
    stress_pos_sP2D_raw_first = stress_pos_sP2D_raw
end

# 无量纲化应力 (σ/E)
stress_neg_sP2D = stress_neg_sP2D_raw_last ./ E_neg
stress_pos_sP2D = stress_pos_sP2D_raw_first ./ E_pos

# 计算整个时间范围的RMSE
voltage_error = voltage_sP2D .- voltage_p2d
rmse_voltage_full = sqrt(mean(voltage_error.^2))

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time_result, voltage_p2d, 
    label="voltage (P2D)", 
    xlabel="time [s]", 
    ylabel="voltage [V]", 
    lw=2, color=:blue)

plot!(time_result, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:red,dpi = 600)

# 添加RMSE标注
annotate!(500, maximum(voltage_p2d) - 0.23,  # 调整位置往下移动
    text("RMSE = $(round(rmse_voltage_full, digits=6)) V", 
    :left, 10, :black))

# 创建浓度图，包括P2D和sP2D的结果
p2 = plot(time_result, concentration_p2d, 
    label="average concentration (P2D)", 
    xlabel="time [s]", 
    ylabel="c/c_max (dimensionless)", 
    lw=2, color=:blue,dpi = 600)

plot!(time_result, concentration_sP2D, 
    label="average concentration (sP2D)", 
    title="average concentration_comparison (P2D vs sP2D) - Normalized",
    titlefontsize=10,
    linestyle=:dash, color=:fuchsia)

# 创建负极应力图，包括P2D和sP2D的结果
p3 = plot(time_result, stress_neg_p2d, 
    label="negative stress (P2D)", 
    xlabel="time [s]", 
    title="negative stress comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:blue,dpi = 600)

plot!(time_result, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    ylabel="σ/E (dimensionless)", 
    linestyle=:dash, color=:darkorange)

# 创建正极应力图，包括P2D和sP2D的结果
p3b = plot(time_result, stress_pos_p2d, 
    label="positive stress (P2D)", 
    xlabel="time [s]", 
    title="positive stress comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:purple,dpi = 600)

plot!(time_result, stress_pos_sP2D, 
    label="positive stress (sP2D)", 
    ylabel="σ/E (dimensionless)", 
    linestyle=:dash, color=:orange)

# 创建合并的正极和负极应力图
p3_combined = plot(time_result, stress_neg_p2d, 
    label="negative stress (P2D)", 
    xlabel="time [s]", 
    ylabel="σ/E (dimensionless)", 
    title="Combined Stress Comparison with Error (P2D vs sP2D) - Normalized",
    titlefontsize=10,
    legend=:outerright,
    size=(1000, 400),
    left_margin=5mm,
    bottom_margin=4mm,
    lw=2, color=:blue,dpi = 600)

plot!(time_result, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    linestyle=:dash, color=:darkorange)

plot!(time_result, stress_pos_p2d, 
    label="positive stress (P2D)", 
    lw=2, color=:purple)

plot!(time_result, stress_pos_sP2D, 
    label="positive stress (sP2D)", 
    linestyle=:dash, color=:orange)

# 添加误差曲线（需要先计算误差）
stress_neg_error_combined = abs.(stress_neg_sP2D) .- abs.(stress_neg_p2d)
stress_pos_error_combined = abs.(stress_pos_sP2D) .- abs.(stress_pos_p2d)

plot!(time_result, stress_neg_error_combined, 
    label="negative stress error", 
    linestyle=:dot, color=:red, lw=1.5)

plot!(time_result, stress_pos_error_combined, 
    label="positive stress error", 
    linestyle=:dot, color=:green, lw=1.5)

# 创建双y轴版本的合并应力图（解决量级差异问题）
p3_combined_dual = plot(time_result, stress_neg_p2d, 
    label="negative stress (P2D)", 
    xlabel="time [s]", 
    ylabel="Stress σ/E (dimensionless)", 
    title="Combined Stress Comparison with Error (Dual Y-axis) - Normalized",
    titlefontsize=10,
    legend=:outertopleft,
    size=(1000, 400),
    left_margin=5mm,
    bottom_margin=4mm,
    lw=2, color=:blue,dpi = 600)

plot!(time_result, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    linestyle=:dash, color=:darkorange)

plot!(time_result, stress_pos_p2d, 
    label="positive stress (P2D)", 
    lw=2, color=:purple)

plot!(time_result, stress_pos_sP2D, 
    label="positive stress (sP2D)", 
    linestyle=:dash, color=:orange)

# 创建右y轴用于显示误差
p3_combined_twin = twinx(p3_combined_dual)
plot!(p3_combined_twin, xlabel="", ylabel="Stress Error σ/E (dimensionless)", 
      legend=:outertopright, framestyle=:box)

# 右y轴：绘制误差分布
plot!(p3_combined_twin, time_result, stress_neg_error_combined, 
    label="negative stress error", 
    linestyle=:dot, color=:red, lw=1.5)

plot!(p3_combined_twin, time_result, stress_pos_error_combined, 
    label="positive stress error", 
    linestyle=:dot, color=:green, lw=1.5)

# 创建电流图
p4 = plot(time_udds, current_udds, 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="UDDS Waveform Current Profile",
    lw=2, size=(1600, 400), color=:orange,dpi = 600)

# 创建0-0.1s电流局部放大图
time_zoom = time_udds[time_udds .<= 0.1]
current_zoom = current_udds[time_udds .<= 0.1]
p4b = plot(time_zoom, current_zoom, 
    label="Current (0-100s)", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="UDDS Waveform Current Profile (0-100s Zoom)",
    lw=2, size=(1600, 400), color=:orange)

# 创建0-100s电压局部放大图
# 获取与时间范围对应的电压数据索引
zoom_indices = findall(x -> x <= 100, time_result)
time_voltage_zoom = time_result[zoom_indices]
voltage_p2d_zoom = voltage_p2d[zoom_indices]
voltage_sp2d_zoom = voltage_sP2D[zoom_indices]

# 计算0-100s区间的RMSE
zoom_voltage_error = voltage_sp2d_zoom .- voltage_p2d_zoom
rmse_voltage_zoom = sqrt(mean(zoom_voltage_error.^2))

p1b = plot(time_voltage_zoom, voltage_p2d_zoom,
    label="voltage (P2D)",
    xlabel="time [s]",
    ylabel="voltage [V]",
    title="Voltage Comparison (0-0.1s Zoom)",
    lw=2, color=:blue, dpi=600,
    # 增大字体大小
    titlefontsize=32,        # 标题字体
    tickfontsize=28,         # 刻度字体
    guidefontsize=28,        # 坐标轴标签字体
    legendfontsize=28,       # 图例字体
    )

# 添加RMSE标注
annotate!(70, minimum(voltage_p2d_zoom) + 0.12,  # x坐标从50改为70，y坐标偏移从0.2改为0.1
    text("RMSE = $(round(rmse_voltage_zoom, digits=6)) V", 
    :left, 28, :black))

plot!(time_voltage_zoom, voltage_sp2d_zoom,
    label="voltage (sP2D)",
    linestyle=:dash, color=:red)

# 调整图表布局
plot!(size=(1600, 1300),
    left_margin=5mm,
    bottom_margin=4mm)

# 保存局部放大图
savefig(p1b, "UDDS_voltage_comparison_zoom_0_100s.png")

# 创建新的组合图表，包含局部放大图
plot_combined_with_zoom = plot(p1, p1b, p4, p4b,
    layout=(4, 1),
    size=(1200, 1000))

# 保存包含局部放大图的组合图
savefig(plot_combined_with_zoom, "UDDS_battery_analysis_with_zoom.png")

println("\n局部放大图已生成：")
println("- UDDS_voltage_comparison_zoom_0_100s.png (0-0.1s电压放大图)")
println("- UDDS_battery_analysis_with_zoom.png (包含局部放大图的组合图)")
println("\n图表显示特点：")
println("  * 时间范围：0-0.1s")
println("  * 蓝色实线：P2D模型电压")
println("  * 红色虚线：sP2D模型电压")
println("  * 分辨率：600dpi")
println("  * 图表尺寸：1600×400像素")

# 组合图表（基本结果）
plot_combined = plot(p1, p2, p3, p3b, p4, 
layout=(5, 1), 
size=(800, 1000))

# 组合图表（包含合并应力图）
plot_combined_with_merged = plot(p1, p2, p3_combined, p4, 
layout=(4, 1), 
size=(1200, 800))

# 组合图表（包含双y轴合并应力图）
plot_combined_with_dual = plot(p1, p2, p3_combined_dual, p4, 
layout=(4, 1), 
size=(1200, 800))

# 保存基本图表
savefig(p1, "UDDS_voltage_comparison.png")
savefig(p2, "UDDS_average_concentration_comparison.png")
savefig(p3, "UDDS_negative_stress_comparison.png")
savefig(p3b, "UDDS_positive_stress_comparison.png")
savefig(p3_combined, "UDDS_combined_stress_comparison.png")
savefig(p3_combined_dual, "UDDS_combined_stress_dual_axis.png")
savefig(p4, "UDDS_current_profile.png")
savefig(p4b, "UDDS_current_profile_zoom_0_100s.png")

savefig(plot_combined, "UDDS_battery_analysis_comparison_super.png")
savefig(plot_combined_with_merged, "UDDS_battery_analysis_with_merged_stress.png")
savefig(plot_combined_with_dual, "UDDS_battery_analysis_with_dual_axis_stress.png")

# 计算绝对值差
voltage_error = abs.(voltage_sP2D) .- abs.(voltage_p2d)
concentration_error = abs.(concentration_sP2D) .- abs.(concentration_p2d)
stress_neg_error = abs.(stress_neg_sP2D) .- abs.(stress_neg_p2d)
stress_pos_error = abs.(stress_pos_sP2D) .- abs.(stress_pos_p2d)

# 计算误差百分比
voltage_error_percentage = (voltage_error ./ voltage_p2d) .* 100
concentration_error_percentage = (concentration_error ./ concentration_p2d) .* 100
stress_neg_error_percentage = (stress_neg_error ./ stress_neg_p2d) .* 100
stress_pos_error_percentage = (stress_pos_error ./ stress_pos_p2d) .* 100

# 获取最大误差
max_voltage_error = maximum(voltage_error)
max_concentration_error = maximum(concentration_error)
max_stress_neg_error = maximum(stress_neg_error)
max_stress_pos_error = maximum(stress_pos_error)

# 获取最大误差百分比
max_voltage_error_percentage = maximum(voltage_error_percentage)
max_concentration_error_percentage = maximum(concentration_error_percentage)
max_stress_neg_error_percentage = maximum(stress_neg_error_percentage)
max_stress_pos_error_percentage = maximum(stress_pos_error_percentage)

# 计算平均绝对误差
mean_voltage_error = mean(abs.(voltage_error))
mean_concentration_error = mean(abs.(concentration_error))
mean_stress_neg_error = mean(abs.(stress_neg_error))
mean_stress_pos_error = mean(abs.(stress_pos_error))

# 计算平均误差百分比（排除除以零的情况）
function safe_mean_percentage(errors, references)
    valid_indices = findall(x -> x != 0, references)
    isempty(valid_indices) ? 0.0 : mean(abs.(errors[valid_indices] ./ references[valid_indices])) * 100
end

mean_voltage_error_percentage = safe_mean_percentage(voltage_error, voltage_p2d)
mean_concentration_error_percentage = safe_mean_percentage(concentration_error, concentration_p2d)
mean_stress_neg_error_percentage = safe_mean_percentage(stress_neg_error, stress_neg_p2d)
mean_stress_pos_error_percentage = safe_mean_percentage(stress_pos_error, stress_pos_p2d)

# 确保所有数组长度一致
min_length = min(length(time), length(voltage_error), 
                length(concentration_error), length(stress_neg_error), length(stress_pos_error))
time_plot = time[1:min_length]
voltage_error_plot = voltage_error[1:min_length]
concentration_error_plot = concentration_error[1:min_length]
stress_neg_error_plot = stress_neg_error[1:min_length]
stress_pos_error_plot = stress_pos_error[1:min_length]

# 创建电压误差绝对值变化图
p5 =scatter(time_plot, voltage_error_plot, 
    label="Voltage Error (V)", 
    xlabel="Time [s]", 
    ylabel="Error (V)",
    title="Voltage Error Comparison (P2D vs sP2D)",markersize=2,
    lw=2, color=:green,dpi = 600)

# 创建浓度误差绝对值变化图
p6 = plot(time_plot, concentration_error_plot, 
    label="Concentration Error (dimensionless)", 
    xlabel="Time [s]", 
    ylabel="Error (c/c_max)",
    title="Concentration Error Comparison (P2D vs sP2D) - Normalized",markersize=2,
    lw=2, color=:blue,dpi = 600)

# 创建负极应力误差绝对值变化图
p7 = scatter(time_plot, stress_neg_error_plot, 
    label="Negative Stress Error (dimensionless)", 
    xlabel="Time [s]", 
    ylabel="Error (σ/E)",
    title="Negative Stress Error Comparison (P2D vs sP2D) - Normalized",markersize=2,
    lw=2, color=:purple,dpi = 600)

# 创建正极应力误差绝对值变化图
p8 = scatter(time_plot, stress_pos_error_plot, 
    label="Positive Stress Error (dimensionless)", 
    xlabel="Time [s]", 
    ylabel="Error (σ/E)",
    title="Positive Stress Error Comparison (P2D vs sP2D) - Normalized",markersize=2,
    lw=2, color=:red,dpi = 600)

# 组合误差图表
plot_error_combined = plot(p5, p6, p7, p8, 
layout=(4, 1), 
size=(800, 800))

# 保存误差图
savefig(p5, "UDDS_battery_error_voltage_absolute.png")
savefig(p6, "UDDS_battery_error_concentration_absolute.png")
savefig(p7, "UDDS_battery_error_negative_stress_absolute.png")
savefig(p8, "UDDS_battery_error_positive_stress_absolute.png")
savefig(plot_error_combined, "UDDS_battery_error_analysis_absolute_super.png")

# 打印最大误差和最大误差百分比
println("最大电压误差: $max_voltage_error V")
println("平均电压误差: $mean_voltage_error V")
println("最大电压误差百分比: $max_voltage_error_percentage %")
println("平均电压误差百分比: $mean_voltage_error_percentage %")
println("最大平均浓度误差: $max_concentration_error (dimensionless)")
println("平均浓度误差: $mean_concentration_error (dimensionless)")
println("最大平均浓度误差百分比: $max_concentration_error_percentage %")
println("平均浓度误差百分比: $mean_concentration_error_percentage %")
println("最大负极应力误差: $max_stress_neg_error (dimensionless)")
println("平均负极应力误差: $mean_stress_neg_error (dimensionless)")
println("最大负极应力误差百分比: $max_stress_neg_error_percentage %")
println("平均负极应力误差百分比: $mean_stress_neg_error_percentage %")
println("最大正极应力误差: $max_stress_pos_error (dimensionless)")
println("平均正极应力误差: $mean_stress_pos_error (dimensionless)")
println("最大正极应力误差百分比: $max_stress_pos_error_percentage %")
println("平均正极应力误差百分比: $mean_stress_pos_error_percentage %")

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time_result,
    "Voltage (P2D) (V)" => voltage_p2d,
    "Voltage (sP2D) (V)" => voltage_sP2D,
    "Voltage Error (V)" => voltage_error,
    "average Voltage Error (V)" => mean_voltage_error,
    "Voltage Error Percentage (%)" => voltage_error_percentage,
    "average Voltage Error Percentage (%)" => mean_voltage_error_percentage,
    "Average Concentration (P2D) (c/c_max)" => concentration_p2d,
    "Average Concentration (sP2D) (c/c_max)" => concentration_sP2D,
    "First Position Concentration (P2D) (c/c_max)" => concentration_p2d_first,
    "First Position Concentration (sP2D) (c/c_max)" => concentration_sP2D_first,
    "Concentration Error (c/c_max)" => concentration_error,
    "average Concentration Error (c/c_max)" => mean_concentration_error,
    "Concentration Error Percentage (%)" => concentration_error_percentage,
    "average Concentration Error Percentage (%)" => mean_concentration_error_percentage,
    "Negative Stress (P2D) (σ/E)" => stress_neg_p2d,
    "Negative Stress (sP2D) (σ/E)" => stress_neg_sP2D,
    "Negative Stress Error (σ/E)" => stress_neg_error,
    "average Negative Stress Error (σ/E)" => mean_stress_neg_error,
    "Negative Stress Error Percentage (%)" => stress_neg_error_percentage,
    "average Negative Stress Error Percentage (%)" => mean_stress_neg_error_percentage,
    "Positive Stress (P2D) (σ/E)" => stress_pos_p2d,
    "Positive Stress (sP2D) (σ/E)" => stress_pos_sP2D,
    "Positive Stress Error (σ/E)" => stress_pos_error,
    "average Positive Stress Error (σ/E)" => mean_stress_pos_error,
    "Positive Stress Error Percentage (%)" => stress_pos_error_percentage,
    "average Positive Stress Error Percentage (%)" => mean_stress_pos_error_percentage,
    "Current (A)" => [current_interp(t) for t in time_result]
)

CSV.write("UDDS_battery_results_comparison.csv", results_df)

# 引用信息
JuBat.Citation()
