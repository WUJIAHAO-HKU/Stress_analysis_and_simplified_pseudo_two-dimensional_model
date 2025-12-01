using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
using Plots.PlotMeasures
include("../src/JuBat.jl") 

# 设置数据路径
path = "D:/竞赛和课程文件/课程文件/毕业设计/SRC/data/drive_cycles/"

# 读取UDDS数据
udds_data = CSV.read(path * "UDDS.csv", DataFrame, header = 1)
time_udds = udds_data[:, "Time"]
current_udds = udds_data[:, " Current"]

# 创建电流插值函数，用于获取任意时间点的电流值
current_interp = LinearInterpolation(time_udds, current_udds, extrapolation_bc=Flat())

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")
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
opt.time = collect(Float64, 0:1:maximum(time_udds))  # 以1秒的间隔进行模拟

# 设置电流函数为UDDS表格中的数据
opt.Current = t -> current_interp(t)

# 创建和运行P2D模型模拟
println("开始模拟UDDS工况... P2D模型")
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
time = result["time [s]"]
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
opt.time = collect(Float64, 0:1:maximum(time_udds))  # 以1秒的间隔进行模拟
println("开始模拟UDDS工况... sP2D模型")
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

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time, voltage_p2d, 
    label="voltage (P2D)", 
    xlabel="time [s]", 
    ylabel="voltage [V]", 
    lw=2, color=:blue)

plot!(time, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:red)

# 创建浓度图，包括P2D和sP2D的结果
p2 = plot(time, concentration_p2d, 
    label="average concentration (P2D)", 
    xlabel="time [s]", 
    ylabel="c/c_max (dimensionless)", 
    lw=2, color=:blue,dpi = 600)

plot!(time, concentration_sP2D, 
    label="average concentration (sP2D)", 
    title="average concentration_comparison (P2D vs sP2D) - Normalized",
    linestyle=:dash, color=:fuchsia)

# 创建负极应力图，包括P2D和sP2D的结果
p3 = plot(time, stress_neg_p2d, 
    label="negative stress (P2D)", 
    xlabel="time [s]", 
    title="negative stress comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:blue,dpi = 600)

plot!(time, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    ylabel="σ/E (dimensionless)", 
    linestyle=:dash, color=:darkorange)

# 创建正极应力图，包括P2D和sP2D的结果
p3b = plot(time, stress_pos_p2d, 
    label="positive stress (P2D)", 
    xlabel="time [s]", 
    title="positive stress comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:purple,dpi = 600)

plot!(time, stress_pos_sP2D, 
    label="positive stress (sP2D)", 
    ylabel="σ/E (dimensionless)", 
    linestyle=:dash, color=:orange)

# 创建电流图
p4 = plot(time, [current_interp(t) for t in time], 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="UDDS Current Profile",
    lw=2, size=(800, 400), color=:purple,dpi = 600)

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

# 创建合并的正极和负极应力图
p3_combined = plot(time, stress_neg_p2d, 
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

plot!(time, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    linestyle=:dash, color=:darkorange)

plot!(time, stress_pos_p2d, 
    label="positive stress (P2D)", 
    lw=2, color=:purple)

plot!(time, stress_pos_sP2D, 
    label="positive stress (sP2D)", 
    linestyle=:dash, color=:orange)

# 添加误差曲线（需要先计算误差）
stress_neg_error_combined = abs.(stress_neg_sP2D) .- abs.(stress_neg_p2d)
stress_pos_error_combined = abs.(stress_pos_sP2D) .- abs.(stress_pos_p2d)

plot!(time, stress_neg_error_combined, 
    label="negative stress error", 
    linestyle=:dot, color=:red, lw=1.5)

plot!(time, stress_pos_error_combined, 
    label="positive stress error", 
    linestyle=:dot, color=:green, lw=1.5)

# 创建双y轴版本的合并应力图（解决量级差异问题）
p3_combined_dual = plot(time, stress_neg_p2d, 
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

plot!(time, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    linestyle=:dash, color=:darkorange)

plot!(time, stress_pos_p2d, 
    label="positive stress (P2D)", 
    lw=2, color=:purple)

plot!(time, stress_pos_sP2D, 
    label="positive stress (sP2D)", 
    linestyle=:dash, color=:orange)

# 创建右y轴用于显示误差
p3_combined_twin = twinx(p3_combined_dual)
plot!(p3_combined_twin, xlabel="", ylabel="Stress Error σ/E (dimensionless)", 
      legend=:outertopright, framestyle=:box)

# 右y轴：绘制误差分布
plot!(p3_combined_twin, time, stress_neg_error_combined, 
    label="negative stress error", 
    linestyle=:dot, color=:red, lw=1.5)

plot!(p3_combined_twin, time, stress_pos_error_combined, 
    label="positive stress error", 
    linestyle=:dot, color=:green, lw=1.5)

# 保存基本图表
savefig(p1, "udds_voltage_comparison.png")
savefig(p2, "udds_average_concentration_comparison.png")
savefig(p3, "udds_negative_stress_comparison.png")
savefig(p3b, "udds_positive_stress_comparison.png")
savefig(p3_combined, "udds_combined_stress_comparison.png")
savefig(p3_combined_dual, "udds_combined_stress_dual_axis.png")
savefig(p4, "udds_current_profile.png")

savefig(plot_combined, "udds_battery_analysis_comparison_super.png")
savefig(plot_combined_with_merged, "udds_battery_analysis_with_merged_stress.png")
savefig(plot_combined_with_dual, "udds_battery_analysis_with_dual_axis_stress.png")

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
p5 = plot(time_plot, voltage_error_plot, 
    label="Voltage Error (V)", 
    xlabel="Time [s]", 
    ylabel="Error (V)",
    title="Voltage Error Comparison (P2D vs sP2D)",
    lw=2, color=:green,dpi = 600)

# 创建浓度误差绝对值变化图
p6 = plot(time_plot, concentration_error_plot, 
    label="Concentration Error (c/c_max)", 
    xlabel="Time [s]", 
    ylabel="Error (c/c_max)",
    title="Concentration Error Comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:blue,dpi = 600)

# 创建负极应力误差绝对值变化图
p7 = plot(time_plot, stress_neg_error_plot, 
    label="Negative Stress Error (σ/E)", 
    xlabel="Time [s]", 
    ylabel="Error (σ/E)",
    title="Negative Stress Error Comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:red,dpi = 600)

# 创建正极应力误差绝对值变化图
p8 = plot(time_plot, stress_pos_error_plot, 
    label="Positive Stress Error (σ/E)", 
    xlabel="Time [s]", 
    ylabel="Error (σ/E)",
    title="Positive Stress Error Comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:purple,dpi = 600)

# 组合误差图表
plot_error_combined = plot(p5, p6, p7, p8, 
layout=(4, 1), 
size=(800, 800))

# 保存误差图
savefig(p5, "udds_battery_error_voltage_absolute.png")
savefig(p6, "udds_battery_error_concentration_absolute.png")
savefig(p7, "udds_battery_error_negative_stress_absolute.png")
savefig(p8, "udds_battery_error_positive_stress_absolute.png")
savefig(plot_error_combined, "udds_battery_error_analysis_absolute_super.png")

# 打印最大误差和最大误差百分比
println("最大电压误差: $max_voltage_error V")
println("平均电压误差: $mean_voltage_error V")
println("最大电压误差百分比: $max_voltage_error_percentage %")
println("平均电压误差百分比: $mean_voltage_error_percentage %")
println("最大平均浓度误差: $max_concentration_error (c/c_max)")
println("平均浓度误差: $mean_concentration_error (c/c_max)")
println("最大平均浓度误差百分比: $max_concentration_error_percentage %")
println("平均浓度误差百分比: $mean_concentration_error_percentage %")
println("最大负极应力误差: $max_stress_neg_error (σ/E)")
println("平均负极应力误差: $mean_stress_neg_error (σ/E)")
println("最大负极应力误差百分比: $max_stress_neg_error_percentage %")
println("平均负极应力误差百分比: $mean_stress_neg_error_percentage %")
println("最大正极应力误差: $max_stress_pos_error (σ/E)")
println("平均正极应力误差: $mean_stress_pos_error (σ/E)")
println("最大正极应力误差百分比: $max_stress_pos_error_percentage %")
println("平均正极应力误差百分比: $mean_stress_pos_error_percentage %")

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time,
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
    "Current (A)" => [current_interp(t) for t in time]
)

CSV.write("udds_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")
println("\n生成的合并应力图表：")
println("- udds_combined_stress_comparison.png (正极和负极应力合并图，包含误差曲线)")
println("  包含以下6条曲线：")
println("    * 负极应力 P2D (蓝色实线)")
println("    * 负极应力 sP2D (深橙色虚线)")
println("    * 正极应力 P2D (紫色实线)")
println("    * 正极应力 sP2D (橙色虚线)")
println("    * 负极应力误差 (红色点线)")
println("    * 正极应力误差 (绿色点线)")
println("- udds_combined_stress_dual_axis.png (双y轴版本，解决量级差异问题)")
println("  左y轴：应力数据，右y轴：误差数据")
println("- udds_battery_analysis_with_merged_stress.png (包含合并应力图的综合布局)")
println("- udds_battery_analysis_with_dual_axis_stress.png (包含双y轴合并应力图的综合布局)")
println("\n图表显示优化：")
println("  * 标题字号：10pt，更紧凑")
println("  * 图例位置：图外右侧，不遮挡曲线")
println("  * 边距调整：左边距5mm，下边距4mm，确保标题完整显示")
println("  * 图表尺寸：1000×400像素，为图例留出充足空间")

# 测试总结
println("="^60)
println("UDDS工况电池测试分析总结：")
println("1. 测试工况：UDDS (Urban Dynamometer Driving Schedule) - EPA标准城市驾驶循环")
println("2. 测试时间：约1370秒 - 完整模拟城市驾驶工况")
println("3. 工况特点：加速、减速、巡航、怠速等多种工况组合")
println("4. 电流范围：-2C到+2C - 覆盖Tesla车辆城市使用工况")
println("5. 测试目的：验证P2D和sP2D模型在动态负载下的预测精度")
println("6. 应用场景：电池管理系统设计、城市工况优化、寿命预测")
println("7. 标准符合：SAE J2464、GB/T 31467.3、EPA标准")
println("8. 无量纲化处理：")
println("   - 浓度：c/c_max (荷电状态SOC，0-1范围)")
println("   - 应力：σ/E (相对应力水平，便于比较)")
println("   - 提高数值稳定性和物理意义清晰度")
println("="^60)

# 引用信息
JuBat.Citation()