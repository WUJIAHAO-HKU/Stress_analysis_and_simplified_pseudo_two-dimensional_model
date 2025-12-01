using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
using Plots.PlotMeasures
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("Enertech")
param_dim.cell.v_h = 4.3

# 定义高频脉冲电流参数
# 基于Tesla电池滥用测试数据选择1000s测试时间：
# 1. Tesla Model S/X电池包容量约85-100kWh，在极端工况下需要测试电池的长期稳定性
# 2. 根据SAE J2464和GB/T 31467.3标准，电池滥用测试通常需要持续15-30分钟
# 3. Tesla超级充电站V3的充电功率可达250kW，对应约2.5C充电倍率
# 4. 在高速行驶+快速充电的复合工况下，电池需要承受持续的高频脉冲负载
# 5. 1000s约16.7分钟，能够充分验证电池在高频脉冲工况下的热管理和电化学稳定性
total_time = 1000      # Total simulation time in seconds (约16.7分钟)
time_step = 0.01        # 时间步长 (1ms，捕捉高频脉冲细节)

# 创建时间数组
time = collect(Float64, 0:time_step:total_time)

# 多级脉冲参数（基于Tesla电池测试工况设计）
# Tesla电池测试特点：
# - 超级充电站V3：250kW功率，对应约2.5C充电倍率
# - 高速行驶工况：持续1-2C放电，间歇性4-6C加速
# - 热管理系统：需要在高频脉冲下维持电池温度在15-45°C
# - 电池管理系统：需要实时监控电压、温度、SOC变化
pulse_interval = 5   # 每5秒一个脉冲（模拟充电站间歇性高功率输出）
pulse_width = 0.5      # 脉冲持续时间500ms（模拟IGBT开关频率）
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
println("开始模拟高频脉冲电流工况... sP2D模型")
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

# Add this synchronization code:
min_length = min(
    length(time_result),
    length(voltage_p2d),
    length(voltage_sP2D),
    length(concentration_p2d),
    length(concentration_sP2D),
    length(stress_neg_p2d),
    length(stress_neg_sP2D),
    length(stress_pos_p2d),
    length(stress_pos_sP2D)
)

# Truncate all arrays to the same length
time_result = time_result[1:min_length]
voltage_p2d = voltage_p2d[1:min_length]
voltage_sP2D = voltage_sP2D[1:min_length]
concentration_p2d = concentration_p2d[1:min_length]
concentration_sP2D = concentration_sP2D[1:min_length]
stress_neg_p2d = stress_neg_p2d[1:min_length]
stress_neg_sP2D = stress_neg_sP2D[1:min_length]
stress_pos_p2d = stress_pos_p2d[1:min_length]
stress_pos_sP2D = stress_pos_sP2D[1:min_length]

# 计算整个时间范围的RMSE
voltage_error = voltage_sP2D .- voltage_p2d
rmse_voltage_full = sqrt(mean(voltage_error.^2))

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time_result, voltage_p2d, 
    label="voltage (P2D)", 
    xlabel="time [s]", 
    ylabel="voltage [V]", 
    lw=2, color=:blue,
    # 增大字体大小
    # 增大字体大小
    titlefontsize=45,        # 标题字体
    tickfontsize=32,         # 刻度字体
    guidefontsize=45,        # 坐标轴标签字体
    legendfontsize=32,       # 图例字体
    # 设置四边页边距
    left_margin=15mm,       # 左边距
    right_margin=15mm,      # 右边距
    top_margin=15mm,        # 上边距
    bottom_margin=15mm,      # 下边距
    size=(1600, 1300))

plot!(time_result, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (Enertech)",
    linestyle=:dash, color=:red,dpi = 600)

# 添加竖直虚线标记放大区域 (0s和0.1s)
vline!([0.0, 100], 
    linestyle=:dash, 
    linewidth=2, 
    color=:black, 
    alpha=0.7, 
    label="")

# 添加标注说明这是放大区域
annotate!(50, minimum(voltage_p2d) + 0.27, 
    text("Zoom\nRegion", :center, 30, :black))

# 添加RMSE标注
annotate!(500, maximum(voltage_p2d) - 0.25,  # 调整位置到图表上方
    text("RMSE = $(round(rmse_voltage_full, digits=6)) V", 
    :left, 40, :black))

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
p4 = plot(time, current, 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile",
    lw=2, size=(1600, 400), color=:orange,dpi = 600)

# 创建0-0.1s电流局部放大图
time_zoom = time[time .<= 0.1]
current_zoom = current[time .<= 0.1]
p4b = plot(time_zoom, current_zoom, 
    label="Current (0-100s)", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile (0-100s Zoom)",
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
    title="Voltage Comparison (0-100s Zoom)",
    lw=2, color=:blue, dpi=600,
    # 增大字体大小
    # 增大字体大小
    titlefontsize=45,        # 标题字体
    tickfontsize=32,         # 刻度字体
    guidefontsize=45,        # 坐标轴标签字体
    legendfontsize=32,       # 图例字体
    # 设置四边页边距
    left_margin=15mm,       # 左边距
    right_margin=15mm,      # 右边距
    top_margin=15mm,        # 上边距
    bottom_margin=15mm,      # 下边距
    size=(1600, 1300)
    )

# 添加RMSE标注
annotate!(50, minimum(voltage_p2d_zoom) + 0.2,
    text("RMSE = $(round(rmse_voltage_zoom, digits=6)) V", 
    :left, 40, :black))

plot!(time_voltage_zoom, voltage_sp2d_zoom,
    label="voltage (sP2D)",
    linestyle=:dash, color=:red)

# 调整图表布局
plot!(size=(1600, 1300),
    left_margin=5mm,
    bottom_margin=4mm)

# 保存局部放大图
savefig(p1b, "pulse_voltage_comparison_zoom_0_100s.png")

# 创建新的组合图表，包含局部放大图
plot_combined_with_zoom = plot(p1, p1b, p4, p4b,
    layout=(4, 1),
    size=(1200, 1000))

# 保存包含局部放大图的组合图
savefig(plot_combined_with_zoom, "pulse_battery_analysis_with_zoom.png")

println("\n局部放大图已生成：")
println("- pulse_voltage_comparison_zoom_0_100s.png (0-0.1s电压放大图)")
println("- pulse_battery_analysis_with_zoom.png (包含局部放大图的组合图)")
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
savefig(p1, "pulse_voltage_comparison.png")
savefig(p2, "pulse_average_concentration_comparison.png")
savefig(p3, "pulse_negative_stress_comparison.png")
savefig(p3b, "pulse_positive_stress_comparison.png")
savefig(p3_combined, "pulse_combined_stress_comparison.png")
savefig(p3_combined_dual, "pulse_combined_stress_dual_axis.png")
savefig(p4, "pulse_current_profile.png")
savefig(p4b, "pulse_current_profile_zoom_0_100s.png")

savefig(plot_combined, "pulse_battery_analysis_comparison_super.png")
savefig(plot_combined_with_merged, "pulse_battery_analysis_with_merged_stress.png")
savefig(plot_combined_with_dual, "pulse_battery_analysis_with_dual_axis_stress.png")

# 确保所有数组长度一致
min_length = min(
    length(time_result),
    length(voltage_p2d),
    length(voltage_sP2D),
    length(concentration_p2d),
    length(concentration_sP2D),
    length(stress_neg_p2d),
    length(stress_neg_sP2D),
    length(stress_pos_p2d),
    length(stress_pos_sP2D)
)

# 截断所有数组到相同长度
time_result = time_result[1:min_length]
voltage_p2d = voltage_p2d[1:min_length]
voltage_sP2D = voltage_sP2D[1:min_length]
concentration_p2d = concentration_p2d[1:min_length]
concentration_sP2D = concentration_sP2D[1:min_length]
stress_neg_p2d = stress_neg_p2d[1:min_length]
stress_neg_sP2D = stress_neg_sP2D[1:min_length]
stress_pos_p2d = stress_pos_p2d[1:min_length]
stress_pos_sP2D = stress_pos_sP2D[1:min_length]

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
savefig(p5, "pulse_battery_error_voltage_absolute.png")
savefig(p6, "pulse_battery_error_concentration_absolute.png")
savefig(p7, "pulse_battery_error_negative_stress_absolute.png")
savefig(p8, "pulse_battery_error_positive_stress_absolute.png")
savefig(plot_error_combined, "pulse_battery_error_analysis_absolute_super.png")

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

# 找出所有数据列的最小长度
min_length = min(
    length(time_result),
    length(voltage_p2d),
    length(voltage_sP2D),
    length(concentration_p2d),
    length(concentration_sP2D),
    length(concentration_p2d_first),  # 添加这个
    length(concentration_sP2D_first), # 添加这个
    length(stress_neg_p2d),
    length(stress_neg_sP2D),
    length(stress_pos_p2d),
    length(stress_pos_sP2D)
)

# 截断所有数组到相同长度
time_result = time_result[1:min_length]
voltage_p2d = voltage_p2d[1:min_length]
voltage_sP2D = voltage_sP2D[1:min_length]
concentration_p2d = concentration_p2d[1:min_length]
concentration_sP2D = concentration_sP2D[1:min_length]
concentration_p2d_first = concentration_p2d_first[1:min_length]  # 添加这个
concentration_sP2D_first = concentration_sP2D_first[1:min_length] # 添加这个
stress_neg_p2d = stress_neg_p2d[1:min_length]
stress_neg_sP2D = stress_neg_sP2D[1:min_length]
stress_pos_p2d = stress_pos_p2d[1:min_length]
stress_pos_sP2D = stress_pos_sP2D[1:min_length]

# 重新计算误差数据
voltage_error = abs.(voltage_sP2D) .- abs.(voltage_p2d)
concentration_error = abs.(concentration_sP2D) .- abs.(concentration_p2d)
stress_neg_error = abs.(stress_neg_sP2D) .- abs.(stress_neg_p2d)
stress_pos_error = abs.(stress_pos_sP2D) .- abs.(stress_pos_p2d)

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time_result,
    "Voltage (P2D) (V)" => voltage_p2d,
    "Voltage (sP2D) (V)" => voltage_sP2D,
    "Voltage Error (V)" => voltage_error,
    "average Voltage Error (V)" => repeat([mean_voltage_error], min_length),
    "Voltage Error Percentage (%)" => voltage_error_percentage,
    "average Voltage Error Percentage (%)" => repeat([mean_voltage_error_percentage], min_length),
    "Average Concentration (P2D) (c/c_max)" => concentration_p2d,
    "Average Concentration (sP2D) (c/c_max)" => concentration_sP2D,
    "First Position Concentration (P2D) (c/c_max)" => concentration_p2d_first,
    "First Position Concentration (sP2D) (c/c_max)" => concentration_sP2D_first,
    "Concentration Error (c/c_max)" => concentration_error,
    "average Concentration Error (c/c_max)" => repeat([mean_concentration_error], min_length),
    "Concentration Error Percentage (%)" => concentration_error_percentage,
    "average Concentration Error Percentage (%)" => repeat([mean_concentration_error_percentage], min_length),
    "Negative Stress (P2D) (σ/E)" => stress_neg_p2d,
    "Negative Stress (sP2D) (σ/E)" => stress_neg_sP2D,
    "Negative Stress Error (σ/E)" => stress_neg_error,
    "average Negative Stress Error (σ/E)" => repeat([mean_stress_neg_error], min_length),
    "Negative Stress Error Percentage (%)" => stress_neg_error_percentage,
    "average Negative Stress Error Percentage (%)" => repeat([mean_stress_neg_error_percentage], min_length),
    "Positive Stress (P2D) (σ/E)" => stress_pos_p2d,
    "Positive Stress (sP2D) (σ/E)" => stress_pos_sP2D,
    "Positive Stress Error (σ/E)" => stress_pos_error,
    "average Positive Stress Error (σ/E)" => repeat([mean_stress_pos_error], min_length),
    "Positive Stress Error Percentage (%)" => stress_pos_error_percentage,
    "average Positive Stress Error Percentage (%)" => repeat([mean_stress_pos_error_percentage], min_length),
    "Current (A)" => [current_interp(t) for t in time_result]
)

CSV.write("pulse_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")
println("\n生成的合并应力图表：")
println("- pulse_combined_stress_comparison.png (正极和负极应力合并图，包含误差曲线)")
println("  包含以下6条曲线：")
println("    * 负极应力 P2D (蓝色实线)")
println("    * 负极应力 sP2D (深橙色虚线)")
println("    * 正极应力 P2D (紫色实线)")
println("    * 正极应力 sP2D (橙色虚线)")
println("    * 负极应力误差 (红色点线)")
println("    * 正极应力误差 (绿色点线)")
println("- pulse_combined_stress_dual_axis.png (双y轴版本，解决量级差异问题)")
println("  左y轴：应力数据，右y轴：误差数据")
println("- pulse_battery_analysis_with_merged_stress.png (包含合并应力图的综合布局)")
println("- pulse_battery_analysis_with_dual_axis_stress.png (包含双y轴合并应力图的综合布局)")
println("\n图表显示优化：")
println("  * 标题字号：10pt，更紧凑")
println("  * 图例位置：图外右侧，不遮挡曲线")
println("  * 边距调整：左边距5mm，下边距4mm，确保标题完整显示")
println("  * 图表尺寸：1000×400像素，为图例留出充足空间")

# 测试总结
println("="^60)
println("Tesla电池高频脉冲测试工况分析总结：")
println("1. 测试时间：1000s (16.7分钟) - 符合SAE J2464标准要求")
println("2. 脉冲频率：每12秒一次 - 模拟超级充电站间歇性高功率输出")
println("3. 脉冲宽度：5ms - 模拟IGBT开关频率和快速响应需求")
println("4. 电流倍率：0.1C-4C - 覆盖Tesla车辆实际使用工况")
println("5. 高频调制：1kHz-5kHz - 模拟电力电子器件的开关噪声")
println("6. 测试目的：验证P2D和sP2D模型在高频脉冲下的预测精度")
println("7. 应用场景：电池管理系统设计、热管理优化、寿命预测")
println("8. 无量纲化处理：")
println("   - 浓度：c/c_max (荷电状态SOC，0-1范围)")
println("   - 应力：σ/E (相对应力水平，便于比较)")
println("   - 提高数值稳定性和物理意义清晰度")
println("="^60)

# 引用信息
JuBat.Citation()

# ==================== 新增图表：宽频噪声和三级阶梯电流周期 ====================

# 1. 宽频噪声图
# 提取一个脉冲周期内的噪声数据
pulse_start_time = 0.0  # 从第一个脉冲开始
pulse_end_time = pulse_width
noise_indices = findall(x -> pulse_start_time <= x <= pulse_end_time, time)
noise_time = time[noise_indices]
noise_data = 0.05 * I_base * randn(length(noise_time))  # 重新生成相同的噪声

# 创建宽频噪声图
p_noise = plot(noise_time .* 1000, noise_data ./ I_base,  # 转换为毫秒显示，电流转换为倍率C
    xlabel="Time [ms]", 
    ylabel="Noise Amplitude [C]", 
    title="Wideband Noise Profile (Single Pulse Period)",
    lw=1.5, 
    color=:red,
    markersize=2,
    marker=:circle,
    dpi=600,
    size=(800, 400),
    legend=false)

# 添加零线参考
hline!([0], color=:black, linestyle=:dash, lw=1, label="")

# 添加统计信息
noise_rms = sqrt(mean(noise_data.^2))
noise_peak = maximum(abs.(noise_data))
noise_rms_c = noise_rms / I_base
noise_peak_c = noise_peak / I_base
annotate!(pulse_width*500, maximum(noise_data./I_base)*0.8, 
    text("RMS: $(round(noise_rms_c, digits=4)) C\nPeak: $(round(noise_peak_c, digits=4)) C", 
         :left, 10, :blue))

# 2. 三级阶梯电流周期图
# 创建一个完整的三级阶梯脉冲周期，两端添加背景电流
background_duration = 0.1  # 两端各添加0.1秒的背景电流
total_cycle_duration = pulse_width + 2 * background_duration
cycle_time = collect(0:time_step:total_cycle_duration)
cycle_current = zeros(length(cycle_time))

for i in eachindex(cycle_time)
    t = cycle_time[i]
    
    if t < background_duration
        # 前端背景电流：1.5C
        cycle_current[i] = 1.5 * I_base
    elseif t < background_duration + pulse_width
        # 脉冲期间：三级阶梯
        pulse_t = t - background_duration
        phase = pulse_t / pulse_width  # 脉冲内相对位置[0,1)
        
        if phase < 0.3
            # 第一级：1.5C上升沿 + 1kHz调制
            cycle_current[i] = 3.0 * I_base * (1 + 0.2*sin(2π*1000*pulse_t))
        elseif phase < 0.6
            # 第二级：2C平台 + 5kHz调制
            cycle_current[i] = 4.0 * I_base * (1 + 0.1*sin(2π*5000*pulse_t))
        else
            # 第三级：4C下降沿 + 高斯衰减
            cycle_current[i] = 5.0 * I_base * exp(-(phase-0.6)^2/0.01)
        end
    else
        # 后端背景电流：1.5C
        cycle_current[i] = 1.5 * I_base
    end
    
    # 添加宽频噪声
    cycle_current[i] += 0.05 * I_base * randn()
end

# 创建三级阶梯电流周期图
p_cycle = plot(cycle_time .* 1000, cycle_current ./ I_base,  # 转换为毫秒显示，电流转换为倍率C
    xlabel="Time [ms]", 
    ylabel="Current [C]", 
    title="Three-Stage Pulse Current Profile with Background Current (Single Period)",
    lw=2, 
    color=:blue,
    dpi=600,
    size=(800, 400),
    left_margin=10mm,
    right_margin=10mm,
    bottom_margin=8mm,
    legend=false)

# 添加阶段分隔线和背景电流区域标识
vline!([background_duration * 1000], color=:green, linestyle=:dash, lw=1, label="")
vline!([(background_duration + 0.3 * pulse_width) * 1000], color=:red, linestyle=:dash, lw=1, label="")
vline!([(background_duration + 0.6 * pulse_width) * 1000], color=:red, linestyle=:dash, lw=1, label="")
vline!([(background_duration + pulse_width) * 1000], color=:green, linestyle=:dash, lw=1, label="")

# 添加背景电流区域标注
annotate!(background_duration * 500, maximum(cycle_current./I_base) * 0.95, 
    text("Background\n1.5C", :center, 10, :green))
annotate!((background_duration + pulse_width + background_duration * 0.5) * 1000, maximum(cycle_current./I_base) * 0.95, 
    text("Background\n1.5C", :center, 10, :green))

# 添加阶段标注
annotate!((background_duration + 0.15 * pulse_width) * 1000, maximum(cycle_current./I_base) * 0.9, 
    text("Stage 1\n3.0C + 1kHz", :center, 10, :darkblue))
annotate!((background_duration + 0.45 * pulse_width) * 1000, maximum(cycle_current./I_base) * 0.9, 
    text("Stage 2\n4.0C + 5kHz", :center, 10, :darkblue))
annotate!((background_duration + 0.8 * pulse_width) * 1000, maximum(cycle_current./I_base) * 0.9, 
    text("Stage 3\n5.0C + Gaussian", :center, 10, :darkblue))

# 添加电流值标注（只计算脉冲期间的各阶段平均值）
pulse_start_idx = Int(round(background_duration / time_step)) + 1
pulse_end_idx = Int(round((background_duration + pulse_width) / time_step))
stage1_end_idx = Int(round((background_duration + 0.3 * pulse_width) / time_step))
stage2_end_idx = Int(round((background_duration + 0.6 * pulse_width) / time_step))

stage1_avg = mean(cycle_current[pulse_start_idx:stage1_end_idx])
stage2_avg = mean(cycle_current[stage1_end_idx+1:stage2_end_idx])
stage3_avg = mean(cycle_current[stage2_end_idx+1:pulse_end_idx])
background_avg = mean(cycle_current[1:pulse_start_idx-1])  # 前端背景电流平均值

annotate!((background_duration + 0.15 * pulse_width) * 1000, stage1_avg/I_base - 0.3, 
    text("$(round(stage1_avg/I_base, digits=1)) C", :center, 8, :red))
annotate!((background_duration + 0.45 * pulse_width) * 1000, stage2_avg/I_base - 0.2, 
    text("$(round(stage2_avg/I_base, digits=1)) C", :center, 8, :red))
annotate!((background_duration + 0.8 * pulse_width) * 1000, stage3_avg/I_base - 0.2, 
    text("$(round(stage3_avg/I_base, digits=1)) C", :center, 8, :red))
annotate!(background_duration * 500, background_avg/I_base - 0.3, 
    text("$(round(background_avg/I_base, digits=1)) C", :center, 8, :green))

# 3. 组合图：噪声和电流周期对比
p_combined_new = plot(p_noise, p_cycle, 
    layout=(2, 1), 
    size=(1000, 600),
    title="High-Frequency Pulse Analysis: Noise and Current Profile",
    left_margin=10mm,
    right_margin=10mm,
    bottom_margin=8mm)

# 保存新图表
savefig(p_noise, "pulse_wideband_noise_profile.png")
savefig(p_cycle, "pulse_three_stage_current_cycle.png")
savefig(p_combined_new, "pulse_noise_and_current_cycle_combined.png")

println("\n" * "="^60)
println("新增图表说明：")
println("1. 宽频噪声图 (pulse_wideband_noise_profile.png):")
println("   - 显示单个脉冲周期内的宽频噪声分布")
println("   - 时间轴：毫秒显示，便于观察高频特征")
println("   - 包含RMS和峰值统计信息")
println("   - 红色圆点标记，突出噪声的随机性")
println()
println("2. 三级阶梯电流周期图 (pulse_three_stage_current_cycle.png):")
println("   - 显示完整的三级阶梯脉冲电流波形，两端包含背景电流")
println("   - 前端背景电流 (0-100ms): 1.5C 稳定电流")
println("   - 阶段1 (100-250ms): 1.5C + 1kHz调制")
println("   - 阶段2 (250-400ms): 2.0C + 5kHz调制")
println("   - 阶段3 (400-600ms): 4.0C + 高斯衰减")
println("   - 后端背景电流 (600-700ms): 1.5C 稳定电流")
println("   - 包含阶段边界线、背景电流区域和平均电流值标注")
println("   - 蓝色实线，红色虚线分隔脉冲阶段，绿色虚线分隔背景区域")
println()
println("3. 组合图 (pulse_noise_and_current_cycle_combined.png):")
println("   - 上下布局，便于对比噪声和电流特征")
println("   - 统一的时间轴和样式")
println("   - 适合用于论文或报告展示")
println("="^60)