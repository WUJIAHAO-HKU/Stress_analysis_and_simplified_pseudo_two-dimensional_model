using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
using Plots.PlotMeasures
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.4

# 参数设置
total_time = 7200      # 总时长(秒)
time_step = 0.01       # 10ms步长
time = collect(0:time_step:total_time)
max_time = maximum(time)

# 设置10个等间距的分析时间0点
analysis_times = collect(range(max_time * 0.1, max_time, length=10))

# 电池参数
I_base = 1            # 1C基准电流(1A)
charge_peaks = [2.4, 1.6, 0.8]  # 三阶段充电峰值电流(A)
discharge_peaks = [4.8, 3.2, 1.6]  # 三阶段放电峰值电流(A)
cycle_duration = 300  # 充放电周期300秒

# 阶段时间划分 (单位:秒)
discharge_stages = [60, 120, 180]  # 放电三阶段截止时间：60s/120s/180s
charge_stages = [240, 285, 300]    # 充电三阶段截止时间：240s/285s/300s

current = zeros(length(time))
#temp = 25.0 .+ zeros(length(time))  # 温度模拟(℃)

for i in 1:length(time)
    t = time[i]
    cycle_phase = mod(t, cycle_duration)
    
    # ===== 放电阶段 (0-180s) =====
    if cycle_phase < discharge_stages[end]
        # 阶段1：强放电 (0-60s)
        if cycle_phase < discharge_stages[1]
            I = discharge_peaks[1] * (1 - 0.1*rand())
            
        # 阶段2：中放电 (60-120s)
        elseif cycle_phase < discharge_stages[2]
            I = discharge_peaks[2] * (1 - 0.05*rand())
            
        # 阶段3：弱放电 (120-180s)
        else
            I = discharge_peaks[3] * (1 - 0.02*rand())
        end
        
    # ===== 充电阶段 (180-300s) =====
    else
        # 阶段1：快充 (180-240s)
        if cycle_phase < charge_stages[1]
            I = -charge_peaks[1] * (1 - 0.1*randn())
            
        # 阶段2：脉冲充 (240-285s)
        elseif cycle_phase < charge_stages[2]
            pulse = charge_peaks[2] * (0.8 + 0.2*sign(sin(2π*10*t))) # 10Hz简化脉冲
            I = -pulse * (1 - 0.05*randn())
            
        # 阶段3：涓流充 (285-300s)
        else
            I = -charge_peaks[3] * (1 - 0.02*randn())
        end
    end
    
    current[i] = I
end

# 添加高频噪声（模拟测量噪声）
current .+= 0.1*randn(length(current))

# 创建电流插值函数
itp = linear_interpolation(time, current, extrapolation_bc=Flat())


# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = time  # 使用与电流数据相同的时间步长

# 设置电流函数
opt.Current = t -> itp(t)  # 正确用法

# 创建和运行P2D模型模拟
println("开始模拟快速充放电电流工况... P2D模型")
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
println("开始模拟快速充放电电流工况... sP2D模型")
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
    lw=2, color=:blue,
        # 增大字体大小
    titlefontsize=45,        # 标题字体
    tickfontsize=32,         # 刻度字体
    guidefontsize=45,        # 坐标轴标签字体
    legendfontsize=32,       # 图例字体
    # 设置四边页边距
    left_margin=7mm,       # 左边距
    right_margin=7mm,      # 右边距
    top_margin=7mm,        # 上边距
    bottom_margin=7mm,      # 下边距
    size=(1600, 1300))

plot!(time_result, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:red, dpi=600)

# 添加RMSE标注
annotate!(4800, maximum(voltage_p2d)-0.185,  # 调整位置往下移动
    text("RMSE = $(round(rmse_voltage_full, digits=6)) V", 
    :left, 32, :black))

# 创建浓度图，包括P2D和sP2D的结果
p2 = plot(time_result, concentration_p2d, 
    label="average concentration (P2D)", 
    xlabel="time [s]", 
    ylabel="c/c_max (dimensionless)", 
    lw=2, color=:blue, dpi=600)

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
    lw=2, color=:blue, dpi=600)

plot!(time_result, stress_neg_sP2D, 
    label="negative stress (sP2D)", 
    ylabel="σ/E (dimensionless)", 
    linestyle=:dash, color=:darkorange)

# 创建正极应力图，包括P2D和sP2D的结果
p3b = plot(time_result, stress_pos_p2d, 
    label="positive stress (P2D)", 
    xlabel="time [s]", 
    title="positive stress comparison (P2D vs sP2D) - Normalized",
    lw=2, color=:purple, dpi=600)

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
    lw=2, color=:blue, dpi=600)

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
    lw=2, color=:blue, dpi=600)

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
    lw=2, size=(1600, 400), color=:black, dpi=600,
    # 设置四边页边距（上边距稍加大以容纳外置图例）
    left_margin=14mm,
    right_margin=10mm,
    top_margin=14mm,
    bottom_margin=14mm,
    # 字体
    titlefontsize=21,
    tickfontsize=14,
    guidefontsize=21,
    legendfontsize=14,
    # 将图例移动到图像上方，避免遮挡数据。:outertop 会放在绘图区之外。
    legend=:outertop,
    )

# 添加分析时间点的垂直虚线和交点标记
for t in analysis_times
    # 添加垂直虚线
    vline!([t], color=:gray, alpha=0.4, linestyle=:dash, linewidth=1, label="")
    
    # 计算交点位置（电流值）
    current_value = itp(t)
    
    # 添加红色叉叉标记
    scatter!([t], [current_value], 
             marker=:xcross, color=:red, 
             markersize=16, markerstrokewidth=2, 
             label=(t == analysis_times[1] ? "Spatial Analysis Points" : ""))
end

# 创建0-100s电流局部放大图
time_zoom = time[time .<= 100]
current_zoom = current[time .<= 100]
p4b = plot(time_zoom, current_zoom, 
    label="Current (0-100s)", 
    color=:black,
    xlabel="Time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile (0-100s Zoom)",
    lw=2, size=(1600, 400),
    legend=:outertop,
    top_margin=14mm)

# 添加位于放大范围内的分析时间点
for t in analysis_times
    if t <= 100
        # 添加垂直虚线
        vline!([t], color=:gray, alpha=0.4, linestyle=:dash, linewidth=1, label="")
        
        # 计算交点位置（电流值）
        current_value = itp(t)
        
        # 添加红色叉叉标记
        scatter!([t], [current_value], 
                marker=:cross, color=:red, 
                markersize=8, markerstrokewidth=2, 
                label=(t == analysis_times[1] ? "Spatial Analysis Points" : ""))
    end
end

# 创建0-330s电压局部放大图
# 获取与时间范围对应的电压数据索引
zoom_indices = findall(x -> x <= 330, time_result)
time_voltage_zoom = time_result[zoom_indices]
voltage_p2d_zoom = voltage_p2d[zoom_indices]
voltage_sp2d_zoom = voltage_sP2D[zoom_indices]

# 计算0-330s区间的RMSE
zoom_voltage_error = voltage_sp2d_zoom .- voltage_p2d_zoom
rmse_voltage_zoom = sqrt(mean(zoom_voltage_error.^2))

p1b = plot(time_voltage_zoom, voltage_p2d_zoom,
    label="voltage (P2D)",
    xlabel="time [s]",
    ylabel="voltage [V]",
    title="Voltage Comparison (0-330s Zoom)",
    lw=2, color=:blue, dpi=600,
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
annotate!(30, minimum(voltage_p2d_zoom)+0.125,  # x坐标从50改为70，y坐标偏移从0.2改为0.1
    text("RMSE = $(round(rmse_voltage_zoom, digits=6)) V", 
    :left, 32, :black))

plot!(time_voltage_zoom, voltage_sp2d_zoom,
    label="voltage (sP2D)",
    linestyle=:dash, color=:red)

# 调整图表布局
plot!(size=(1600, 1300),
    left_margin=5mm,
    bottom_margin=4mm)

# 创建新的组合图表，包含局部放大图
plot_combined_with_zoom = plot(p1, p1b, p4, p4b,
    layout=(4, 1),
    size=(1200, 1000))

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
savefig(p1, "dramatic_voltage_comparison.png")
savefig(p2, "dramatic_concentration_comparison.png")
savefig(p3, "dramatic_negative_stress_comparison.png")
savefig(p3b, "dramatic_positive_stress_comparison.png")
savefig(p3_combined, "dramatic_combined_stress_comparison.png")
savefig(p3_combined_dual, "dramatic_combined_stress_dual_axis.png")
savefig(p4, "dramatic_current_profile.png")
savefig(p4b, "dramatic_current_profile_zoom_0_100s.png")

savefig(plot_combined, "dramatic_battery_analysis_comparison_super.png")
savefig(plot_combined_with_merged, "dramatic_battery_analysis_with_merged_stress.png")
savefig(plot_combined_with_dual, "dramatic_battery_analysis_with_dual_axis_stress.png")
savefig(p1b, "dramatic_voltage_comparison_zoom_0_100s.png")
savefig(plot_combined_with_zoom, "dramatic_battery_analysis_with_zoom.png")

println("\n局部放大图已生成：")
println("- dramatic_voltage_comparison_zoom_0_100s.png (0-100s电压放大图)")
println("- dramatic_battery_analysis_with_zoom.png (包含局部放大图的组合图)")
println("\n图表显示特点：")
println("  * 时间范围：0-100s")
println("  * 蓝色实线：P2D模型电压")
println("  * 红色虚线：sP2D模型电压")
println("  * 分辨率：600dpi")
println("  * 图表尺寸：1600×400像素")

# 计算绝对值差
voltage_error = abs.(voltage_sP2D .- voltage_p2d)
concentration_error = abs.(concentration_sP2D .- concentration_p2d)
stress_neg_error = abs.(stress_neg_sP2D .- stress_neg_p2d)
stress_pos_error = abs.(stress_pos_sP2D .- stress_pos_p2d)

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
p5 = scatter(time_plot, voltage_error_plot, 
    label="Voltage Error (V)", 
    xlabel="Time [s]", 
    ylabel="Error (V)",
    title="Voltage Error Comparison (P2D vs sP2D)",
    markersize=2,
    lw=2, color=:green, dpi=600)

# 创建浓度误差绝对值变化图
p6 = plot(time_plot, concentration_error_plot, 
    label="Concentration Error (dimensionless)", 
    xlabel="Time [s]", 
    ylabel="Error (c/c_max)",
    title="Concentration Error Comparison (P2D vs sP2D) - Normalized",
    markersize=2,
    lw=2, color=:blue, dpi=600)

# 创建负极应力误差绝对值变化图
p7 = scatter(time_plot, stress_neg_error_plot, 
    label="Negative Stress Error (dimensionless)", 
    xlabel="Time [s]", 
    ylabel="Error (σ/E)",
    title="Negative Stress Error Comparison (P2D vs sP2D) - Normalized",
    markersize=2,
    lw=2, color=:purple, dpi=600)

# 创建正极应力误差绝对值变化图
p8 = scatter(time_plot, stress_pos_error_plot, 
    label="Positive Stress Error (dimensionless)", 
    xlabel="Time [s]", 
    ylabel="Error (σ/E)",
    title="Positive Stress Error Comparison (P2D vs sP2D) - Normalized",
    markersize=2,
    lw=2, color=:red, dpi=600)

# 组合误差图表
plot_error_combined = plot(p5, p6, p7, p8, 
layout=(4, 1), 
size=(800, 800))

# 保存误差图
savefig(p5, "dramatic_battery_error_voltage_absolute.png")
savefig(p6, "dramatic_battery_error_concentration_absolute.png")
savefig(p7, "dramatic_battery_error_negative_stress_absolute.png")
savefig(p8, "dramatic_battery_error_positive_stress_absolute.png")
savefig(plot_error_combined, "dramatic_battery_error_analysis_absolute_super.png")

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
    "Current (A)" => [itp(t) for t in time_result]
)

CSV.write("dramatic_battery_results_comparison.csv", results_df)

# 打印分析时间点
println("分析时间点: ", analysis_times, " 秒")
println("分析完成，结果已保存")

# 引用信息
JuBat.Citation()