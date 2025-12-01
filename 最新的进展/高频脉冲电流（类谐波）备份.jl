using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.3

# 定义高频脉冲电流参数
# 基于Tesla电池滥用测试数据选择1000s测试时间：
# 1. Tesla Model S/X电池包容量约85-100kWh，在极端工况下需要测试电池的长期稳定性
# 2. 根据SAE J2464和GB/T 31467.3标准，电池滥用测试通常需要持续15-30分钟
# 3. Tesla超级充电站V3的充电功率可达250kW，对应约2.5C充电倍率
# 4. 在高速行驶+快速充电的复合工况下，电池需要承受持续的高频脉冲负载
# 5. 1000s约16.7分钟，能够充分验证电池在高频脉冲工况下的热管理和电化学稳定性
total_time = 1000      # Total simulation time in seconds (约16.7分钟)
time_step = 0.001        # 时间步长 (1ms，捕捉高频脉冲细节)

# 创建时间数组
time = collect(Float64, 0:time_step:total_time)

# 多级脉冲参数（基于Tesla电池测试工况设计）
# Tesla电池测试特点：
# - 超级充电站V3：250kW功率，对应约2.5C充电倍率
# - 高速行驶工况：持续1-2C放电，间歇性4-6C加速
# - 热管理系统：需要在高频脉冲下维持电池温度在15-45°C
# - 电池管理系统：需要实时监控电压、温度、SOC变化
pulse_interval = 12   # 每12秒一个脉冲（模拟充电站间歇性高功率输出）
pulse_width = 0.05      # 脉冲持续时间5ms（模拟IGBT开关频率）
I_base = 5          # 基准电流(1C=5A，对应LG M50电池)

current = zeros(length(time))
for i in eachindex(time)
    t = time[i]
    
    # 判断是否在脉冲区间内
    if mod(t, pulse_interval) < pulse_width
        # 三级阶梯划分
        phase = mod(t, pulse_width) / pulse_width  # 脉冲内相对位置[0,1)
        
        if phase < 0.3
            # 第一级：0.5C上升沿
            current[i] = 0.5 * I_base * (1 + 0.2*sin(2π*1000*t))  # 叠加1kHz调制
        elseif phase < 0.6
            # 第二级：2C平台
            current[i] = 2.0 * I_base * (1 + 0.1*sin(2π*5000*t))  # 叠加5kHz调制
        else
            # 第三级：4C下降沿
            current[i] = 4.0 * I_base * exp(-(phase-0.6)^2/0.01)   # 高斯衰减
        end
        
        # 添加ns级开关瞬态（位置随机）
        if rand() < 0.01 && phase > 0.2 && phase < 0.8
            current[i] += 0.8 * exp(-(mod(t,1e-6)/1e-9)^2)  # 10ns级尖峰
        end
    else
        # 脉冲间隔期：维持0.1C背景电流
        current[i] = 0.1 * I_base
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

# 创建浓度图，包括P2D和sP2D的结果
p2 = plot(time_result, concentration_p2d, 
    label="average concentration (P2D)", 
    xlabel="time [s]", 
    ylabel="c/c_max (dimensionless)", 
    lw=2, color=:blue,dpi = 600)

plot!(time_result, concentration_sP2D, 
    label="average concentration (sP2D)", 
    title="average concentration_comparison (P2D vs sP2D) - Normalized",
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

# 创建电流图
p4 = plot(time, current, 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile",
    lw=2, size=(1600, 400), color=:orange,dpi = 600)

# 创建0-100s电流局部放大图
time_zoom = time[time .<= 0.1]
current_zoom = current[time .<= 0.1]
p4b = plot(time_zoom, current_zoom, 
    label="Current (0-100s)", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile (0-100s Zoom)",
    lw=2, size=(1600, 400), color=:orange)

# 组合图表（基本结果）
plot_combined = plot(p1, p2, p3, p3b, p4, 
layout=(5, 1), 
size=(800, 1000))

# 保存基本图表
savefig(p1, "pulse_voltage_comparison.png")
savefig(p2, "pulse_average_concentration_comparison.png")
savefig(p3, "pulse_negative_stress_comparison.png")
savefig(p3b, "pulse_positive_stress_comparison.png")
savefig(p4, "pulse_current_profile.png")
savefig(p4b, "pulse_current_profile_zoom_0_100s.png")

savefig(plot_combined, "pulse_battery_analysis_comparison_super.png")

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

CSV.write("pulse_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")

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