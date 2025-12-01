using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
include("../src/JuBat.jl") 

# ======================
# 1. 电池参数与初始化
# ======================
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.3  # 设置截止电压

# ======================
# 2. 定义火星昼夜周期参数 (基于Perseverance任务数据)
# ======================
const mars_day_seconds = 88775  # 火星日长度(s)
sim_duration = 1 * mars_day_seconds  # 先模拟0.5个火星日（测试稳定性）
time_step = 2  # 增大时间步长提高稳定性
time = 0:time_step:sim_duration

# ======================
# 3. 温度-电流耦合模型(简化火星昼夜温度曲线)
# ======================
function get_mars_temp(t)
    """ 火星昼夜温度模型 (-90℃夜间到0℃白天) """
    -45 - 45*cos(2π*t/mars_day_seconds) + 3*randn()  # 添加随机波动
end

function smooth_transition(t, t_start, t_end, val_start, val_end)
    """ 余弦平滑过渡函数 """
    k = clamp((t - t_start) / (t_end - t_start), 0.0, 1.0)
    val_start + (val_end - val_start) * (1 - cos(π*k))/2
end

function generate_mars_current(t, voltage_history=[])
    """ 生成符合物理约束的火星电流工况 """
    current_scale_factor = 0.05  # 经验值，根据模拟调整
    temp = get_mars_temp(t)
    
    # 温度相关电流限制
    max_charge = temp > -10 ? 1.0 : 0.0  # -10℃以下禁止充电
    max_discharge = temp < -60 ? 0.05 : 0.2
    
    # 基础电流模式
    if temp < -70  # 极寒模式
        base = 0.03 * (1 + 0.03randn())  # 极小电流维持
    elseif -70 ≤ temp < -40  # 过渡区1
        base = smooth_transition(t, findlast(t -> get_mars_temp(t) < -70, time), 
                                findfirst(t -> get_mars_temp(t) ≥ -40, time), 
                                0.03, max_discharge)
    elseif -40 ≤ temp < -10  # 过渡区2
        base = smooth_transition(t, findlast(t -> get_mars_temp(t) < -40, time), 
                                findfirst(t -> get_mars_temp(t) ≥ -10, time), 
                                max_discharge, -max_charge)
    else  # 正常工作区
        base = -max_charge * (1 + 0.02randn())
    end
    

    # 添加设备脉冲（幅度随温度调整）
    pulse_amp = temp < -60 ? 0.01 : min(0.2, 0.05 + 0.15*(temp + 90)/90)
    if mod(t, 1800) < 10  # 每30分钟工作10秒
        base -= pulse_amp * sin(π*t/5)^2  # 缓启动脉冲
    end
    
    # 温度-阻抗耦合效应
    (base * exp(0.03*(25 - temp)) + 0.01randn()) * current_scale_factor  # Arrhenius修正+噪声
end

# 生成电流序列
current = [generate_mars_current(t) for t in time]
current_interp = LinearInterpolation(time, current, extrapolation_bc=Flat())

# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = time  # 使用与电流数据相同的时间步长

# 设置电流函数
opt.Current = t -> current_interp(t)

# 创建和运行P2D模型模拟
println("开始模拟NASA火星探测器电流工况... P2D模型")
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
concentration_p2d = result["negative particle surface lithium concentration [mol/m^3]"]
stress_p2d = result["negative particle surface tangential stress[Pa]"]

# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d = concentration_p2d[1, :]
end

if ndims(stress_p2d) > 1
    stress_p2d = stress_p2d[1, :]
end

# 设置sP2D模型
opt.model = "sP2D"  # 使用sP2D模型
println("开始模拟NASA火星探测器电流工况... sP2D模型")
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

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D = concentration_sP2D[1, :]
end

if ndims(stress_sP2D) > 1
    stress_sP2D = stress_sP2D[1, :]
end

# Find the minimum length among all result vectors
min_length = min(length(time_result), 
                length(voltage_p2d), length(voltage_sP2D),
                length(concentration_p2d), length(concentration_sP2D),
                length(stress_p2d), length(stress_sP2D))

# Truncate all vectors to this minimum length
time_result = time_result[1:min_length]
voltage_p2d = voltage_p2d[1:min_length]
voltage_sP2D = voltage_sP2D[1:min_length]
concentration_p2d = concentration_p2d[1:min_length]
concentration_sP2D = concentration_sP2D[1:min_length]
stress_p2d = stress_p2d[1:min_length]
stress_sP2D = stress_sP2D[1:min_length]

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time_result, voltage_p2d, 
    label="voltage (P2D)", 
    xlabel="time [s]", 
    ylabel="voltage [V]", 
    lw=2, color=:red)

plot!(time_result, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:purple)

# 创建温度轴（右侧Y轴）
temp_plot = plot!(twinx(), time_result, get_mars_temp.(time_result),
label="Temperature", 
ylabel="Temperature [℃]",
color=:jet, lw=1.5, linestyle=:dot,
legend=:topright,
grid=false) # 禁用第二条Y轴的网格线

# 创建浓度图，包括P2D和sP2D的结果
p2 = plot(time_result, concentration_p2d, 
    label="concentration (P2D)", 
    xlabel="time [s]", 
    ylabel="concentration [mol/m³]", 
    lw=2, color=:red)

plot!(time_result, concentration_sP2D, 
    label="concentration (sP2D)", 
    title="concentration_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:fuchsia)

# 创建温度轴（右侧Y轴）
temp_plot = plot!(twinx(), time_result, get_mars_temp.(time_result),
label="Temperature", 
ylabel="Temperature [℃]",
color=:jet, lw=1.5, linestyle=:dot,
legend=:topright,
grid=false) # 禁用第二条Y轴的网格线

# 创建应力图，包括P2D和sP2D的结果
p3 = plot(time_result, stress_p2d, 
    label="stress (P2D)", 
    xlabel="time [s]", 
    title="stress_comparison (P2D vs sP2D)",
    lw=2, color=:red)

plot!(time_result, stress_sP2D, 
    label="stress (sP2D)", 
    ylabel="stress [Pa]", 
    linestyle=:dash, color=:darkorange)

# 创建温度轴（右侧Y轴）
temp_plot = plot!(twinx(), time_result, get_mars_temp.(time_result),
label="Temperature", 
ylabel="Temperature [℃]",
color=:jet, lw=1.5, linestyle=:dot,
legend=:topright,
grid=false) # 禁用第二条Y轴的网格线

current_plot_length = min(length(time), length(current))

# 创建电流图
p4 = plot(time[1:current_plot_length], current[1:current_plot_length], 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile",
    lw=2, size=(1600, 400), color=:orange)

# 创建温度轴（右侧Y轴）
temp_plot = plot!(twinx(), time_result, get_mars_temp.(time_result),
label="Temperature", 
ylabel="Temperature [℃]",
color=:jet, lw=1.5, linestyle=:dot,
legend=:topright,
grid=false) # 禁用第二条Y轴的网格线

# 组合图表
plot_combined = plot(p1, p2, p3, p4, 
layout=(4, 1), 
size=(800, 800))

# 保存图表
savefig(p1, "Temperature_voltage_comparison.pdf")
savefig(p2, "Temperature_concentration_comparison.pdf")
savefig(p3, "Temperature_stress_comparison.pdf")
savefig(p4, "Temperature_current_profile.pdf")

savefig(plot_combined, "Temperature_battery_analysis_comparison_super.pdf")

# 计算绝对值差
voltage_error = abs.((voltage_sP2D) .- (voltage_p2d))
concentration_error = abs.((concentration_sP2D) .- (concentration_p2d))
stress_error = abs.((stress_sP2D) .- (stress_p2d))

# 计算误差百分比
voltage_error_percentage = (voltage_error ./ voltage_p2d) .* 100
concentration_error_percentage = (concentration_error ./ concentration_p2d) .* 100
stress_error_percentage = (stress_error ./ stress_p2d) .* 100

# 获取最大误差
max_voltage_error = maximum(voltage_error)
max_concentration_error = maximum(concentration_error)
max_stress_error = maximum(stress_error)

# 获取最大误差百分比
max_voltage_error_percentage = maximum(voltage_error_percentage)
max_concentration_error_percentage = maximum(concentration_error_percentage)
max_stress_error_percentage = maximum(stress_error_percentage)

# 计算平均绝对误差
mean_voltage_error = mean(abs.(voltage_error))
mean_concentration_error = mean(abs.(concentration_error))
mean_stress_error = mean(abs.(stress_error))

# 计算平均误差百分比（排除除以零的情况）
function safe_mean_percentage(errors, references)
    valid_indices = findall(x -> x != 0, references)
    isempty(valid_indices) ? 0.0 : mean(abs.(errors[valid_indices] ./ references[valid_indices])) * 100
end

mean_voltage_error_percentage = safe_mean_percentage(voltage_error, voltage_p2d)
mean_concentration_error_percentage = safe_mean_percentage(concentration_error, concentration_p2d)
mean_stress_error_percentage = safe_mean_percentage(stress_error, stress_p2d)

# 确保所有数组长度一致
#min_length = min(length(time), length(voltage_error), 
                #length(concentration_error), length(stress_error))
#time_plot = time_result[1:min_length]
#voltage_error_plot = voltage_error[1:min_length]
#concentration_error_plot = concentration_error[1:min_length]
#stress_error_plot = stress_error[1:min_length]

# 创建误差百分比变化图
p5 = scatter(time_result, voltage_error, 
    label="Voltage Error", 
    xlabel="Time [s]", 
    ylabel="Error[V]",
    title="Model Voltage Error Comparison (P2D vs sP2D)",
    markershape=:diamond, markercolor=:green,
    markersize=2)  # 调整坐标范围

p6 = scatter(time_result, concentration_error,
    label="Concentation Error",
    xlabel="Time [s]", 
    ylabel="Error[mol/m³] ",
    title="Model Conc-Error Comparison (P2D vs sP2D)", 
    markershape=:diamond, markercolor=:blue,
    markersize=2)

# 添加应力误差散点（仅有效点）
p7 = scatter(time_result, stress_error,
    label="Stress Error",
    xlabel="Time [s]", 
    ylabel="Error[Pa]",
    title="Model Stress Error Comparison (P2D vs sP2D)", 
    markershape=:diamond, markercolor=:red, 
    markersize=2)

savefig(p5, "Temperature_battery_error_voltage_comparison.pdf")
savefig(p6, "Temperature_battery_error__C_comparison.pdf")
savefig(p7, "Temperature_battery_error__stress_comparison.pdf")

# 打印最大误差和最大误差百分比
println("最大电压误差: $max_voltage_error V")
println("平均电压误差: $mean_voltage_error V")
println("最大电压误差百分比: $max_voltage_error_percentage %")
println("平均电压误差百分比: $mean_voltage_error_percentage %")
println("最大浓度误差: $max_concentration_error mol/m³")
println("平均浓度误差: $mean_concentration_error mol/m³")
println("最大浓度误差百分比: $max_concentration_error_percentage %")
println("平均浓度误差百分比: $mean_concentration_error_percentage %")
println("最大应力误差: $max_stress_error Pa")
println("平均应力误差: $mean_stress_error Pa")
println("最大应力误差百分比: $max_stress_error_percentage %")
println("平均应力误差百分比: $mean_stress_error_percentage %")

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time,
    "Voltage (P2D) (V)" => voltage_p2d,
    "Voltage (sP2D) (V)" => voltage_sP2D,
    "Voltage Error (V)" => voltage_error,
    "average Voltage Error (V)" => mean_voltage_error,
    "Voltage Error Percentage (%)" => voltage_error_percentage,
    "average Voltage Error Percentage (%)" => mean_voltage_error_percentage,
    "Concentration (P2D) (mol/m³)" => concentration_p2d,
    "Concentration (sP2D) (mol/m³)" => concentration_sP2D,
    "Concentration Error (mol/m³)" => concentration_error,
    "average Concentration Error (mol/m³)" => mean_concentration_error,
    "Concentration Error Percentage (%)" => concentration_error_percentage,
    "average Concentration Error Percentage (%)" => mean_concentration_error_percentage,
    "Stress (P2D) (Pa)" => stress_p2d,
    "Stress (sP2D) (Pa)" => stress_sP2D,
    "Stress Error (Pa)" => stress_error,
    "average Stress Error (Pa)" => mean_stress_error,
    "Stress Error Percentage (%)" => stress_error_percentage,
    "average Stress Error Percentage (%)" => mean_stress_error_percentage,
    "Current (A)" => [current_interp(t) for t in time]
)

CSV.write("Temperature_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")

# 引用信息
JuBat.Citation()