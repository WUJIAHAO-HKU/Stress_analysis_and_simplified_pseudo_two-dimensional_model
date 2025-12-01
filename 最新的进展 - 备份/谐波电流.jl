using Plots, CSV, DataFrames, Interpolations
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.3

total_time = 1000
time_step = 3  # 更小步长以捕捉高频
time = collect(0:time_step:total_time)

# 谐波电流参数
f_base = 50      # 基频(Hz)
I_1C = 5.0       # 1C电流(A)
harmonics = [
    (3, 0.12),   # 3次谐波，幅值12%
    (5, 0.08),   # 5次谐波，8%
    (7, 0.05)    # 7次谐波，5%
]

# 生成含谐波的电流
current = I_1C .* (sin.(2π*f_base.*time) .+ 
           sum(amp.*sin.(2π*order*f_base.*time) for (order,amp) in harmonics))

# 2. 进行THD标准化校准
target_THD = 0.15  # 15%
current = current ./ sqrt(1 + sum(amp^2 for (_,amp) in harmonics)) .* target_THD

# 添加随机噪声（可选）
current .+= 0.05*I_1C .* randn(length(time))

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
println("开始模拟高谐波电流工况... P2D模型")
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
println("P2D模型模拟完成")

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
println("开始模拟谐波电流工况... sP2D模型")
case_sP2D = JuBat.SetCase(param_dim, opt)
result_sP2D = JuBat.Solve(case_sP2D)
println("sP2D模型模拟完成")

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

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time_result, voltage_p2d, 
    label="voltage (P2D)", 
    xlabel="time [s]", 
    ylabel="voltage [V]", 
    lw=2, color=:blue)

plot!(time_result, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:red)

# 创建浓度图，包括P2D和sP2D的结果
p2 = plot(time_result, concentration_p2d, 
    label="concentration (P2D)", 
    xlabel="time [s]", 
    ylabel="concentration [mol/m³]", 
    lw=2, color=:blue)

plot!(time_result, concentration_sP2D, 
    label="concentration (sP2D)", 
    title="concentration_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:green)

# 创建应力图，包括P2D和sP2D的结果
p3 = plot(time_result, stress_p2d, 
    label="stress (P2D)", 
    xlabel="time [s]", 
    title="stress_comparison (P2D vs sP2D)",
    lw=2, color=:blue)

plot!(time_result, stress_sP2D, 
    label="stress (sP2D)", 
    ylabel="stress [Pa]", 
    linestyle=:dash, color=:green)

# 创建电流图
p4 = plot(time, current, 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile",
    lw=2, size=(1600, 400), color=:purple)

# 组合图表
plot_combined = plot(p1, p2, p3, p4, 
layout=(4, 1), 
size=(800, 800))

# 保存图表
savefig(p1, "Harmonic_voltage_comparison.pdf")
savefig(p2, "Harmonic_concentration_comparison.pdf")
savefig(p3, "Harmonic_stress_comparison.pdf")
savefig(p4, "Harmonic_current_profile.pdf")

savefig(plot_combined, "Harmonic_battery_analysis_comparison_super.pdf")

# 计算绝对值差
voltage_error = abs.(voltage_sP2D) .- abs.(voltage_p2d)
concentration_error = abs.(concentration_sP2D) .- abs.(concentration_p2d)
stress_error = abs.(stress_sP2D) .- abs.(stress_p2d)

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

# 打印最大误差和最大误差百分比
println("最大电压误差: $max_voltage_error V")
println("最大电压误差百分比: $max_voltage_error_percentage %")
println("最大浓度误差: $max_concentration_error mol/m³")
println("最大浓度误差百分比: $max_concentration_error_percentage %")
println("最大应力误差: $max_stress_error Pa")
println("最大应力误差百分比: $max_stress_error_percentage %")

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time_result,
    "Voltage (P2D) (V)" => voltage_p2d,
    "Voltage (sP2D) (V)" => voltage_sP2D,
    "Voltage Error (V)" => voltage_error,
    "Voltage Error Percentage (%)" => voltage_error_percentage,
    "Concentration (P2D) (mol/m³)" => concentration_p2d,
    "Concentration (sP2D) (mol/m³)" => concentration_sP2D,
    "Concentration Error (mol/m³)" => concentration_error,
    "Concentration Error Percentage (%)" => concentration_error_percentage,
    "Stress (P2D) (Pa)" => stress_p2d,
    "Stress (sP2D) (Pa)" => stress_sP2D,
    "Stress Error (Pa)" => stress_error,
    "Stress Error Percentage (%)" => stress_error_percentage,
    "Current (A)" => current
)
CSV.write("Harmonic_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")

# 引用信息
JuBat.Citation()