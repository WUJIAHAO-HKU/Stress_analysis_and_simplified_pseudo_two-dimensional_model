using Plots, CSV, DataFrames, Interpolations
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

# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = collect(Float64, 0:1:maximum(time_udds))  # 以1秒的间隔进行模拟

# 设置电流函数为UDDS表格中的数据
opt.Current = t -> current_interp(t)

# 创建和运行模拟
println("开始模拟UDDS工况... P2D模型")
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
println("P2D模型模拟完成")

# 提取P2D模型结果数据
time = result["time [s]"]
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

# 创建P2D模型的子图
p1 = plot(time, voltage_p2d, 
    label="电压 (P2D)", 
    xlabel="时间 [s]", 
    ylabel="电压 [V]", 
    title="UDDS工况",
    lw=2)

p2 = plot(time, concentration_p2d, 
    label="浓度 (P2D)", 
    xlabel="时间 [s]", 
    ylabel="浓度 [mol/m³]", 
    lw=2)

p3 = plot(time, stress_p2d, 
    label="应力 (P2D)", 
    xlabel="时间 [s]", 
    ylabel="应力 [Pa]", 
    lw=2)

# 创建电流图
p7 = plot(time, [current_interp(t) for t in time], 
label="Current", 
xlabel="time [s]", 
ylabel="Current [A]", 
title="UDDS Current Profile",
lw=2, size=(800, 400), color=:mediumslateblue)

# 组合P2D模型图表
plot_combined_p2d = plot(p1, p2, p3, p7,  
    layout=(4, 1), 
    size=(800, 600))

# 设置sP2D模型
opt.model = "sP2D"  # 使用sP2D模型
println("开始模拟UDDS工况... sP2D模型")
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

# 创建sP2D模型的子图
p4 = plot(time, voltage_sP2D, 
    label="电压 (sP2D)", 
    xlabel="时间 [s]", 
    ylabel="电压 [V]", 
    title="UDDS工况",
    lw=2, linecolor=:red)

p5 = plot(time, concentration_sP2D, 
    label="浓度 (sP2D)", 
    xlabel="时间 [s]", 
    ylabel="浓度 [mol/m³]", 
    lw=2, linecolor=:green)

p6 = plot(time, stress_sP2D, 
    label="应力 (sP2D)", 
    xlabel="时间 [s]", 
    ylabel="应力 [Pa]", 
    lw=2, linecolor=:blue)

# 组合sP2D模型图表
plot_combined_sP2D = plot(p4, p5, p6, p7, 
    layout=(4, 1), 
    size=(800, 600))

# 组合P2D和sP2D模型的结果图表
plot_combined_all = plot(plot_combined_p2d, plot_combined_sP2D, 
    layout=(1, 2), 
    size=(1000, 800))

# 保存图表
savefig(plot_combined_p2d, "udds_battery_analysis_p2d.pdf")
savefig(plot_combined_sP2D, "udds_battery_analysis_sP2D.pdf")
savefig(p7, "udds_current_profile.pdf")
savefig(plot_combined_all, "udds_battery_analysis_comparison.pdf")

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
    "Time (s)" => time,
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
    "Current (A)" => [current_interp(t) for t in time]
)
CSV.write("udds_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")

# 引用信息
JuBat.Citation()