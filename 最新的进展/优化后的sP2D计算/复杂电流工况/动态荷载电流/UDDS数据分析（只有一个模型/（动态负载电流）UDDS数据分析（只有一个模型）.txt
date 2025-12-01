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
println("开始模拟UDDS工况...")
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
println("模拟完成")

# 提取结果数据
time = result["time [s]"]
voltage = result["cell voltage [V]"]
concentration = result["negative particle surface lithium concentration [mol/m^3]"]
stress = result["negative particle surface tangential stress[Pa]"]

# 处理结果的维度
if ndims(concentration) > 1
    concentration = concentration[1, :]
end

if ndims(stress) > 1
    stress = stress[1, :]
end

# 创建子图
p1 = plot(time, voltage, 
    label="电压", 
    xlabel="时间 [s]", 
    ylabel="电压 [V]", 
    title="UDDS工况",
    lw=2)

p2 = plot(time, concentration, 
    label="浓度", 
    xlabel="时间 [s]", 
    ylabel="浓度 [mol/m³]", 
    lw=2)

p3 = plot(time, stress, 
    label="应力", 
    xlabel="时间 [s]", 
    ylabel="应力 [Pa]", 
    lw=2)

# 组合图表
plot_combined = plot(p1, p2, p3, 
    layout=(3, 1), 
    size=(800, 600))

# 保存图表
savefig(plot_combined, "udds_battery_analysis.pdf")

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time,
    "Voltage (V)" => voltage,
    "Concentration (mol/m³)" => concentration,
    "Stress (Pa)" => stress,
    "Current (A)" => [current_interp(t) for t in time]
)
CSV.write("udds_battery_results.csv", results_df)

println("分析完成，结果已保存")

# 引用信息
JuBat.Citation()