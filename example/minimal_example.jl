using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("Enertech")

opt = JuBat.Option()
opt.mechanicalmodel = "full"
Crate = 1
i= 2.28*Crate
opt.Current = x-> i
opt.time = [0 3600]
opt.model = "P2D" # choose model, other options are "SPM" or "SPMe"
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

# 提取数据
time = result["time [s]"]
voltage = result["cell voltage [V]"]
concentration = result["negative particle surface lithium concentration [mol/m^3]"][1, :]
stress = result["negative particle surface tangential stress[Pa]"][1, :]

# 创建三个子图
p1 = plot(time, voltage, label="voltage", xlabel="time [s]", ylabel="voltage [V]", title="5C discharging")
p2 = plot(time, concentration, label="concentration", xlabel="time [s]", ylabel="concentration [mol/m³]")
p3 = plot(time, stress, label="stress", xlabel="time [s]", ylabel="stress [Pa]")

# 组合图表
plot_combined = plot(p1, p2, p3, layout=(3, 1), size=(800, 600))

# 保存图表
savefig(plot_combined, "battery_analysis_p2d.pdf")

# 引用信息
JuBat.Citation()