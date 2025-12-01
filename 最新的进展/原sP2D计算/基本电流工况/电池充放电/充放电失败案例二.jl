using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
rates = [0.5, 1, 2]
ts = [3700, 3700, 3600, 3600]

# 为放电和充电分别绘制电压曲线
pV_discharge = plot(xlabel="time [s]", ylabel="cell voltage [V]", title="Discharge vs Charge", lw=1)
pV_charge = plot(xlabel="time [s]", ylabel="cell voltage [V]", title="Discharge vs Charge", lw=1)

for j = 1:3
    Crate = rates[j]
    opt.dt = [1, 20] / Crate
    opt.dtType = "fixed"
    
    # 放电仿真
    i = 5 * Crate
    opt.Current = x -> i  # 负电流表示放电
    opt.time = [0 ts[j] / Crate]
    opt.model = "sP2D"
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # 放电仿真 P2D
    opt.model = "P2D"
    case2 = JuBat.SetCase(param_dim, opt)
    result2 = JuBat.Solve(case2)

    # 绘制放电结果
    plot!(pV_discharge, result["time [s]"], result["cell voltage [V]"], label="sP2D (discharge)", lw=1)
    plot!(pV_discharge, result2["time [s]"], result2["cell voltage [V]"], label="P2D (discharge)", linestyle=:dot, linecolor=:black, lw=2.5)

    # 充电仿真
    opt.Current = x -> -i  # 正电流表示充电
    opt.time = [0 ts[j] / Crate]
    case_charge = JuBat.SetCase(param_dim, opt)
    result_charge = JuBat.Solve(case_charge)

    # 充电仿真 P2D
    opt.model = "P2D"
    case_charge2 = JuBat.SetCase(param_dim, opt)
    result_charge2 = JuBat.Solve(case_charge2)

    # 绘制充电结果
    plot!(pV_charge, result_charge["time [s]"], result_charge["cell voltage [V]"], label="sP2D (charge)", lw=1)
    plot!(pV_charge, result_charge2["time [s]"], result_charge2["cell voltage [V]"], label="P2D (charge)", linestyle=:dot, linecolor=:black, lw=2.5)
end

# 保存放电和充电的比较图
savefig(pV_discharge, "discharge_comparison.pdf")
savefig(pV_charge,"charge_comparison.pdf")

JuBat.Citation(["ai2024b"])