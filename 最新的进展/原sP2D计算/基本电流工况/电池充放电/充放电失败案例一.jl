using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
rates = [ 0.5, 1, 2]
ts = [3700, 3700, 3600, 3600]

# Loop through rates for discharging and charging
for j = 1:3
    Crate = rates[j]
    opt.dt = [1, 20]/Crate
    opt.dtType = "fixed"
    i = 5*Crate
    opt.Current = x -> i
    opt.time = [0 ts[j]/Crate]
    opt.model = "sP2D" # Discharge model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    opt.model = "P2D"  # Alternative model
    case2 = JuBat.SetCase(param_dim, opt)
    result2 = JuBat.Solve(case2)

    # Plotting discharging results
    plot!(result["time [s]"], result["cell voltage [V]"], label="sP2D Discharge", xlabel="Time [s]", ylabel="Cell Voltage [V]", lw=1)
    plot!(result2["time [s]"], result2["cell voltage [V]"], label="P2D Discharge", linestyle=:dot, linecolor=:black, lw=2.5)
    
    # Adding charging case (positive current)
    opt.Current = x -> -i  # Negative current for charging
    opt.time = [0 ts[j]/Crate]
    case3 = JuBat.SetCase(param_dim, opt)
    result3 = JuBat.Solve(case3)
    
    opt.model = "P2D"  # Alternative charging model
    case4 = JuBat.SetCase(param_dim, opt)
    result4 = JuBat.Solve(case4)

    # Plotting charging results
    plot!(result3["time [s]"], result3["cell voltage [V]"], label="sP2D Charge", linestyle=:dash, linecolor=:blue, lw=1.5)
    plot!(result4["time [s]"], result4["cell voltage [V]"], label="P2D Charge", linestyle=:dashdot, linecolor=:red, lw=2)

end

savefig("sP2D_P2D_validation_en.pdf")  # Save the updated plot with English labels
JuBat.Citation(["ai2024b"])
