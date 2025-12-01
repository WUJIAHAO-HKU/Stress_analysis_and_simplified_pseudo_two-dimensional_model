using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# Define the rates, time steps for discharging and charging
rates = [0.5, 1, 2]
rates2 = [2.5, 5, 10]
ts_discharge = [3700, 3700, 3600, 3600]
ts_charge = [3700, 3700, 3600, 3600]  # Time for charging

# For discharging simulation
for j = 1:3
    Crate = rates[j]
    opt.dt = [1, 20]/Crate
    opt.dtType = "fixed"
    i = 5 * Crate
    opt.Current = x -> i  # Current for discharging
    opt.time = [0 ts_discharge[j]/Crate]
    opt.model = "sP2D"
    case = JuBat.SetCase(param_dim, opt)
    result_discharge = JuBat.Solve(case)

    opt.model = "P2D"
    case2 = JuBat.SetCase(param_dim, opt)
    result2_discharge = JuBat.Solve(case2)

    # Save the final state after discharge
    final_voltage = result_discharge["cell voltage [V]"][-1]

    # Plotting discharging results
    plot!(result_discharge["time [s]"], result_discharge["cell voltage [V]"], label="sP2D Discharge", xlabel="time [s]", ylabel="cell voltage [V]", lw=1)
    plot!(result2_discharge["time [s]"], result2_discharge["cell voltage [V]"], label="P2D Discharge", linestyle=:dot, linecolor=:black, lw=2.5)

    # For charging simulation, reset the battery state to the final voltage after discharge
    for charge_rate in rates2
        opt.dt = [1, 20]/Crate
        opt.dtType = "fixed"
        opt.Current = x -> charge_rate  # Set current for charging
        opt.time = [0 ts_charge[j]/Crate]
        
        # Initialize the charge simulation with the final voltage state from discharge
        opt.initial_voltage = final_voltage  # Assuming this field is used to set the initial voltage
        opt.model = "sP2D"
        case_charge = JuBat.SetCase(param_dim, opt)
        result_charge = JuBat.Solve(case_charge)

        opt.model = "P2D"
        case2_charge = JuBat.SetCase(param_dim, opt)
        result2_charge = JuBat.Solve(case2_charge)

        # Plotting charging results
        plot!(result_charge["time [s]"], result_charge["cell voltage [V]"], label="sP2D Charge", linestyle=:dash, lw=2.5)
        plot!(result2_charge["time [s]"], result2_charge["cell voltage [V]"], label="P2D Charge", linestyle=:dashdot, linecolor=:blue, lw=2)
    end
end

# Save the figure as a PDF
savefig("sP2D_and_P2D_voltage_comparison.pdf")
JuBat.Citation(["ai2024b"])