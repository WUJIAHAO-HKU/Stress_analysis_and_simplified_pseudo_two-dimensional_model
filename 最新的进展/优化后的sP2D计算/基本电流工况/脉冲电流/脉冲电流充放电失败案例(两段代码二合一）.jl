using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# Define pulse current parameters
pulse_current = 4.0     # Peak pulse current in Amperes
low_current = 0.5       # Low current in Amperes
pulse_width = 200       # Pulse width in seconds
pulse_interval = 1000   # Pulse interval in seconds
total_time = 3000       # Total simulation time in seconds

# 充放电电流参数
charging_current = 4.0  # 充电电流（安培）
discharging_current = -2.0  # 放电电流（安培）
charge_time = 2000  # 充电时间（秒）
discharge_time = 2000  # 放电时间（秒）
total_time = charge_time + discharge_time  # 总时间

opt.dtType = "fixed"
opt.time = [0, total_time]  # Set total simulation time

# Function to define the pulse current
opt.Current = (t) -> begin
    # Determine if the current is in a pulse period
    if t % (pulse_width + pulse_interval) < pulse_width
        return pulse_current  # Apply pulse current during pulse period
    else
        return low_current  # Apply low current outside pulse period
    end
end

# 设置电流随时间变化：充电和放电阶段
opt.Current = (t) -> begin
    if t < charge_time
        return charging_current  # 充电阶段
    else
        return discharging_current  # 放电阶段
    end
end

# Run simulations for both models
for model in ["sP2D", "P2D"]
    opt.model = model
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)

    # Plotting results
    plot!(result["time [s]"], result["cell voltage [V]"], label=model, xlabel="time [s]", ylabel="cell voltage [V]", lw=2)
end

# Save the plot
savefig("pulse_current_validation.pdf")

# Reference for citation
JuBat.Citation(["ai2024b"])
