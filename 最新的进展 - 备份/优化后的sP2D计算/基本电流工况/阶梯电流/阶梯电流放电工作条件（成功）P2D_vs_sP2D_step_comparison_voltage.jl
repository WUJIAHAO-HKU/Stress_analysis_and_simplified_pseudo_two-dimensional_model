using Plots, CSV, DataFrames
include("../src/JuBat.jl") 

param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()

# Define current steps and corresponding times
current_steps = [0.5, 1.0, 2.0, 1.5]  # Current steps (A)
time_steps = [0, 1000, 2000, 3000]    # Time steps (s) when current changes

opt.dtType = "fixed"
opt.time = [0, time_steps[end]]  # Set total simulation time

# Function to define the current as a step profile
opt.Current = (t) -> begin
    idx = findfirst(x -> t < x, time_steps)
    idx = isnothing(idx) ? length(current_steps) : idx - 1
    return current_steps[idx]
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
savefig("step_current_validation.pdf")

# Reference for citation
JuBat.Citation(["ai2024b"])
