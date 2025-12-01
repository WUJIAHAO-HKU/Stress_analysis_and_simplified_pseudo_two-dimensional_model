using Plots
include("src/JuBat.jl") 

println("="^60)
println("Testing sP2D thermal model fix")
println("="^60)

param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_l = 2.5

opt = JuBat.Option()
opt.thermalmodel = "lumped"
opt.dtType = "auto"
opt.jacobi = "update"
opt.Current = x -> 5.0  # 1C discharge
opt.dt = [1, 10]
opt.time = [0, 100]  # 短时间测试
opt.model = "sP2D"

println("\n1. Setting up case...")
case1 = JuBat.SetCase(param_dim, opt)
println("   Case setup complete.")

println("\n2. Solving with sP2D model...")
try
    result = JuBat.Solve(case1)
    println("   ✓ Solution successful!")
    println("\n3. Results:")
    println("   Initial voltage: ", round(result["cell voltage [V]"][1], digits=4), " V")
    println("   Final voltage:   ", round(result["cell voltage [V]"][end], digits=4), " V")
    println("   Initial temp:    ", round(result["temperature [K]"][1], digits=2), " K")
    println("   Final temp:      ", round(result["temperature [K]"][end], digits=2), " K")
    println("\n", "="^60)
    println("✓ All tests passed!")
    println("="^60)
catch e
    println("   ✗ Error occurred:")
    println("   ", e)
    rethrow(e)
end
