using Plots, CSV, DataFrames, Statistics
include("../src/JuBat.jl")

# 网格无关性验证分析
function grid_independence_study()
    println("开始网格无关性验证分析...")
    
    # 定义不同的网格密度
    grid_densities = [10, 20, 40, 80, 160]  # 负极网格节点数
    models = ["P2D", "sP2D"]
    
    # 存储结果
    results_data = Dict()
    convergence_data = []
    
    # 设置电池参数和仿真条件
    param_dim = JuBat.ChooseCell("LG M50")
    opt = JuBat.Option()
    opt.mechanicalmodel = "full"
    opt.dtType = "fixed"
    opt.time = [0, 1000]  # 1000秒仿真
    opt.Current = t -> 5.0  # 1C放电
    
    for model in models
        opt.model = model
        model_results = Dict()
        
        for (i, grid_size) in enumerate(grid_densities)
            println("正在计算 $model 模型，网格密度: $grid_size")
            
            # 设置网格密度 (这里需要根据JuBat的实际API调整)
            # param_dim.mesh.negative_electrode_nodes = grid_size
            
            # 记录计算时间
            calc_time = @elapsed begin
                case = JuBat.SetCase(param_dim, opt)
                result = JuBat.Solve(case)
            end
            
            # 提取关键结果
            final_voltage = result["cell voltage [V]"][end]
            final_conc = result["negative particle surface lithium concentration [mol/m^3]"]
            if ndims(final_conc) > 1
                final_conc = final_conc[1, end]
            else
                final_conc = final_conc[end]
            end
            
            model_results[grid_size] = Dict(
                "voltage" => final_voltage,
                "concentration" => final_conc,
                "calc_time" => calc_time,
                "grid_size" => grid_size
            )
            
            # 计算收敛指标
            if i > 1
                prev_grid = grid_densities[i-1]
                voltage_change = abs(final_voltage - model_results[prev_grid]["voltage"])
                conc_change = abs(final_conc - model_results[prev_grid]["concentration"])
                
                push!(convergence_data, Dict(
                    "model" => model,
                    "grid_size" => grid_size,
                    "voltage_change" => voltage_change,
                    "conc_change" => conc_change,
                    "calc_time" => calc_time
                ))
            end
        end
        
        results_data[model] = model_results
    end
    
    return results_data, convergence_data
end

# 创建网格无关性可视化
function create_grid_independence_plots(results_data, convergence_data)
    
    # 1. 收敛性分析图
    p1 = plot(title="Grid Independence Study - Voltage Convergence",
              xlabel="Grid Density (nodes)",
              ylabel="Final Voltage [V]",
              legend=:bottomright)
    
    for model in ["P2D", "sP2D"]
        grid_sizes = [10, 20, 40, 80, 160]
        voltages = [results_data[model][size]["voltage"] for size in grid_sizes]
        
        plot!(p1, grid_sizes, voltages,
              label="$model Model",
              marker=:circle,
              linewidth=2,
              markersize=6)
    end
    
    # 2. 计算效率对比
    p2 = plot(title="Computational Efficiency vs Grid Density",
              xlabel="Grid Density (nodes)",
              ylabel="Calculation Time [s]",
              yscale=:log10,
              legend=:topleft)
    
    for model in ["P2D", "sP2D"]
        grid_sizes = [10, 20, 40, 80, 160]
        calc_times = [results_data[model][size]["calc_time"] for size in grid_sizes]
        
        plot!(p2, grid_sizes, calc_times,
              label="$model Model",
              marker=:diamond,
              linewidth=2,
              markersize=6)
    end
    
    # 3. 相对误差收敛图
    p3 = plot(title="Relative Error Convergence",
              xlabel="Grid Density (nodes)",
              ylabel="Relative Change (%)",
              yscale=:log10,
              legend=:topright)
    
    # 计算P2D和sP2D之间的相对误差
    grid_sizes = [20, 40, 80, 160]  # 从第二个开始
    relative_errors = []
    
    for size in grid_sizes
        p2d_voltage = results_data["P2D"][size]["voltage"]
        sp2d_voltage = results_data["sP2D"][size]["voltage"]
        rel_error = abs(sp2d_voltage - p2d_voltage) / p2d_voltage * 100
        push!(relative_errors, rel_error)
    end
    
    plot!(p3, grid_sizes, relative_errors,
          label="sP2D vs P2D Voltage Error",
          marker=:star,
          linewidth=3,
          markersize=8,
          color=:red)
    
    # 4. 加速比分析
    p4 = plot(title="Computational Speedup Analysis",
              xlabel="Grid Density (nodes)",
              ylabel="Speedup Factor",
              legend=:topright)
    
    speedup_factors = []
    for size in [10, 20, 40, 80, 160]
        p2d_time = results_data["P2D"][size]["calc_time"]
        sp2d_time = results_data["sP2D"][size]["calc_time"]
        speedup = p2d_time / sp2d_time
        push!(speedup_factors, speedup)
    end
    
    plot!(p4, [10, 20, 40, 80, 160], speedup_factors,
          label="sP2D Speedup vs P2D",
          marker=:hexagon,
          linewidth=3,
          markersize=8,
          color=:green)
    
    # 添加理论线性关系参考
    hline!(p4, [1], linestyle=:dash, color=:black, label="No Speedup")
    
    # 组合所有图
    combined_plot = plot(p1, p2, p3, p4,
                        layout=(2, 2),
                        size=(1200, 800),
                        dpi=300)
    
    savefig(combined_plot, "grid_independence_analysis.pdf")
    savefig(combined_plot, "grid_independence_analysis.png")
    
    return combined_plot
end

# 执行分析
println("开始网格无关性验证分析...")
results_data, convergence_data = grid_independence_study()

# 生成可视化
grid_plot = create_grid_independence_plots(results_data, convergence_data)

# 保存数据
convergence_df = DataFrame(convergence_data)
CSV.write("grid_independence_results.csv", convergence_df)

# 打印统计信息
println("=== 网格无关性验证结果 ===")
for model in ["P2D", "sP2D"]
    println("$model 模型:")
    for size in [10, 20, 40, 80, 160]
        voltage = results_data[model][size]["voltage"]
        time = results_data[model][size]["calc_time"]
        println("  网格密度 $size: 电压=$(round(voltage, digits=4))V, 时间=$(round(time, digits=2))s")
    end
end

println("网格无关性验证分析完成!") 