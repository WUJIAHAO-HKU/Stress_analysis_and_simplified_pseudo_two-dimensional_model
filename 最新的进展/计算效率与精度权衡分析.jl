using Plots, CSV, DataFrames, Statistics, BenchmarkTools
include("../src/JuBat.jl")

# 计算效率与精度权衡分析
function efficiency_accuracy_tradeoff_analysis()
    println("开始计算效率与精度权衡分析...")
    
    # 定义不同的仿真场景
    scenarios = Dict(
        "快速充电" => Dict(
            "current_func" => t -> t < 1800 ? 10.0 : 0.0,  # 0.5小时10A快充
            "time_range" => [0, 3600],
            "description" => "Fast Charging (2C)"
        ),
        "标准放电" => Dict(
            "current_func" => t -> 5.0,  # 1C恒流放电
            "time_range" => [0, 3600],
            "description" => "Standard Discharge (1C)"
        ),
        "脉冲功率" => Dict(
            "current_func" => t -> 15.0 * (mod(t, 120) < 30 ? 1 : 0),  # 脉冲电流
            "time_range" => [0, 1800],
            "description" => "Pulse Power Test"
        ),
        "动态负载" => Dict(
            "current_func" => t -> 3.0 + 2.0 * sin(2π * t / 300),  # 动态电流
            "time_range" => [0, 2400],
            "description" => "Dynamic Load Profile"
        ),
        "多阶梯" => Dict(
            "current_func" => t -> t < 600 ? 2.0 : (t < 1200 ? 5.0 : (t < 1800 ? 8.0 : 3.0)),
            "time_range" => [0, 2400],
            "description" => "Multi-step Current"
        )
    )
    
    # 定义不同的时间步长设置（精度级别）
    precision_levels = Dict(
        "粗糙" => Dict("dt" => [5.0, 50.0], "description" => "Coarse (Fast)"),
        "中等" => Dict("dt" => [2.0, 20.0], "description" => "Medium (Balanced)"),
        "精细" => Dict("dt" => [1.0, 10.0], "description" => "Fine (Accurate)"),
        "极精细" => Dict("dt" => [0.5, 5.0], "description" => "Very Fine (High Accuracy)")
    )
    
    models = ["P2D", "sP2D"]
    results = Dict()
    
    # 首先运行基准仿真（极精细P2D）
    reference_results = run_reference_simulations(scenarios)
    
    for scenario_name in keys(scenarios)
        println("分析场景: $scenario_name")
        scenario_results = Dict()
        
        for precision_name in keys(precision_levels)
            println("  精度级别: $precision_name")
            precision_results = Dict()
            
            for model in models
                println("    模型: $model")
                
                # 运行仿真
                result = run_efficiency_simulation(
                    scenarios[scenario_name], 
                    precision_levels[precision_name], 
                    model
                )
                
                # 计算精度指标
                accuracy_metrics = calculate_accuracy_metrics(
                    result, 
                    reference_results[scenario_name],
                    scenario_name
                )
                
                # 合并结果
                combined_result = merge(result, accuracy_metrics)
                combined_result["scenario"] = scenario_name
                combined_result["precision"] = precision_name
                combined_result["model"] = model
                
                precision_results[model] = combined_result
            end
            
            scenario_results[precision_name] = precision_results
        end
        
        results[scenario_name] = scenario_results
    end
    
    return results
end

# 运行基准仿真
function run_reference_simulations(scenarios)
    println("运行基准仿真（极精细P2D）...")
    reference_results = Dict()
    
    for (scenario_name, scenario_config) in scenarios
        println("  基准场景: $scenario_name")
        
        param_dim = JuBat.ChooseCell("LG M50")
        opt = JuBat.Option()
        opt.mechanicalmodel = "full"
        opt.model = "P2D"
        opt.dtType = "fixed"
        opt.dt = [0.1, 1.0]  # 极精细时间步长
        opt.time = scenario_config["time_range"]
        opt.Current = scenario_config["current_func"]
        
        case = JuBat.SetCase(param_dim, opt)
        result = JuBat.Solve(case)
        
        reference_results[scenario_name] = result
    end
    
    return reference_results
end

# 运行效率仿真
function run_efficiency_simulation(scenario_config, precision_config, model)
    param_dim = JuBat.ChooseCell("LG M50")
    opt = JuBat.Option()
    opt.mechanicalmodel = "full"
    opt.model = model
    opt.dtType = "fixed"
    opt.dt = precision_config["dt"]
    opt.time = scenario_config["time_range"]
    opt.Current = scenario_config["current_func"]
    
    # 使用BenchmarkTools精确测量计算时间
    benchmark_result = @benchmark begin
        local case = JuBat.SetCase($param_dim, $opt)
        local result = JuBat.Solve(case)
        result
    end samples=3 evals=1
    
    # 获取仿真结果
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    
    # 提取计算性能指标
    calc_time = minimum(benchmark_result.times) / 1e9  # 转换为秒
    memory_usage = benchmark_result.memory / 1024^2    # 转换为MB
    
    return Dict(
        "simulation_result" => result,
        "calc_time" => calc_time,
        "memory_usage" => memory_usage,
        "time_data" => result["time [s]"],
        "voltage_data" => result["cell voltage [V]"],
        "current_data" => result["cell current [A]"]
    )
end

# 计算精度指标
function calculate_accuracy_metrics(test_result, reference_result, scenario_name)
    test_voltage = test_result["simulation_result"]["cell voltage [V]"]
    ref_voltage = reference_result["cell voltage [V]"]
    
    test_time = test_result["simulation_result"]["time [s]"]
    ref_time = reference_result["time [s]"]
    
    # 插值到相同时间网格
    common_time = range(0, stop=minimum([maximum(test_time), maximum(ref_time)]), length=1000)
    
    # 使用线性插值
    test_voltage_interp = linear_interpolation(test_time, test_voltage, common_time)
    ref_voltage_interp = linear_interpolation(ref_time, ref_voltage, common_time)
    
    # 计算误差指标
    voltage_diff = test_voltage_interp .- ref_voltage_interp
    
    # 均方根误差 (RMSE)
    rmse = sqrt(mean(voltage_diff.^2))
    
    # 平均绝对误差 (MAE)
    mae = mean(abs.(voltage_diff))
    
    # 最大绝对误差
    max_error = maximum(abs.(voltage_diff))
    
    # 相对误差
    mean_ref_voltage = mean(ref_voltage_interp)
    relative_rmse = rmse / mean_ref_voltage * 100
    
    # R²决定系数
    ss_res = sum(voltage_diff.^2)
    ss_tot = sum((ref_voltage_interp .- mean(ref_voltage_interp)).^2)
    r_squared = 1 - ss_res / ss_tot
    
    # 容量误差
    test_capacity = calculate_capacity(test_result["simulation_result"])
    ref_capacity = calculate_capacity(reference_result)
    capacity_error = abs(test_capacity - ref_capacity) / ref_capacity * 100
    
    return Dict(
        "rmse" => rmse,
        "mae" => mae,
        "max_error" => max_error,
        "relative_rmse" => relative_rmse,
        "r_squared" => r_squared,
        "capacity_error" => capacity_error,
        "test_capacity" => test_capacity,
        "ref_capacity" => ref_capacity
    )
end

# 线性插值函数
function linear_interpolation(x_original, y_original, x_new)
    y_new = zeros(length(x_new))
    for (i, x_val) in enumerate(x_new)
        if x_val <= x_original[1]
            y_new[i] = y_original[1]
        elseif x_val >= x_original[end]
            y_new[i] = y_original[end]
        else
            # 找到插值点
            idx = findfirst(x -> x >= x_val, x_original)
            if idx > 1
                x1, x2 = x_original[idx-1], x_original[idx]
                y1, y2 = y_original[idx-1], y_original[idx]
                y_new[i] = y1 + (y2 - y1) * (x_val - x1) / (x2 - x1)
            else
                y_new[i] = y_original[idx]
            end
        end
    end
    return y_new
end

# 计算容量
function calculate_capacity(result)
    time_data = result["time [s]"]
    current_data = result["cell current [A]"]
    return sum(abs.(current_data) .* diff([0; time_data])) / 3600
end

# 创建效率精度权衡可视化
function create_efficiency_accuracy_plots(results)
    # 1. 帕累托前沿图
    p1 = create_pareto_frontier_plot(results)
    
    # 2. 精度vs时间散点图
    p2 = create_accuracy_time_scatter(results)
    
    # 3. 场景对比雷达图
    p3 = create_scenario_radar_plot(results)
    
    # 4. 效率提升分析
    p4 = create_efficiency_improvement_plot(results)
    
    # 组合图
    efficiency_plot = plot(p1, p2, p3, p4,
                          layout=(2, 2),
                          size=(1400, 1000),
                          dpi=300)
    
    savefig(efficiency_plot, "efficiency_accuracy_tradeoff.pdf")
    savefig(efficiency_plot, "efficiency_accuracy_tradeoff.png")
    
    return efficiency_plot
end

# 创建帕累托前沿图
function create_pareto_frontier_plot(results)
    p = plot(title="Pareto Frontier: Efficiency vs Accuracy",
             xlabel="Calculation Time [s]",
             ylabel="Accuracy (R²)",
             xscale=:log10,
             legend=:bottomright,
             grid=true)
    
    colors = Dict("P2D" => :blue, "sP2D" => :red)
    markers = Dict("P2D" => :circle, "sP2D" => :diamond)
    
    for model in ["P2D", "sP2D"]
        calc_times = []
        r_squared_values = []
        labels = []
        
        for scenario_name in keys(results)
            for precision_name in keys(results[scenario_name])
                result = results[scenario_name][precision_name][model]
                push!(calc_times, result["calc_time"])
                push!(r_squared_values, result["r_squared"])
                push!(labels, "$(scenario_name[1:2])-$(precision_name[1])")
            end
        end
        
        scatter!(p, calc_times, r_squared_values,
                label="$model Model",
                marker=markers[model],
                color=colors[model],
                markersize=6,
                alpha=0.7)
    end
    
    # 添加理想区域
    plot!(p, xlims=(0.1, 1000), ylims=(0.8, 1.0))
    
    return p
end

# 创建精度时间散点图
function create_accuracy_time_scatter(results)
    p = plot(title="Accuracy vs Computation Time by Scenario",
             xlabel="RMSE [V]",
             ylabel="Calculation Time [s]",
             yscale=:log10,
             legend=:topright,
             grid=true)
    
    scenario_colors = [:blue, :red, :green, :orange, :purple]
    
    for (i, scenario_name) in enumerate(keys(results))
        for model in ["P2D", "sP2D"]
            rmse_values = []
            calc_times = []
            
            for precision_name in keys(results[scenario_name])
                result = results[scenario_name][precision_name][model]
                push!(rmse_values, result["rmse"])
                push!(calc_times, result["calc_time"])
            end
            
            linestyle = model == "P2D" ? :solid : :dash
            
            plot!(p, rmse_values, calc_times,
                  label="$scenario_name ($model)",
                  color=scenario_colors[i],
                  linestyle=linestyle,
                  marker=:circle,
                  linewidth=2,
                  markersize=4)
        end
    end
    
    return p
end

# 创建场景雷达图
function create_scenario_radar_plot(results)
    # 计算各场景下的平均性能指标
    scenarios = collect(keys(results))
    n_scenarios = length(scenarios)
    
    # 为P2D和sP2D分别计算平均加速比
    speedup_ratios = []
    accuracy_scores = []
    
    for scenario in scenarios
        # 计算加速比（相对于最精细的P2D）
        p2d_fine_time = results[scenario]["极精细"]["P2D"]["calc_time"]
        sp2d_medium_time = results[scenario]["中等"]["sP2D"]["calc_time"]
        speedup = p2d_fine_time / sp2d_medium_time
        
        # 计算精度分数（基于R²）
        accuracy = results[scenario]["中等"]["sP2D"]["r_squared"]
        
        push!(speedup_ratios, speedup)
        push!(accuracy_scores, accuracy)
    end
    
    # 标准化数据到[0,1]
    norm_speedup = (speedup_ratios .- minimum(speedup_ratios)) ./ (maximum(speedup_ratios) - minimum(speedup_ratios))
    norm_accuracy = (accuracy_scores .- minimum(accuracy_scores)) ./ (maximum(accuracy_scores) - minimum(accuracy_scores))
    
    angles = range(0, 2π, length=n_scenarios+1)[1:end-1]
    
    p = plot(title="Performance Radar Chart",
             projection=:polar,
             legend=:topright)
    
    # sP2D性能
    plot!(p, angles, norm_speedup,
          label="Speedup (sP2D vs P2D)",
          linewidth=2,
          marker=:circle,
          fill=true,
          alpha=0.3)
    
    plot!(p, angles, norm_accuracy,
          label="Accuracy (R²)",
          linewidth=2,
          marker=:diamond,
          fill=true,
          alpha=0.3)
    
    return p
end

# 创建效率提升分析图
function create_efficiency_improvement_plot(results)
    p = plot(title="Computational Efficiency Improvement",
             xlabel="Scenarios",
             ylabel="Speedup Factor (P2D/sP2D)",
             legend=:topright,
             xrotation=45)
    
    scenarios = collect(keys(results))
    precision_levels = ["粗糙", "中等", "精细", "极精细"]
    colors = [:lightblue, :blue, :darkblue, :navy]
    
    for (i, precision) in enumerate(precision_levels)
        speedups = []
        
        for scenario in scenarios
            p2d_time = results[scenario][precision]["P2D"]["calc_time"]
            sp2d_time = results[scenario][precision]["sP2D"]["calc_time"]
            speedup = p2d_time / sp2d_time
            push!(speedups, speedup)
        end
        
        x_positions = (1:length(scenarios)) .+ (i-2.5)*0.15
        
        plot!(p, x_positions, speedups,
              label="$precision Precision",
              marker=:circle,
              linewidth=2,
              markersize=6,
              color=colors[i])
    end
    
    plot!(p, xticks=(1:length(scenarios), scenarios))
    hline!(p, [1], linestyle=:dash, color=:black, label="No Improvement")
    
    return p
end

# 创建详细性能分析表
function create_performance_table(results)
    table_data = []
    
    for scenario_name in keys(results)
        for precision_name in keys(results[scenario_name])
            for model in ["P2D", "sP2D"]
                result = results[scenario_name][precision_name][model]
                
                push!(table_data, Dict(
                    "Scenario" => scenario_name,
                    "Precision" => precision_name,
                    "Model" => model,
                    "Calc_Time_s" => round(result["calc_time"], digits=3),
                    "Memory_MB" => round(result["memory_usage"], digits=1),
                    "RMSE_V" => round(result["rmse"], digits=6),
                    "R_Squared" => round(result["r_squared"], digits=4),
                    "Capacity_Error_%" => round(result["capacity_error"], digits=2)
                ))
            end
        end
    end
    
    return DataFrame(table_data)
end

# 执行效率精度权衡分析
println("开始计算效率与精度权衡分析...")
efficiency_results = efficiency_accuracy_tradeoff_analysis()

# 生成可视化
efficiency_plot = create_efficiency_accuracy_plots(efficiency_results)

# 创建性能表
performance_table = create_performance_table(efficiency_results)
CSV.write("efficiency_accuracy_performance_table.csv", performance_table)

# 打印关键发现
println("=== 效率精度权衡分析结果 ===")

# 计算平均加速比
total_speedups = []
for scenario_name in keys(efficiency_results)
    for precision_name in keys(efficiency_results[scenario_name])
        p2d_time = efficiency_results[scenario_name][precision_name]["P2D"]["calc_time"]
        sp2d_time = efficiency_results[scenario_name][precision_name]["sP2D"]["calc_time"]
        speedup = p2d_time / sp2d_time
        push!(total_speedups, speedup)
    end
end

avg_speedup = mean(total_speedups)
max_speedup = maximum(total_speedups)
min_speedup = minimum(total_speedups)

println("平均加速比: $(round(avg_speedup, digits=2))x")
println("最大加速比: $(round(max_speedup, digits=2))x")
println("最小加速比: $(round(min_speedup, digits=2))x")

# 计算平均精度损失
total_accuracy_loss = []
for scenario_name in keys(efficiency_results)
    for precision_name in keys(efficiency_results[scenario_name])
        p2d_r2 = efficiency_results[scenario_name][precision_name]["P2D"]["r_squared"]
        sp2d_r2 = efficiency_results[scenario_name][precision_name]["sP2D"]["r_squared"]
        accuracy_loss = (p2d_r2 - sp2d_r2) * 100
        push!(total_accuracy_loss, accuracy_loss)
    end
end

avg_accuracy_loss = mean(total_accuracy_loss)
println("平均精度损失: $(round(avg_accuracy_loss, digits=2))%")

println("效率精度权衡分析完成!") 