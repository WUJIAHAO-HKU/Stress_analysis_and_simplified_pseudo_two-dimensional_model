using Plots, CSV, DataFrames, Statistics, LinearAlgebra
include("../src/JuBat.jl")

# 敏感性分析与参数空间探索
function sensitivity_analysis()
    println("开始敏感性分析...")
    
    # 定义要分析的参数及其变化范围
    param_variations = Dict(
        "扩散系数负极" => [0.5, 0.75, 1.0, 1.25, 1.5, 2.0],  # 相对基准值的倍数
        "扩散系数正极" => [0.5, 0.75, 1.0, 1.25, 1.5, 2.0],
        "反应速率常数负极" => [0.1, 0.5, 1.0, 2.0, 5.0, 10.0],
        "反应速率常数正极" => [0.1, 0.5, 1.0, 2.0, 5.0, 10.0],
        "电极孔隙率" => [0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
        "固相粒子半径" => [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    )
    
    models = ["P2D", "sP2D"]
    sensitivity_results = Dict()
    
    # 基准仿真
    base_results = run_base_simulation()
    
    for param_name in keys(param_variations)
        param_results = Dict()
        
        for model in models
            model_results = []
            
            for multiplier in param_variations[param_name]
                println("分析参数: $param_name, 倍数: $multiplier, 模型: $model")
                
                # 运行参数变化仿真
                result = run_parameter_variation(param_name, multiplier, model)
                
                # 计算敏感性指标
                voltage_sensitivity = calculate_voltage_sensitivity(result, base_results[model])
                capacity_sensitivity = calculate_capacity_sensitivity(result, base_results[model])
                stress_sensitivity = calculate_stress_sensitivity(result, base_results[model])
                
                push!(model_results, Dict(
                    "parameter" => param_name,
                    "multiplier" => multiplier,
                    "model" => model,
                    "voltage_sensitivity" => voltage_sensitivity,
                    "capacity_sensitivity" => capacity_sensitivity,
                    "stress_sensitivity" => stress_sensitivity,
                    "final_voltage" => result["cell voltage [V]"][end],
                    "calc_time" => result["calc_time"]
                ))
            end
            
            param_results[model] = model_results
        end
        
        sensitivity_results[param_name] = param_results
    end
    
    return sensitivity_results
end

# 运行基准仿真
function run_base_simulation()
    base_results = Dict()
    
    for model in ["P2D", "sP2D"]
        param_dim = JuBat.ChooseCell("LG M50")
        opt = JuBat.Option()
        opt.mechanicalmodel = "full"
        opt.model = model
        opt.dtType = "fixed"
        opt.dt = [2.0, 20.0]
        opt.time = [0, 3600]  # 1小时
        opt.Current = t -> 5.0  # 1C放电
        
        calc_time = @elapsed begin
            case = JuBat.SetCase(param_dim, opt)
            result = JuBat.Solve(case)
        end
        
        result["calc_time"] = calc_time
        base_results[model] = result
    end
    
    return base_results
end

# 运行参数变化仿真
function run_parameter_variation(param_name, multiplier, model)
    param_dim = JuBat.ChooseCell("LG M50")
    
    # 根据参数名称修改相应参数
    if param_name == "扩散系数负极"
        param_dim.NE.Ds = param_dim.NE.Ds * multiplier
    elseif param_name == "扩散系数正极"
        param_dim.PE.Ds = param_dim.PE.Ds * multiplier
    elseif param_name == "反应速率常数负极"
        param_dim.NE.k = param_dim.NE.k * multiplier
    elseif param_name == "反应速率常数正极"
        param_dim.PE.k = param_dim.PE.k * multiplier
    elseif param_name == "电极孔隙率"
        param_dim.NE.eps = multiplier
        param_dim.PE.eps = multiplier
    elseif param_name == "固相粒子半径"
        param_dim.NE.rs = param_dim.NE.rs * multiplier
        param_dim.PE.rs = param_dim.PE.rs * multiplier
    end
    
    opt = JuBat.Option()
    opt.mechanicalmodel = "full"
    opt.model = model
    opt.dtType = "fixed"
    opt.dt = [2.0, 20.0]
    opt.time = [0, 3600]
    opt.Current = t -> 5.0
    
    calc_time = @elapsed begin
        case = JuBat.SetCase(param_dim, opt)
        result = JuBat.Solve(case)
    end
    
    result["calc_time"] = calc_time
    return result
end

# 计算电压敏感性
function calculate_voltage_sensitivity(varied_result, base_result)
    varied_voltage = varied_result["cell voltage [V]"]
    base_voltage = base_result["cell voltage [V]"]
    
    # 计算均方根相对差
    min_length = min(length(varied_voltage), length(base_voltage))
    voltage_diff = abs.(varied_voltage[1:min_length] .- base_voltage[1:min_length])
    rms_error = sqrt(mean(voltage_diff.^2))
    
    return rms_error
end

# 计算容量敏感性
function calculate_capacity_sensitivity(varied_result, base_result)
    # 计算容量 (通过电流积分)
    varied_time = varied_result["time [s]"]
    base_time = base_result["time [s]"]
    
    varied_capacity = sum(abs.(varied_result["cell current [A]"]) .* diff([0; varied_time])) / 3600
    base_capacity = sum(abs.(base_result["cell current [A]"]) .* diff([0; base_time])) / 3600
    
    return abs(varied_capacity - base_capacity) / base_capacity
end

# 计算应力敏感性
function calculate_stress_sensitivity(varied_result, base_result)
    # 提取应力数据
    # 检查两种可能的键名
    stress_keys = ["negative electrode stress [Pa]", "negative electrode tangential stress [Pa]"]
    stress_key = first([k for k in keys(varied_result) if any(sk -> occursin(sk, k), stress_keys)], nothing)
    
    if stress_key === nothing
        @warn "找不到应力数据键，使用默认空数组"
        return 0.0
    end
    
    varied_stress = varied_result[stress_key]
    base_stress = base_result[stress_key]
    
    if ndims(varied_stress) > 1
        varied_stress_final = varied_stress[:, end]
        base_stress_final = base_stress[:, end]
    else
        varied_stress_final = [varied_stress[end]]
        base_stress_final = [base_stress[end]]
    end
    
    # 计算相对差
    stress_diff = abs.(varied_stress_final .- base_stress_final)
    max_stress_diff = maximum(stress_diff)
    max_base_stress = maximum(abs.(base_stress_final))
    
    return max_stress_diff / max_base_stress
end

# 创建敏感性分析可视化
function create_sensitivity_plots(sensitivity_results)
    param_names = collect(keys(sensitivity_results))
    
    # 1. 电压敏感性热图
    p1 = create_sensitivity_heatmap(sensitivity_results, "voltage_sensitivity", "Voltage Sensitivity")
    
    # 2. 容量敏感性热图
    p2 = create_sensitivity_heatmap(sensitivity_results, "capacity_sensitivity", "Capacity Sensitivity")
    
    # 3. 应力敏感性热图
    p3 = create_sensitivity_heatmap(sensitivity_results, "stress_sensitivity", "Stress Sensitivity")
    
    # 4. 模型差异敏感性分析
    p4 = create_model_difference_plot(sensitivity_results)
    
    # 组合图
    sensitivity_plot = plot(p1, p2, p3, p4,
                           layout=(2, 2),
                           size=(1400, 1000),
                           dpi=300)
    
    savefig(sensitivity_plot, "sensitivity_analysis.pdf")
    savefig(sensitivity_plot, "sensitivity_analysis.png")
    
    return sensitivity_plot
end

# 创建敏感性热图
function create_sensitivity_heatmap(sensitivity_results, metric, title)
    param_names = collect(keys(sensitivity_results))
    models = ["P2D", "sP2D"]
    
    # 创建数据矩阵
    n_params = length(param_names)
    n_variations = 6  # 每个参数6个变化值
    
    p2d_matrix = zeros(n_params, n_variations)
    sp2d_matrix = zeros(n_params, n_variations)
    
    for (i, param) in enumerate(param_names)
        for (j, result) in enumerate(sensitivity_results[param]["P2D"])
            p2d_matrix[i, j] = result[metric]
        end
        for (j, result) in enumerate(sensitivity_results[param]["sP2D"])
            sp2d_matrix[i, j] = result[metric]
        end
    end
    
    # 创建P2D热图
    p_p2d = heatmap(p2d_matrix,
                    title="$title - P2D Model",
                    xlabel="Parameter Variation Index",
                    ylabel="Parameters",
                    yticks=(1:n_params, param_names),
                    color=:viridis,
                    aspect_ratio=:auto)
    
    # 创建sP2D热图
    p_sp2d = heatmap(sp2d_matrix,
                     title="$title - sP2D Model",
                     xlabel="Parameter Variation Index",
                     ylabel="Parameters",
                     yticks=(1:n_params, param_names),
                     color=:viridis,
                     aspect_ratio=:auto)
    
    # 组合P2D和sP2D
    combined = plot(p_p2d, p_sp2d,
                   layout=(1, 2),
                   size=(800, 400))
    
    return combined
end

# 创建模型差异图
function create_model_difference_plot(sensitivity_results)
    param_names = collect(keys(sensitivity_results))
    
    p = plot(title="Model Sensitivity Differences",
             xlabel="Parameters",
             ylabel="Relative Sensitivity Difference",
             legend=:topright,
             xrotation=45)
    
    voltage_diffs = []
    capacity_diffs = []
    stress_diffs = []
    
    for param in param_names
        # 计算每个参数的平均敏感性差异
        p2d_results = sensitivity_results[param]["P2D"]
        sp2d_results = sensitivity_results[param]["sP2D"]
        
        voltage_diff = mean([abs(sp2d["voltage_sensitivity"] - p2d["voltage_sensitivity"]) 
                           for (p2d, sp2d) in zip(p2d_results, sp2d_results)])
        capacity_diff = mean([abs(sp2d["capacity_sensitivity"] - p2d["capacity_sensitivity"]) 
                            for (p2d, sp2d) in zip(p2d_results, sp2d_results)])
        stress_diff = mean([abs(sp2d["stress_sensitivity"] - p2d["stress_sensitivity"]) 
                          for (p2d, sp2d) in zip(p2d_results, sp2d_results)])
        
        push!(voltage_diffs, voltage_diff)
        push!(capacity_diffs, capacity_diff)
        push!(stress_diffs, stress_diff)
    end
    
    x_positions = 1:length(param_names)
    
    plot!(p, x_positions .- 0.2, voltage_diffs,
          label="Voltage Sensitivity",
          marker=:circle,
          bar=true,
          alpha=0.7,
          color=:blue)
    
    plot!(p, x_positions, capacity_diffs,
          label="Capacity Sensitivity",
          marker=:diamond,
          bar=true,
          alpha=0.7,
          color=:red)
    
    plot!(p, x_positions .+ 0.2, stress_diffs,
          label="Stress Sensitivity",
          marker=:star,
          bar=true,
          alpha=0.7,
          color=:green)
    
    plot!(p, xticks=(x_positions, param_names))
    
    return p
end

# 创建参数空间探索图
function create_parameter_space_plot(sensitivity_results)
    # 选择两个最重要的参数进行2D参数空间分析
    param1 = "扩散系数负极"
    param2 = "反应速率常数负极"
    
    variations1 = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    variations2 = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    
    # 创建2D网格
    n1, n2 = length(variations1), length(variations2)
    voltage_surface_p2d = zeros(n1, n2)
    voltage_surface_sp2d = zeros(n1, n2)
    
    for (i, var1) in enumerate(variations1)
        for (j, var2) in enumerate(variations2)
            # 运行双参数变化仿真
            result_p2d = run_dual_parameter_variation(param1, var1, param2, var2, "P2D")
            result_sp2d = run_dual_parameter_variation(param1, var1, param2, var2, "sP2D")
            
            voltage_surface_p2d[i, j] = result_p2d["cell voltage [V]"][end]
            voltage_surface_sp2d[i, j] = result_sp2d["cell voltage [V]"][end]
        end
    end
    
    # 创建3D表面图
    p1 = surface(variations1, variations2, voltage_surface_p2d',
                title="P2D Model - Parameter Space",
                xlabel=param1,
                ylabel=param2,
                zlabel="Final Voltage [V]")
    
    p2 = surface(variations1, variations2, voltage_surface_sp2d',
                title="sP2D Model - Parameter Space",
                xlabel=param1,
                ylabel=param2,
                zlabel="Final Voltage [V]")
    
    # 创建差异图
    voltage_diff = abs.(voltage_surface_sp2d - voltage_surface_p2d)
    p3 = surface(variations1, variations2, voltage_diff',
                title="Voltage Difference (sP2D - P2D)",
                xlabel=param1,
                ylabel=param2,
                zlabel="Voltage Difference [V]")
    
    space_plot = plot(p1, p2, p3,
                     layout=(1, 3),
                     size=(1500, 500),
                     dpi=300)
    
    savefig(space_plot, "parameter_space_exploration.pdf")
    savefig(space_plot, "parameter_space_exploration.png")
    
    return space_plot
end

# 运行双参数变化仿真
function run_dual_parameter_variation(param1, var1, param2, var2, model)
    param_dim = JuBat.ChooseCell("LG M50")
    
    # 修改第一个参数
    if param1 == "扩散系数负极"
        param_dim.NE.Ds = param_dim.NE.Ds * var1
    elseif param1 == "反应速率常数负极"
        param_dim.NE.k = param_dim.NE.k * var1
    end
    
    # 修改第二个参数
    if param2 == "扩散系数负极"
        param_dim.NE.Ds = param_dim.NE.Ds * var2
    elseif param2 == "反应速率常数负极"
        param_dim.NE.k = param_dim.NE.k * var2
    end
    
    opt = JuBat.Option()
    opt.mechanicalmodel = "full"
    opt.model = model
    opt.dtType = "fixed"
    opt.dt = [2.0, 20.0]
    opt.time = [0, 1800]  # 30分钟（减少计算时间）
    opt.Current = t -> 5.0
    
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
    
    return result
end

# 执行敏感性分析
println("开始敏感性分析与参数空间探索...")
sensitivity_results = sensitivity_analysis()

# 生成可视化
sensitivity_plot = create_sensitivity_plots(sensitivity_results)
parameter_space_plot = create_parameter_space_plot(sensitivity_results)

# 保存数据
sensitivity_data = []
for param in keys(sensitivity_results)
    for model in ["P2D", "sP2D"]
        for result in sensitivity_results[param][model]
            push!(sensitivity_data, result)
        end
    end
end

sensitivity_df = DataFrame(sensitivity_data)
CSV.write("sensitivity_analysis_results.csv", sensitivity_df)

# 统计分析
println("=== 敏感性分析结果统计 ===")
for param in keys(sensitivity_results)
    p2d_voltage_sens = [r["voltage_sensitivity"] for r in sensitivity_results[param]["P2D"]]
    sp2d_voltage_sens = [r["voltage_sensitivity"] for r in sensitivity_results[param]["sP2D"]]
    
    p2d_mean = mean(p2d_voltage_sens)
    sp2d_mean = mean(sp2d_voltage_sens)
    sensitivity_diff = abs(sp2d_mean - p2d_mean) / p2d_mean * 100
    
    println("参数: $param")
    println("  P2D平均敏感性: $(round(p2d_mean, digits=6))")
    println("  sP2D平均敏感性: $(round(sp2d_mean, digits=6))")
    println("  相对差异: $(round(sensitivity_diff, digits=2))%")
end

println("敏感性分析完成!") 