using Plots, CSV, DataFrames, Statistics, Distributions, Random
include("../src/JuBat.jl")

# 不确定性量化与鲁棒性分析
function uncertainty_quantification_analysis()
    println("开始不确定性量化与鲁棒性分析...")
    
    # 设置随机种子确保可重复性
    Random.seed!(1234)
    
    # 定义参数的不确定性分布
    parameter_distributions = Dict(
        "扩散系数负极" => Normal(1.0, 0.2),      # 正态分布，均值1，标准差0.2
        "扩散系数正极" => Normal(1.0, 0.15),
        "反应速率常数负极" => LogNormal(0.0, 0.3),  # 对数正态分布
        "反应速率常数正极" => LogNormal(0.0, 0.25),
        "电极孔隙率" => Beta(5, 5),              # Beta分布，映射到[0.3, 0.8]
        "固相粒子半径" => Normal(1.0, 0.1),
        "电解液电导率" => Normal(1.0, 0.1),
        "温度" => Normal(298.15, 5.0)            # 温度变化±5K
    )
    
    # 蒙特卡洛仿真参数
    n_samples = 100  # 蒙特卡洛样本数
    models = ["P2D", "sP2D"]
    
    # 存储结果
    mc_results = Dict()
    
    for model in models
        println("正在进行 $model 模型的蒙特卡洛分析...")
        model_samples = []
        
        for i in 1:n_samples
            println("进度: $i/$n_samples")
            
            # 生成参数样本
            param_sample = generate_parameter_sample(parameter_distributions)
            
            # 运行仿真
            try
                result = run_uncertain_simulation(param_sample, model)
                
                # 提取关键性能指标
                performance_metrics = extract_performance_metrics(result, param_sample)
                performance_metrics["model"] = model
                performance_metrics["sample_id"] = i
                
                push!(model_samples, performance_metrics)
            catch e
                println("样本 $i 仿真失败: $e")
            end
        end
        
        mc_results[model] = model_samples
    end
    
    return mc_results
end

# 生成参数样本
function generate_parameter_sample(distributions)
    sample = Dict()
    
    for (param_name, dist) in distributions
        if param_name == "电极孔隙率"
            # Beta分布映射到[0.3, 0.8]
            beta_sample = rand(dist)
            sample[param_name] = 0.3 + (0.8 - 0.3) * beta_sample
        elseif param_name == "温度"
            sample[param_name] = rand(dist)
        else
            sample[param_name] = rand(dist)
        end
    end
    
    return sample
end

# 运行不确定性仿真
function run_uncertain_simulation(param_sample, model)
    param_dim = JuBat.ChooseCell("LG M50")
    
    # 应用参数样本
    param_dim.negative_electrode.Ds *= param_sample["扩散系数负极"]
    param_dim.positive_electrode.Ds *= param_sample["扩散系数正极"]
    param_dim.negative_electrode.k0 *= param_sample["反应速率常数负极"]
    param_dim.positive_electrode.k0 *= param_sample["反应速率常数正极"]
    param_dim.negative_electrode.eps_s = param_sample["电极孔隙率"]
    param_dim.positive_electrode.eps_s = param_sample["电极孔隙率"]
    param_dim.negative_electrode.Rs *= param_sample["固相粒子半径"]
    param_dim.positive_electrode.Rs *= param_sample["固相粒子半径"]
    
    # 设置仿真选项
    opt = JuBat.Option()
    opt.mechanicalmodel = "full"
    opt.model = model
    opt.dtType = "fixed"
    opt.dt = [2.0, 20.0]
    opt.time = [0, 3600]  # 1小时
    opt.Current = t -> 5.0  # 1C放电
    
    # 温度设置
    opt.T_amb = param_sample["温度"]
    
    calc_time = @elapsed begin
        case = JuBat.SetCase(param_dim, opt)
        result = JuBat.Solve(case)
    end
    
    result["calc_time"] = calc_time
    return result
end

# 提取性能指标
function extract_performance_metrics(result, param_sample)
    time_data = result["time [s]"]
    voltage_data = result["cell voltage [V]"]
    current_data = result["cell current [A]"]
    
    # 计算容量
    capacity = sum(abs.(current_data) .* diff([0; time_data])) / 3600
    
    # 计算能量
    energy = sum(voltage_data[1:end-1] .* abs.(current_data[1:end-1]) .* diff(time_data)) / 3600
    
    # 电压统计
    final_voltage = voltage_data[end]
    min_voltage = minimum(voltage_data)
    voltage_drop = voltage_data[1] - final_voltage
    voltage_std = std(voltage_data)
    
    # 电压平台长度 (电压在3.0V以上的时间)
    platform_time = sum(time_data[voltage_data .> 3.0]) / maximum(time_data)
    
    # 应力分析
    stress_data = result["negative electrode stress [Pa]"]
    if ndims(stress_data) > 1
        max_stress = maximum(abs.(stress_data[:, end]))
        stress_std = std(stress_data[:, end])
    else
        max_stress = abs(stress_data[end])
        stress_std = 0.0
    end
    
    metrics = Dict(
        "capacity" => capacity,
        "energy" => energy,
        "final_voltage" => final_voltage,
        "min_voltage" => min_voltage,
        "voltage_drop" => voltage_drop,
        "voltage_std" => voltage_std,
        "platform_time" => platform_time,
        "max_stress" => max_stress,
        "stress_std" => stress_std,
        "calc_time" => result["calc_time"]
    )
    
    # 添加参数信息
    for (param, value) in param_sample
        metrics["param_$(param)"] = value
    end
    
    return metrics
end

# 创建不确定性分析可视化
function create_uncertainty_plots(mc_results)
    # 1. 性能指标的不确定性分布
    p1 = create_performance_distribution_plot(mc_results)
    
    # 2. 模型鲁棒性对比
    p2 = create_robustness_comparison_plot(mc_results)
    
    # 3. 参数敏感性雷达图
    p3 = create_sensitivity_radar_plot(mc_results)
    
    # 4. 不确定性传播分析
    p4 = create_uncertainty_propagation_plot(mc_results)
    
    # 组合图
    uncertainty_plot = plot(p1, p2, p3, p4,
                           layout=(2, 2),
                           size=(1400, 1000),
                           dpi=300)
    
    savefig(uncertainty_plot, "uncertainty_quantification_analysis.pdf")
    savefig(uncertainty_plot, "uncertainty_quantification_analysis.png")
    
    return uncertainty_plot
end

# 创建性能分布图
function create_performance_distribution_plot(mc_results)
    metrics = ["capacity", "energy", "final_voltage", "platform_time"]
    
    plots_array = []
    
    for metric in metrics
        p = plot(title="$(uppercase(metric)) Distribution",
                xlabel="$metric",
                ylabel="Probability Density",
                legend=:topright)
        
        for model in ["P2D", "sP2D"]
            values = [sample[metric] for sample in mc_results[model]]
            
            # 创建直方图
            histogram!(p, values,
                      label="$model Model",
                      alpha=0.6,
                      normalize=:pdf,
                      bins=20)
            
            # 添加均值线
            mean_val = mean(values)
            vline!(p, [mean_val],
                   label="$model Mean",
                   linewidth=2,
                   linestyle=:dash)
        end
        
        push!(plots_array, p)
    end
    
    combined = plot(plots_array...,
                   layout=(2, 2),
                   size=(800, 600))
    
    return combined
end

# 创建鲁棒性对比图
function create_robustness_comparison_plot(mc_results)
    metrics = ["capacity", "energy", "final_voltage", "platform_time"]
    
    p = plot(title="Model Robustness Comparison",
             xlabel="Performance Metrics",
             ylabel="Coefficient of Variation (%)",
             legend=:topright,
             xrotation=45)
    
    p2d_cv = []
    sp2d_cv = []
    
    for metric in metrics
        # 计算变异系数 (CV = std/mean * 100%)
        p2d_values = [sample[metric] for sample in mc_results["P2D"]]
        sp2d_values = [sample[metric] for sample in mc_results["sP2D"]]
        
        p2d_cv_val = std(p2d_values) / mean(p2d_values) * 100
        sp2d_cv_val = std(sp2d_values) / mean(sp2d_values) * 100
        
        push!(p2d_cv, p2d_cv_val)
        push!(sp2d_cv, sp2d_cv_val)
    end
    
    x_positions = 1:length(metrics)
    
    plot!(p, x_positions .- 0.2, p2d_cv,
          label="P2D Model",
          marker=:circle,
          bar=true,
          alpha=0.7,
          color=:blue)
    
    plot!(p, x_positions .+ 0.2, sp2d_cv,
          label="sP2D Model",
          marker=:diamond,
          bar=true,
          alpha=0.7,
          color=:red)
    
    plot!(p, xticks=(x_positions, metrics))
    
    return p
end

# 创建敏感性雷达图
function create_sensitivity_radar_plot(mc_results)
    # 计算每个参数对性能指标的相关性
    param_names = ["扩散系数负极", "扩散系数正极", "反应速率常数负极", 
                   "反应速率常数正极", "电极孔隙率", "固相粒子半径"]
    
    # 对P2D模型进行分析
    p2d_data = mc_results["P2D"]
    correlations = []
    
    for param in param_names
        param_values = [sample["param_$(param)"] for sample in p2d_data]
        capacity_values = [sample["capacity"] for sample in p2d_data]
        
        # 计算相关系数
        correlation = cor(param_values, capacity_values)
        push!(correlations, abs(correlation))
    end
    
    # 创建雷达图
    angles = range(0, 2π, length=length(param_names)+1)[1:end-1]
    
    p = plot(title="Parameter Sensitivity Radar",
             projection=:polar,
             legend=:topright)
    
    # P2D模型
    plot!(p, angles, correlations,
          label="P2D Sensitivity",
          linewidth=2,
          marker=:circle,
          fill=true,
          alpha=0.3)
    
    # sP2D模型
    sp2d_data = mc_results["sP2D"]
    sp2d_correlations = []
    
    for param in param_names
        param_values = [sample["param_$(param)"] for sample in sp2d_data]
        capacity_values = [sample["capacity"] for sample in sp2d_data]
        correlation = cor(param_values, capacity_values)
        push!(sp2d_correlations, abs(correlation))
    end
    
    plot!(p, angles, sp2d_correlations,
          label="sP2D Sensitivity",
          linewidth=2,
          marker=:diamond,
          fill=true,
          alpha=0.3)
    
    return p
end

# 创建不确定性传播图
function create_uncertainty_propagation_plot(mc_results)
    # 分析输入不确定性如何传播到输出
    p = plot(title="Uncertainty Propagation Analysis",
             xlabel="Input Parameter Uncertainty (CV%)",
             ylabel="Output Performance Uncertainty (CV%)",
             legend=:topright,
             grid=true)
    
    # 计算输入参数的变异系数
    param_names = ["扩散系数负极", "扩散系数正极", "反应速率常数负极"]
    
    for model in ["P2D", "sP2D"]
        input_cv = []
        output_cv = []
        
        for param in param_names
            # 输入不确定性
            param_values = [sample["param_$(param)"] for sample in mc_results[model]]
            input_cv_val = std(param_values) / mean(param_values) * 100
            
            # 输出不确定性 (容量)
            capacity_values = [sample["capacity"] for sample in mc_results[model]]
            output_cv_val = std(capacity_values) / mean(capacity_values) * 100
            
            push!(input_cv, input_cv_val)
            push!(output_cv, output_cv_val)
        end
        
        plot!(p, input_cv, output_cv,
              label="$model Model",
              marker=:circle,
              linewidth=2,
              markersize=8)
        
        # 添加参数标签
        for (i, param) in enumerate(param_names)
            annotate!(p, input_cv[i], output_cv[i], 
                     text(param[1:4], 8, :bottom))
        end
    end
    
    # 添加理想线性关系参考
    max_cv = maximum([maximum(input_cv) for input_cv in [[], []]])  # 需要实际计算
    plot!(p, [0, max_cv], [0, max_cv],
          linestyle=:dash,
          color=:black,
          label="1:1 Propagation")
    
    return p
end

# 统计分析和置信区间
function statistical_analysis(mc_results)
    println("=== 不确定性量化统计分析 ===")
    
    for model in ["P2D", "sP2D"]
        println("$model 模型统计:")
        
        # 容量统计
        capacities = [sample["capacity"] for sample in mc_results[model]]
        cap_mean = mean(capacities)
        cap_std = std(capacities)
        cap_cv = cap_std / cap_mean * 100
        
        # 95%置信区间
        cap_ci_lower = quantile(capacities, 0.025)
        cap_ci_upper = quantile(capacities, 0.975)
        
        println("  容量: $(round(cap_mean, digits=3)) ± $(round(cap_std, digits=3)) Ah")
        println("  变异系数: $(round(cap_cv, digits=2))%")
        println("  95%置信区间: [$(round(cap_ci_lower, digits=3)), $(round(cap_ci_upper, digits=3))] Ah")
        
        # 计算时间统计
        calc_times = [sample["calc_time"] for sample in mc_results[model]]
        time_mean = mean(calc_times)
        time_std = std(calc_times)
        
        println("  计算时间: $(round(time_mean, digits=2)) ± $(round(time_std, digits=2)) s")
        println()
    end
    
    # 模型对比
    p2d_caps = [sample["capacity"] for sample in mc_results["P2D"]]
    sp2d_caps = [sample["capacity"] for sample in mc_results["sP2D"]]
    
    # 配对t检验
    cap_differences = sp2d_caps .- p2d_caps
    mean_diff = mean(cap_differences)
    std_diff = std(cap_differences)
    
    println("模型对比:")
    println("  平均容量差异(sP2D-P2D): $(round(mean_diff, digits=4)) Ah")
    println("  差异标准差: $(round(std_diff, digits=4)) Ah")
    println("  相对差异: $(round(abs(mean_diff)/mean(p2d_caps)*100, digits=2))%")
end

# 执行不确定性量化分析
println("开始不确定性量化与鲁棒性分析...")
mc_results = uncertainty_quantification_analysis()

# 生成可视化
uncertainty_plot = create_uncertainty_plots(mc_results)

# 统计分析
statistical_analysis(mc_results)

# 保存蒙特卡洛结果
all_mc_data = []
for model in ["P2D", "sP2D"]
    for sample in mc_results[model]
        push!(all_mc_data, sample)
    end
end

mc_df = DataFrame(all_mc_data)
CSV.write("uncertainty_quantification_results.csv", mc_df)

println("不确定性量化分析完成!") 