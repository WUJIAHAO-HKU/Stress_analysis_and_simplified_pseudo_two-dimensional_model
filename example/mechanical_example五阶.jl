using Plots, CSV, DataFrames 
include("../src/JuBat.jl") 
include("../src/parameter_optimization.jl")  # 添加参数优化模块

# 原有的模型设置
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
path = "D:\\竞赛和课程文件\\课程文件\\毕业设计\\src\\data\\"
Crate = 1
i = 5*Crate
opt.Current = x-> i
opt.time = [0 3600]
opt.dt = [0.01, 2]
opt.model = "sP2D"
opt.mechanicalmodel = "full"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)

# 修改后的参数优化函数
function optimize_phie_parameters(case1, result1)
        # 从case1和result1中提取需要的参数
        mesh_el = case1.mesh["electrolyte"]
        x = mesh_el.node  # 空间网格
        Li = [case1.param.NE.thickness, case1.param.SP.thickness, case1.param.PE.thickness]
        
        # 获取电导率参数
        kappa_n_eff = case1.param.EL.kappa(1.0) * case1.param.NE.eps^case1.param.NE.brugg
        kappa_s_eff = case1.param.EL.kappa(1.0) * case1.param.SP.eps^case1.param.SP.brugg
        kappa_p_eff = case1.param.EL.kappa(1.0) * case1.param.PE.eps^case1.param.PE.brugg
        
        ki = [kappa_n_eff, kappa_s_eff, kappa_p_eff]
        
        # 获取电解质电位数据（使用第一个时间步的数据）
        phie_n = result1["electrolyte potential in negative electrode [V]"][:,1]
        phie_p = result1["electrolyte potential in positive electrode [V]"][:,1]
        phie_full = result1["electrolyte potential [V]"][:,1]
        
        # 获取初始电位
        phie0 = phie_n[1]
        
        # 构建参考解
        reference_solution = phie_full / case1.param.scale.phi  # 转换回无量纲形式
    
        println("Starting parameter optimization...")
        best_params, best_error = optimize_parameters(x, Li, ki, phie0, reference_solution)
    
        # 输出优化结果
        println("\nOptimization Results:")
        println("Best parameters:")
        for (key, value) in best_params
            println("$key = $value")
        end
        println("Best error: $best_error")
    
        # 保存优化参数
        param_df = DataFrame(
            parameter = collect(keys(best_params)),
            value = collect(values(best_params))
        )
        CSV.write(path * "optimal_parameters.csv", param_df)
    
        # 可视化比较
        optimized_solution = phie_fit_with_params(x, Li, ki, phie0, 
                                            best_params["δ"], 
                                            best_params["corr"], 
                                            best_params["poly3"], 
                                            best_params["poly4"],
                                            best_params["poly5"],
                                            best_params["trans"],
                                            best_params["harm7"])
        
        # 绘制对比图（转换回有量纲形式进行比较）
        p = plot(x, reference_solution * case1.param.scale.phi, 
                label="Reference Solution", 
                title="Electrolyte Potential Optimization Results",
                xlabel="Position (m)",
                ylabel="Potential (V)")
        plot!(x, optimized_solution * case1.param.scale.phi, 
              label="Optimized Fit")
        savefig(p, path * "optimization_results.png")
    
        return best_params, best_error
    end
    
    # 运行参数优化
    best_params, best_error = optimize_phie_parameters(case1, result1)
    
    # 原有的应力分析代码
    stress_t_vector = vec(result1["negative particle surface tangential stress[Pa]"])
    maxt_df = maximum(stress_t_vector)
    stress_t_df = DataFrame(stress_t = stress_t_vector)
    println("Max tangential stress: ", maxt_df)
    CSV.write(path * "sig_tsurf_test1111.csv", stress_t_df)
    