using LinearAlgebra

function diagnose_potential_fit(case1, result1, x, optimized_solution)
    # 1. 检查量纲转换
    println("Scaling factor: ", case1.param.scale.phi)
    
    # 2. 检查电位范围
    ref_potential = result1["electrolyte potential [V]"][:,1]
    println("\nPotential Ranges (with units):")
    println("Reference: ", extrema(ref_potential))
    println("Optimized: ", extrema(optimized_solution * case1.param.scale.phi))
    
    # 3. 检查无量纲值
    ref_dimensionless = ref_potential / case1.param.scale.phi
    println("\nDimensionless Ranges:")
    println("Reference: ", extrema(ref_dimensionless))
    println("Optimized: ", extrema(optimized_solution))
    
    # 4. 检查初始值
    println("\nInitial Values:")
    println("Reference: ", ref_potential[1])
    println("Optimized: ", optimized_solution[1] * case1.param.scale.phi)
end

# 使用特定参数的拟合函数
function phie_fit_with_params(x::Array{Float64}, Li::Array{Float64}, ki::Array{Float64}, phie0::Float64,
    δ::Float64, corr::Float64, poly3::Float64, poly4::Float64, trans::Float64,
    typef::String="sin")
Ln, Ls, Lp = Li
kn, ks, kp = ki
L = Ln + Ls + Lp
phie = zeros(length(x))
phie .+= phie0

δn = δ * Ln
δp = δ * Lp

if typef == "sin"
for i in eachindex(x)
if x[i] <= Ln
t = x[i]/Ln
phie[i] =  phie0 + kn * (0.5*t^2 - poly3*t^3 + poly4*t^4) - kn * Ln/2 - ks * Ls/2

if x[i] >= Ln - δn
β = (x[i] - (Ln - δn))/δn
correction = corr * kn * δn * (β^2 * (3 - 2β) + trans*β*(1 - β)*(1 - 2β))
phie[i] += correction
end

elseif x[i] < Ln + Ls
t = (x[i] - Ln)/Ls
base = ks * (x[i] - Ln - Ls/2)
transition = t^2 * (3 - 2t) + trans*t*(1 - t)*(1 - 2t)
correction = corr * ks * Ls * (transition - 0.5)
phie[i] = base + correction + phie0

else
t = (x[i] - (Ln + Ls))/Lp
phie[i] = -kp * (0.5*t^2 - poly3*t^3 + poly4*t^4) + kp * Lp/2 + ks * Ls/2 + phie0

if x[i] <= Ln + Ls + δp
β = (x[i] - (Ln + Ls))/δp
correction = corr * kp * δp * (β^2 * (3 - 2β) + trans*β*(1 - β)*(1 - 2β))
phie[i] += correction
end
end
end
elseif typef == "quad"
pi = acos(-1.0)
for i in eachindex(x)
if x[i] <= Ln
t = x[i]/Ln
phie[i] = -2 * kn * Ln/pi * (
cos(pi*t/2) + 
poly3*cos(3*pi*t/2) + 
poly4*cos(5*pi*t/2)
) - ks * Ls/2 + phie0

if x[i] >= Ln - δn
β = (x[i] - (Ln - δn))/δn
correction = corr * kn * δn * (sin(pi*β) + trans*sin(2*pi*β)*(1 - β))
phie[i] += correction
end

elseif x[i] < Ln + Ls
t = (x[i] - Ln)/Ls
base = ks * (x[i] - Ln - Ls/2)
transition = 0.5 * (1 - cos(pi*t)) + trans * sin(2*pi*t)*(1 - t)*t
correction = corr * ks * Ls * (transition - 0.5)
phie[i] = base + correction + phie0

else
t = (x[i] - (Ln + Ls))/Lp
phie[i] = 2 * kp * Lp/pi * (
sin(pi*t/2) + 
poly3*sin(3*pi*t/2) + 
poly4*sin(5*pi*t/2)
) + ks * Ls/2 + phie0

if x[i] <= Ln + Ls + δp
β = (x[i] - (Ln + Ls))/δp
correction = corr * kp * δp * (sin(pi*β) + trans*sin(2*pi*β)*(1 - β))
phie[i] += correction
end
end
end
end

# 边界处理
for i in eachindex(x)
if abs(x[i] - Ln) < 1e-10 || abs(x[i] - (Ln + Ls)) < 1e-10
if i > 2 && i < length(x) - 1
w1, w2, w3, w4 = 0.1165, 0.3835, 0.3835, 0.1165
phie[i] = w1*phie[i-2] + w2*phie[i-1] + w3*phie[i+1] + w4*phie[i+2]
end
end
end

return phie
end



# 参数优化函数
function optimize_parameters(x::Array{Float64}, Li::Array{Float64}, ki::Array{Float64}, phie0::Float64, 
                           reference_solution::Array{Float64}, typef::String="sin")
    # 定义搜索范围（以当前最优值为中心的小范围）
    δ_range = 0.0575:0.0001:0.0577
    corr_range = 0.1580:0.0001:0.1590
    poly3_range = 0.3450:0.0001:0.3460
    poly4_range = 0.1050:0.0001:0.1060
    trans_range = 0.1200:0.0001:0.1220
    
    best_error = Inf
    best_params = Dict()
    
    # 记录优化进度
    total_iterations = length(δ_range) * length(corr_range) * length(poly3_range) * 
                      length(poly4_range) * length(trans_range)
    current_iteration = 0
    
    # 网格搜索
    for δ in δ_range
        for corr in corr_range
            for poly3 in poly3_range
                for poly4 in poly4_range
                    for trans in trans_range
                        current_iteration += 1
                        if current_iteration % 1000 == 0
                            progress = 100 * current_iteration / total_iterations
                            println("Optimization progress: $(round(progress, digits=2))%")
                        end
                        
                        # 使用当前参数计算结果
                        result = phie_fit_with_params(x, Li, ki, phie0, δ, corr, poly3, poly4, trans, typef)
                        
                        # 计算误差
                        error = norm(result - reference_solution)
                        # 更新最优参数
                        if error < best_error
                            best_error = error
                            best_params = Dict(
                                "δ" => δ,
                                "corr" => corr,
                                "poly3" => poly3,
                                "poly4" => poly4,
                                "trans" => trans
                            )
                            println("New best error: $best_error")
                            println("Parameters: ", best_params)
                        end
                    end
                end
            end
        end
    end
    
    return best_params, best_error
end

# 使用示例
function test_optimization()
    # 这里需要提供实际的测试数据
    x = collect(0.0:0.01:1.0)  # 示例空间网格
    Li = [0.3, 0.4, 0.3]      # 示例长度参数
    ki = [1.0, 1.0, 1.0]      # 示例电导率参数
    phie0 = 0.0               # 示例初始电位
    
    # 生成参考解（这里需要替换为实际的参考解）
    reference_solution = zeros(length(x))  # 替换为实际的参考解
    
    # 运行优化
    best_params, best_error = optimize_parameters(x, Li, ki, phie0, reference_solution)
    
    println("Optimization completed!")
    println("Best parameters: ", best_params)
    println("Best error: ", best_error)
    
    return best_params, best_error
end