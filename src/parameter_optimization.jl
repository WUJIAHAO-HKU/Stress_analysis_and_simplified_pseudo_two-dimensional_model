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
# 适配新的三阶模型的参数化函数
function phie_fit_with_params(x::Array{Float64}, Li::Array{Float64}, ki::Array{Float64}, phie0::Float64,
    cubic_coef::Float64, sep_coef::Float64, harm_coef::Float64, 
    typef::String="sin")
Ln, Ls, Lp = Li
kn, ks, kp = ki
L = Ln + Ls + Lp
phie = zeros(length(x))

if typef == "sin"
for i in eachindex(x)
if x[i] <= Ln
# 负极区域：改进的二次函数
t = x[i]/Ln
# 添加三次项，系数可调整
phie[i] = kn * (0.5 * t^2 - cubic_coef * t^3) - kn * Ln/2 - ks * Ls/2 + phie0
elseif x[i] < Ln + Ls
# 隔膜区域：改进的线性插值
t = (x[i] - Ln)/Ls
# 添加平滑过渡，系数可调整
phie[i] = ks * (x[i] - Ln - Ls/2) + 
sep_coef * ks * Ls * (t^2 - t) + phie0
else
# 正极区域：改进的二次函数
t = (x[i] - (Ln + Ls))/Lp
# 添加三次项，系数可调整
phie[i] = -kp * (0.5 * t^2 - cubic_coef * t^3) + kp * Lp/2 + ks * Ls/2 + phie0
end
end
elseif typef == "quad"
pi = acos(-1.0)
for i in eachindex(x)
if x[i] <= Ln
# 负极区域：改进的余弦函数
t = x[i]/Ln
# 添加高阶谐波，系数可调整
phie[i] = -2 * kn * Ln/pi * (cos(pi*t/2) + harm_coef*cos(3*pi*t/2)) - 
ks * Ls/2 + phie0
elseif x[i] < Ln + Ls
# 隔膜区域：改进的线性插值
t = (x[i] - Ln)/Ls
# 添加平滑过渡，系数可调整
phie[i] = ks * (x[i] - Ln - Ls/2) + 
sep_coef * ks * Ls * (t^2 - t) + phie0
else
# 正极区域：改进的正弦函数
t = (x[i] - (Ln + Ls))/Lp
# 添加高阶谐波，系数可调整
phie[i] = 2 * kp * Lp/pi * (sin(pi*t/2) + harm_coef*sin(3*pi*t/2)) + 
ks * Ls/2 + phie0
end
end
end

# 边界处理以确保连续性
for i in eachindex(x)
if abs(x[i] - Ln) < 1e-10 || abs(x[i] - (Ln + Ls)) < 1e-10
# 在界面处使用相邻点的平均值
if i > 1 && i < length(x)
phie[i] = (phie[i-1] + phie[i+1])/2
end
end
end

return phie
end

# 三阶模型的参数优化函数 - 确保参数名称一致
function optimize_parameters(x::Array{Float64}, Li::Array{Float64}, ki::Array{Float64}, phie0::Float64, 
    reference_solution::Array{Float64}, typef::String="sin")
# 为简化后的三阶模型定义搜索范围
cubic_coef_range = 0.30:0.01:0.35  # 三次项系数范围
sep_coef_range = 0.08:0.01:0.12    # 隔膜过渡系数范围
harm_coef_range = 0.08:0.01:0.12   # 谐波系数范围

best_error = Inf
# 初始化一个包含所有必要参数的字典，避免KeyError
best_params = Dict(
"cubic_coef" => 0.30,  # 默认值，与原始函数一致
"sep_coef" => 0.12,    # 默认值，与原始函数一致
"harm_coef" => 0.08    # 默认值，与原始函数一致
)

# 计算总迭代次数
total_iterations = length(cubic_coef_range) * length(sep_coef_range) * length(harm_coef_range)
current_iteration = 0
last_progress_print = 0

# 网格搜索
for cubic_coef in cubic_coef_range
for sep_coef in sep_coef_range
for harm_coef in harm_coef_range
current_iteration += 1

# 计算当前参数组合的结果
result = phie_fit_with_params(x, Li, ki, phie0, 
                    cubic_coef, sep_coef, harm_coef, typef)

# 计算误差指标
error = norm(result - reference_solution)

# 更新最优参数
if error < best_error
best_error = error
best_params = Dict(
 "cubic_coef" => cubic_coef,
 "sep_coef" => sep_coef,
 "harm_coef" => harm_coef
)
println("New best error: $best_error")
println("Parameters: ", best_params)
end

# 打印进度
progress = 100.0 * current_iteration / total_iterations
if floor(Int, progress) > last_progress_print
println("Optimization progress: $(round(progress, digits=2))%")
last_progress_print = floor(Int, progress)
end
end
end
end

# 确保返回值包含所有必要的键
if !haskey(best_params, "cubic_coef")
best_params["cubic_coef"] = 1/3
end
if !haskey(best_params, "sep_coef")
best_params["sep_coef"] = 0.1
end
if !haskey(best_params, "harm_coef")
best_params["harm_coef"] = 0.1
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