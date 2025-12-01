using Plots, CSV, DataFrames, Interpolations, Statistics, StatsBase
include("../src/JuBat.jl") 

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")
param_dim.cell.v_h = 4.3

# 参数设置
total_time = 3000      # 总时长(秒)
time_step = 0.1       # 100ms步长
time = collect(0:time_step:total_time)

# 电池参数
I_base = 1            # 1C基准电流(1A)
charge_peaks = [1.2, 0.8, 0.4]  # 三阶段充电峰值电流(A)
discharge_peaks = [2.4, 1.6, 0.8]  # 三阶段放电峰值电流(A)
cycle_duration = 300  # 充放电周期300秒

# 阶段时间划分 (单位:秒)
discharge_stages = [60, 120, 180]  # 放电三阶段截止时间：60s/120s/180s
charge_stages = [240, 285, 300]    # 充电三阶段截止时间：240s/285s/300s

current = zeros(length(time))
#temp = 25.0 .+ zeros(length(time))  # 温度模拟(℃)

for i in 1:length(time)
    t = time[i]
    cycle_phase = mod(t, cycle_duration)
    
    # ===== 放电阶段 (0-180s) =====
    if cycle_phase < discharge_stages[end]
        # 阶段1：强放电 (0-60s)
        if cycle_phase < discharge_stages[1]
            I = discharge_peaks[1] * (1 - 0.1*rand())
            
        # 阶段2：中放电 (60-120s)
        elseif cycle_phase < discharge_stages[2]
            I = discharge_peaks[2] * (1 - 0.05*rand())
            
        # 阶段3：弱放电 (120-180s)
        else
            I = discharge_peaks[3] * (1 - 0.02*rand())
        end
        
    # ===== 充电阶段 (180-300s) =====
    else
        # 阶段1：快充 (180-240s)
        if cycle_phase < charge_stages[1]
            I = -charge_peaks[1] * (1 - 0.1*randn())
            
        # 阶段2：脉冲充 (240-285s)
        elseif cycle_phase < charge_stages[2]
            pulse = charge_peaks[2] * (0.8 + 0.2*sign(sin(2π*10*t))) # 10Hz简化脉冲
            I = -pulse * (1 - 0.05*randn())
            
        # 阶段3：涓流充 (285-300s)
        else
            I = -charge_peaks[3] * (1 - 0.02*randn())
        end
    end
    
    current[i] = I
end

# 添加高频噪声（模拟测量噪声）
current .+= 0.1*randn(length(current))

# 创建电流插值函数
itp = linear_interpolation(time, current, extrapolation_bc=Flat())


# 设置模拟选项
opt = JuBat.Option()
opt.mechanicalmodel = "full"
opt.model = "P2D"  # 使用P2D模型
opt.time = time  # 使用与电流数据相同的时间步长

# 设置电流函数
opt.Current = t -> itp(t)  # 正确用法

# 创建和运行P2D模型模拟
println("开始模拟快速充放电电流工况... P2D模型")
case = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case)
println("P2D模型模拟完成")

# 测量P2D模型计算时间
p2d_time = @elapsed begin
    case = JuBat.SetCase(param_dim, opt)
    result = JuBat.Solve(case)
end
println("P2D模型计算时间: ", p2d_time, " 秒")

# 提取P2D模型结果数据
time_result = result["time [s]"]
x_grid_n = case.mesh["negative electrode"].node[:, 1] # 使用 node 字段
x_grid_p = case.mesh["positive electrode"].node[:, 1]  # 正极网格，单位：米 (m)

#设置计算你的时间点
t_index1 = argmin(abs.(time_result .- 60))
t_index2 = argmin(abs.(time_result .- 120))
t_index3 = argmin(abs.(time_result .- 180))
t_index4 = argmin(abs.(time_result .- 240))
t_index5 = argmin(abs.(time_result .- 300))
t_index6 = argmin(abs.(time_result .- 360))
t_index7 = argmin(abs.(time_result .- 420))
t_index8 = argmin(abs.(time_result .- 480))
t_index9 = argmin(abs.(time_result .- 540))
t_index10 = argmin(abs.(time_result .- 600))

voltage_p2d = result["cell voltage [V]"]
concentration_p2d = result["negative particle surface lithium concentration [mol/m^3]"]
stress_p2d = result["negative particle surface tangential stress[Pa]"]
concentration_p_p2d = result["positive particle surface lithium concentration [mol/m^3]"]
stress_p_p2d = result["positive particle surface tangential stress[Pa]"]

# 设置sP2D模型
opt.model = "sP2D"  # 使用sP2D模型
println("开始模拟快速充放电电流工况... sP2D模型")
case_sP2D = JuBat.SetCase(param_dim, opt)
result_sP2D = JuBat.Solve(case_sP2D)
println("sP2D模型模拟完成")

# 测量sP2D模型计算时间
sp2d_time = @elapsed begin
    case_sP2D = JuBat.SetCase(param_dim, opt)
    result_sP2D = JuBat.Solve(case_sP2D)
end

println("P2D模型计算时间: ", p2d_time, " 秒")
println("sP2D模型计算时间: ", sp2d_time, " 秒")

# 提取sP2D模型结果数据
voltage_sP2D = result_sP2D["cell voltage [V]"]
concentration_sP2D = result_sP2D["negative particle surface lithium concentration [mol/m^3]"]
stress_sP2D = result_sP2D["negative particle surface tangential stress[Pa]"]
concentration_p_sP2D = result_sP2D["positive particle surface lithium concentration [mol/m^3]"]
stress_p_sP2D = result_sP2D["positive particle surface tangential stress[Pa]"]

#60s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d1 = concentration_p2d[:, t_index1]
end

if ndims(stress_p2d) > 1
    stress_p2d1 = stress_p2d[:, t_index1]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d1 = concentration_p_p2d[:, t_index1]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d1 = stress_p_p2d[:, t_index1]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D1 = concentration_sP2D[:, t_index1]
end

if ndims(stress_sP2D) > 1
    stress_sP2D1 = stress_sP2D[:, t_index1]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D1 = concentration_p_sP2D[:, t_index1]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D1 = stress_p_sP2D[:, t_index1]
end

#120s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d2 = concentration_p2d[:, t_index2]
end

if ndims(stress_p2d) > 1
    stress_p2d2 = stress_p2d[:, t_index2]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d2 = concentration_p_p2d[:, t_index2]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d2 = stress_p_p2d[:, t_index2]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D2 = concentration_sP2D[:, t_index2]
end

if ndims(stress_sP2D) > 1
    stress_sP2D2 = stress_sP2D[:, t_index2]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D2 = concentration_p_sP2D[:, t_index2]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D2 = stress_p_sP2D[:, t_index2]
end

#180s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d3 = concentration_p2d[:, t_index3]
end

if ndims(stress_p2d) > 1
    stress_p2d3 = stress_p2d[:, t_index3]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d3 = concentration_p_p2d[:, t_index3]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d3 = stress_p_p2d[:, t_index3]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D3 = concentration_sP2D[:, t_index3]
end

if ndims(stress_sP2D) > 1
    stress_sP2D3 = stress_sP2D[:, t_index3]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D3 = concentration_p_sP2D[:, t_index3]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D3 = stress_p_sP2D[:, t_index3]
end

#240s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d4 = concentration_p2d[:, t_index4]
end

if ndims(stress_p2d) > 1
    stress_p2d4 = stress_p2d[:, t_index4]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d4 = concentration_p_p2d[:, t_index4]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d4 = stress_p_p2d[:, t_index4]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D4 = concentration_sP2D[:, t_index4]
end

if ndims(stress_sP2D) > 1
    stress_sP2D4 = stress_sP2D[:, t_index4]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D4 = concentration_p_sP2D[:, t_index4]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D4 = stress_p_sP2D[:, t_index4]
end

#300s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d5 = concentration_p2d[:, t_index5]
end

if ndims(stress_p2d) > 1
    stress_p2d5 = stress_p2d[:, t_index5]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d5 = concentration_p_p2d[:, t_index5]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d5 = stress_p_p2d[:, t_index5]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D5 = concentration_sP2D[:, t_index5]
end

if ndims(stress_sP2D) > 1
    stress_sP2D5 = stress_sP2D[:, t_index5]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D5 = concentration_p_sP2D[:, t_index5]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D5 = stress_p_sP2D[:, t_index5]
end

#360s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d6 = concentration_p2d[:, t_index6]
end

if ndims(stress_p2d) > 1
    stress_p2d6 = stress_p2d[:, t_index6]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d6 = concentration_p_p2d[:, t_index6]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d6 = stress_p_p2d[:, t_index6]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D6 = concentration_sP2D[:, t_index6]
end

if ndims(stress_sP2D) > 1
    stress_sP2D6 = stress_sP2D[:, t_index6]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D6 = concentration_p_sP2D[:, t_index6]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D6 = stress_p_sP2D[:, t_index6]
end

#420s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d7 = concentration_p2d[:, t_index7]
end

if ndims(stress_p2d) > 1
    stress_p2d7 = stress_p2d[:, t_index7]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d7 = concentration_p_p2d[:, t_index7]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d7 = stress_p_p2d[:, t_index7]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D7 = concentration_sP2D[:, t_index7]
end

if ndims(stress_sP2D) > 1
    stress_sP2D7 = stress_sP2D[:, t_index7]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D7 = concentration_p_sP2D[:, t_index7]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D7 = stress_p_sP2D[:, t_index7]
end

#480s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d8 = concentration_p2d[:, t_index8]
end

if ndims(stress_p2d) > 1
    stress_p2d8 = stress_p2d[:, t_index8]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d8 = concentration_p_p2d[:, t_index8]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d8 = stress_p_p2d[:, t_index8]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D8 = concentration_sP2D[:, t_index8]
end

if ndims(stress_sP2D) > 1
    stress_sP2D8 = stress_sP2D[:, t_index8]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D8 = concentration_p_sP2D[:, t_index8]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D8 = stress_p_sP2D[:, t_index8]
end

#540s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d9 = concentration_p2d[:, t_index9]
end

if ndims(stress_p2d) > 1
    stress_p2d9 = stress_p2d[:, t_index9]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d9 = concentration_p_p2d[:, t_index9]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d9 = stress_p_p2d[:, t_index9]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D9 = concentration_sP2D[:, t_index9]
end

if ndims(stress_sP2D) > 1
    stress_sP2D9 = stress_sP2D[:, t_index9]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D9 = concentration_p_sP2D[:, t_index9]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D9 = stress_p_sP2D[:, t_index9]
end

#600s
# 处理P2D模型结果的维度
if ndims(concentration_p2d) > 1
    concentration_p2d10 = concentration_p2d[:, t_index10]
end

if ndims(stress_p2d) > 1
    stress_p2d10 = stress_p2d[:, t_index10]
end

if ndims(concentration_p_p2d) > 1
    concentration_p_p2d10 = concentration_p_p2d[:, t_index10]
end

if ndims(stress_p_p2d) > 1
    stress_p_p2d10 = stress_p_p2d[:, t_index10]
end

# 处理sP2D模型结果的维度
if ndims(concentration_sP2D) > 1
    concentration_sP2D10 = concentration_sP2D[:, t_index10]
end

if ndims(stress_sP2D) > 1
    stress_sP2D10 = stress_sP2D[:, t_index10]
end

if ndims(concentration_p_sP2D) > 1
    concentration_p_sP2D10 = concentration_p_sP2D[:, t_index10]
end

if ndims(stress_p_sP2D) > 1
    stress_p_sP2D10 = stress_p_sP2D[:, t_index10]
end

# 创建电压图，包括P2D和sP2D的结果
p1 = plot(time_result, voltage_p2d, 
    label="voltage (P2D)", 
    xlabel="time [s]", 
    ylabel="voltage [V]", 
    lw=2, color=:blue)

plot!(time_result, voltage_sP2D, 
    label="voltage (sP2D)", 
    title="voltage_comparison (P2D vs sP2D)",
    linestyle=:dash, color=:red)

# 创建浓度图
#1
p2 = plot(x_grid_p, concentration_p_p2d1, 
    label="t=60s", 
    legend=:outerright,
    xlabel="x [m]", 
    ylabel="Concentration [mol/m³]", 
    lw=2, color=:orange, alpha=0.1,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D1, 
    label="", 
    title="Electrode Concentration(P2D & SP2D)", 
    linestyle=:dot, 
    shape =:diamond, alpha=0.1,
    color=:red)
plot!(x_grid_n, concentration_p2d1, 
    label="", alpha=0.1,
    lw=2, color=:orange)
plot!(x_grid_n, concentration_sP2D1, 
    label="", alpha=0.1,
    linestyle=:dot, 
    shape =:diamond,
    color=:red)

#2
plot!(x_grid_p, concentration_p_p2d2, 
    label = "t=120s",alpha=0.2,
    lw=2, color=:orange,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D2, label="", alpha=0.2,
    linestyle=:dot,  shape =:diamond,color=:red)
plot!(x_grid_n, concentration_p2d2,  label="", alpha=0.2,
    lw=2, color=:orange)
plot!(x_grid_n, concentration_sP2D2,label="",  alpha=0.2, 
    linestyle=:dot,  shape =:diamond,color=:red)

#3
plot!(x_grid_p, concentration_p_p2d3, label = "t=180s", alpha=0.3,
    lw=2, color=:orange,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D3,  label="", alpha=0.3,
    linestyle=:dot,  shape =:diamond,color=:red)
plot!(x_grid_n, concentration_p2d3, label="", alpha=0.3,
    lw=2, color=:orange)
plot!(x_grid_n, concentration_sP2D3, label="",  alpha=0.3,
    linestyle=:dot,  shape =:diamond,color=:red)

#4
plot!(x_grid_p, concentration_p_p2d4, label = "t=240s",alpha=0.4,
    lw=2, color=:blue,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D4, label="", alpha=0.4,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)
plot!(x_grid_n, concentration_p2d4,  label="", alpha=0.4,
    lw=2, color=:blue)
plot!(x_grid_n, concentration_sP2D4,  label="", alpha=0.4,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)

#5
plot!(x_grid_p, concentration_p_p2d5,  label = "t=300s",alpha=0.5,
    lw=2, color=:blue,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D5,  label="", alpha=0.5,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)
plot!(x_grid_n, concentration_p2d5, label="", alpha=0.5,
    lw=2, color=:blue)
plot!(x_grid_n, concentration_sP2D5,  label="", alpha=0.5,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)

#6
plot!(x_grid_p, concentration_p_p2d6,  label = "t=360s",alpha=0.6,
    lw=2, color=:orange,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D6, label="",  alpha=0.6,
    linestyle=:dot,  shape =:diamond,color=:red)
plot!(x_grid_n, concentration_p2d6,label="",  alpha=0.6,
    lw=2, color=:orange)
plot!(x_grid_n, concentration_sP2D6, label="",  alpha=0.6,
    linestyle=:dot,  shape =:diamond,color=:red)

#7
plot!(x_grid_p, concentration_p_p2d7,  label = "t=420s",alpha=0.7,
    lw=2, color=:orange,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D7, label="",  alpha=0.7,
    linestyle=:dot,  shape =:diamond,color=:red)
plot!(x_grid_n, concentration_p2d7, label="", alpha=0.7,
    lw=2, color=:orange)
plot!(x_grid_n, concentration_sP2D7, label="", alpha=0.7,
    linestyle=:dot,  shape =:diamond,color=:red)

#8
plot!(x_grid_p, concentration_p_p2d8, label = "t=480s",alpha=0.8,
    lw=2, color=:orange,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D8,  label="", alpha=0.8,
    linestyle=:dot,  shape =:diamond,color=:red)
plot!(x_grid_n, concentration_p2d8, label="", alpha=0.8,
    lw=2, color=:orange)
plot!(x_grid_n, concentration_sP2D8,  label="", alpha=0.8,
    linestyle=:dot,  shape =:diamond,color=:red)

#9
plot!(x_grid_p, concentration_p_p2d9,  label = "t=540s",alpha=0.9,
    lw=2, color=:blue,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D9, label="",  alpha=0.9,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)
plot!(x_grid_n, concentration_p2d9,  label="", alpha=0.9,
    lw=2, color=:blue)
plot!(x_grid_n, concentration_sP2D9, label="",  alpha=0.9,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)

#10
plot!(x_grid_p, concentration_p_p2d10,  label = "t=600s",alpha=1,
    lw=2, color=:blue,xlims=(0,1))
plot!(x_grid_p, concentration_p_sP2D10,  label="", alpha=1,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)
plot!(x_grid_n, concentration_p2d10,  label="", alpha=1,
    lw=2, color=:blue)
plot!(x_grid_n, concentration_sP2D10,  label="", alpha=1,
    linestyle=:dot,  shape =:diamond,color=:fuchsia)

######################################################################################
# 创建应力图
#1
p3 = plot(x_grid_p, stress_p_p2d1, 
    label="t=60s", 
    legend=:outerright,
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    title="Electrode Stress(P2D & SP2D)", 
    lw=2, color=:darkorange,alpha=0.1,xlims=(0,1))
plot!(x_grid_p, stress_p_sP2D1, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.1,color=:yellow)
plot!(x_grid_n, stress_sP2D1, 
    label="", 
    alpha=0.1,color=:darkorange)
plot!(x_grid_n, stress_sP2D1, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.1,color=:yellow)

#2
plot!(x_grid_p, stress_p_p2d2, 
    label="t=120s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.2,color=:darkorange)
plot!(x_grid_p, stress_p_sP2D2, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.2,color=:yellow)
plot!(x_grid_n, stress_sP2D2, 
    label="", 
    alpha=0.2,color=:darkorange)
plot!(x_grid_n, stress_sP2D2, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.2,color=:yellow)

#3
plot!(x_grid_p, stress_p_p2d3, 
    label="t=180s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.3,color=:darkorange)
plot!(x_grid_p, stress_p_sP2D3, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.3,color=:yellow)
plot!(x_grid_n, stress_sP2D3, 
    label="", 
    alpha=0.3,color=:darkorange)
plot!(x_grid_n, stress_sP2D3, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.3,color=:yellow)

#4
plot!(x_grid_p, stress_p_p2d4, 
    label="t=240s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.4,color=:"#4db8ff")
plot!(x_grid_p, stress_p_sP2D4, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.4,color=:"#3c9d4d")
plot!(x_grid_n, stress_sP2D4, 
    label="", 
    alpha=0.4,color=:"#4db8ff")
plot!(x_grid_n, stress_sP2D4, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.4,color=:"#3c9d4d")

#5
plot!(x_grid_p, stress_p_p2d5, 
    label="t=300s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.5,color=:"#4db8ff")
plot!(x_grid_p, stress_p_sP2D5, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.5,color=:"#3c9d4d")
plot!(x_grid_n, stress_sP2D5, 
    label="", 
    alpha=0.5,color=:"#4db8ff")
plot!(x_grid_n, stress_sP2D5, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.5,color=:"#3c9d4d")

#6
plot!(x_grid_p, stress_p_p2d6, 
    label="t=360s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.6,color=:darkorange)
plot!(x_grid_p, stress_p_sP2D6, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.6,color=:yellow)
plot!(x_grid_n, stress_sP2D6, 
    label="", 
    alpha=0.6,color=:darkorange)
plot!(x_grid_n, stress_sP2D6, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.6,color=:yellow)

#7
plot!(x_grid_p, stress_p_p2d7, 
    label="t=420s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.7,color=:darkorange)
plot!(x_grid_p, stress_p_sP2D7, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.7,color=:yellow)
plot!(x_grid_n, stress_sP2D7, 
    label="", 
    alpha=0.7,color=:darkorange)
plot!(x_grid_n, stress_sP2D7, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.7,color=:yellow)

#8
plot!(x_grid_p, stress_p_p2d8, 
    label="t=480s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.8,color=:darkorange)
plot!(x_grid_p, stress_p_sP2D8, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.8,color=:yellow)
plot!(x_grid_n, stress_sP2D8, 
    label="", 
    alpha=0.8,color=:darkorange)
plot!(x_grid_n, stress_sP2D8, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.8,color=:yellow)

#9
plot!(x_grid_p, stress_p_p2d9, 
    label="t=540s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=0.9,color=:"#4db8ff")
plot!(x_grid_p, stress_p_sP2D9, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.9,color=:"#3c9d4d")
plot!(x_grid_n, stress_sP2D9, 
    label="", 
    alpha=0.9,color=:"#4db8ff")
plot!(x_grid_n, stress_sP2D9, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=0.9,color=:"#3c9d4d")

#10
plot!(x_grid_p, stress_p_p2d10, 
    label="t=600s", 
    xlabel="x [m]", 
    ylabel="Stress [Pa]", 
    lw=2, alpha=1,color=:"#4db8ff")
plot!(x_grid_p, stress_p_sP2D10, 
    label="", 
    linestyle=:dot, shape=:star,alpha=1,color=:"#3c9d4d")
plot!(x_grid_n, stress_sP2D10, 
    label="", 
    alpha=1,color=:"#4db8ff")
plot!(x_grid_n, stress_sP2D10, 
    label="", 
    linestyle=:dot,  shape=:star,alpha=1,color=:"#3c9d4d")

# 创建电流图
p4 = plot(time, current, 
    label="Current", 
    xlabel="time [s]", 
    ylabel="Current [A]", 
    title="Pulse Waveform Current Profile",
    lw=2, size=(1600, 400), color=:orange)

# 组合图表
plot_combined = plot(p1, p2, p3, p4, 
layout=(4, 1), 
size=(800, 800))

# 保存图表
savefig(p1, "dramatic_voltage_comparison.pdf")
savefig(p2, "dramatic_concentration_comparison.pdf")
savefig(p3, "dramatic_stress_comparison.pdf")
savefig(p4, "dramatic_current_profile.pdf")

savefig(plot_combined, "dramatic_battery_analysis_comparison_super.pdf")

# 计算绝对值差
voltage_error = abs.((voltage_sP2D) .- (voltage_p2d))
concentration_error = abs.((concentration_sP2D) .- (concentration_p2d))
stress_error = abs.((stress_sP2D) .- (stress_p2d))

# 计算误差百分比
voltage_error_percentage = (voltage_error ./ voltage_p2d) .* 100
concentration_error_percentage = (concentration_error ./ concentration_p2d) .* 100
stress_error_percentage = (stress_error ./ stress_p2d) .* 100

# 获取最大误差
max_voltage_error = maximum(voltage_error)
max_concentration_error = maximum(concentration_error)
max_stress_error = maximum(stress_error)

# 获取最大误差百分比
max_voltage_error_percentage = maximum(voltage_error_percentage)
max_concentration_error_percentage = maximum(concentration_error_percentage)
max_stress_error_percentage = maximum(stress_error_percentage)

max_current = maximum(current)

# 计算平均绝对误差
mean_voltage_error = mean(abs.(voltage_error))
mean_concentration_error = mean(abs.(concentration_error))
mean_stress_error = mean(abs.(stress_error))

# 计算平均误差百分比（排除除以零的情况）
function safe_mean_percentage(errors, references)
    valid_indices = findall(x -> x != 0, references)
    isempty(valid_indices) ? 0.0 : mean(abs.(errors[valid_indices] ./ references[valid_indices])) * 100
end

mean_voltage_error_percentage = safe_mean_percentage(voltage_error, voltage_p2d)
mean_concentration_error_percentage = safe_mean_percentage(concentration_error, concentration_p2d)
mean_stress_error_percentage = safe_mean_percentage(stress_error, stress_p2d)

# 确保所有数组长度一致
min_length = min(length(time_result), length(voltage_error), 
                length(concentration_error), length(stress_error))
time_plot = time_result[1:min_length]
voltage_error_plot = voltage_error[1:min_length]
concentration_error_plot = concentration_error[1:min_length]
stress_error_plot = stress_error[1:min_length]

# 创建误差百分比变化图
p5 = plot(time_plot, voltage_error_plot, 
    label="Voltage Error", 
    xlabel="Time [s]", 
    ylabel="Error",
    title="Model Error Comparison (P2D vs sP2D)",
    lw=2, color=:green)  # 调整坐标范围

plot!(time_plot, concentration_error_plot, 
    label="Concentration Error", 
    lw=2, color=:blue)

# 添加应力误差散点
scatter!(time_plot, stress_error_plot, 
    label="Stress Error", 
    markershape=:diamond, markercolor=:red, 
    markersize=2)

savefig(p5, "dramatic_battery_error_comparison.pdf")

# 打印最大误差和最大误差百分比
println("最大电压误差: $max_voltage_error V")
println("平均电压误差: $mean_voltage_error V")
println("最大电压误差百分比: $max_voltage_error_percentage %")
println("平均电压误差百分比: $mean_voltage_error_percentage %")
println("最大浓度误差: $max_concentration_error mol/m³")
println("平均浓度误差: $mean_concentration_error mol/m³")
println("最大浓度误差百分比: $max_concentration_error_percentage %")
println("平均浓度误差百分比: $mean_concentration_error_percentage %")
println("最大应力误差: $max_stress_error Pa")
println("平均应力误差: $mean_stress_error Pa")
println("最大应力误差百分比: $max_stress_error_percentage %")
println("平均应力误差百分比: $mean_stress_error_percentage %")

# 保存数据到CSV
results_df = DataFrame(
    "Time (s)" => time,
    "Voltage (P2D) (V)" => voltage_p2d,
    "Voltage (sP2D) (V)" => voltage_sP2D,
    "Voltage Error (V)" => voltage_error,
    "average Voltage Error (V)" => mean_voltage_error,
    "Voltage Error Percentage (%)" => voltage_error_percentage,
    "average Voltage Error Percentage (%)" => mean_voltage_error_percentage,
    "Concentration (P2D) (mol/m³)" => concentration_p2d,
    "Concentration (sP2D) (mol/m³)" => concentration_sP2D,
    "Concentration Error (mol/m³)" => concentration_error,
    "average Concentration Error (mol/m³)" => mean_concentration_error,
    "Concentration Error Percentage (%)" => concentration_error_percentage,
    "average Concentration Error Percentage (%)" => mean_concentration_error_percentage,
    "Stress (P2D) (Pa)" => stress_p2d,
    "Stress (sP2D) (Pa)" => stress_sP2D,
    "Stress Error (Pa)" => stress_error,
    "average Stress Error (Pa)" => mean_stress_error,
    "Stress Error Percentage (%)" => stress_error_percentage,
    "average Stress Error Percentage (%)" => mean_stress_error_percentage,
    "Current (A)" => [itp(t) for t in time]
)

CSV.write("D:\\竞赛和课程文件\\课程文件\\毕业设计\\dramatic_battery_results_comparison.csv", results_df)

println("分析完成，结果已保存")

# 引用信息
JuBat.Citation()