using Plots, CSV, DataFrames, Interpolations, SparseArrays, Statistics
include("D:/竞赛和课程文件/课程文件/毕业设计/SRC/JuBat.jl") 
# 设置数据路径
path = "D:/竞赛和课程文件/课程文件/毕业设计/SRC/data/drive_cycles/"

# 预处理数据函数
function preprocess_data(df, time_col, current_col)
    df = dropmissing(df, [time_col, current_col])
    df = sort(df, time_col)
    df = combine(groupby(df, time_col), current_col => mean => current_col)
    return df
end

# 加载并预处理 UDDS 数据
udds_data = CSV.read(path * "UDDS.csv", DataFrame)
println("UDDS 数据列名: ", names(udds_data))
udds_data = preprocess_data(udds_data, "Time", " Current")
time_udds = udds_data[:, "Time"]
current_udds = udds_data[:, " Current"]

# 加载并预处理 US06 数据
us06_data = CSV.read(path * "US06.csv", DataFrame)
println("US06 数据列名: ", names(us06_data))
us06_data = preprocess_data(us06_data, "time", "current")
time_us06 = us06_data[:, "time"]
current_us06 = us06_data[:, "current"]

# 分段模拟函数
function run_segmented_simulation(param_dim, time_array, current_interp, segment_size=200)
    # 计算需要的段数
    num_segments = ceil(Int, length(time_array) / segment_size)
    
    # 存储各段结果
    all_times = Float64[]
    all_voltages = Float64[]
    all_concentrations = Array{Float64}(undef, 0)
    all_stresses = Array{Float64}(undef, 0)
    
    # 设置电压限制
    v_min = 2.5  # 最小电压限制，通常约2.5-3.0V
    v_max = 4.3  # 最大电压限制，通常约4.2-4.3V
    
    println("电压限制范围: $(v_min)V - $(v_max)V")
    
    # 定义安全电流函数
    function safe_current(t, prev_voltage=nothing)
        # 获取原始电流
        current = current_interp(t)
        
        # 如果没有上一次电压数据，直接返回电流
        if prev_voltage === nothing
            return current
        end
        
        # 如果电压接近边界，降低电流幅度
        if prev_voltage < v_min + 0.1 && current < 0  # 放电电流为负
            # 如果电压已经低，减少放电电流
            reduced_current = current * 1.00
            println("电压低 ($(prev_voltage)V)，降低放电电流：$(current) -> $(reduced_current)")
            return reduced_current
        elseif prev_voltage > v_max - 0.1 && current > 0  # 充电电流为正
            # 如果电压已经高，减少充电电流
            reduced_current = current * 1.00
            println("电压高 ($(prev_voltage)V)，降低充电电流：$(current) -> $(reduced_current)")
            return reduced_current
        else
            return current
        end
    end
    
    # 初始状态
    last_voltage = nothing
    last_time = 0.0
    
    # 按段进行模拟
    for seg in 1:num_segments
        println("开始模拟段 $seg / $num_segments")
        
        # 计算当前段的时间范围
        start_idx = (seg-1) * segment_size + 1
        end_idx = min(seg * segment_size, length(time_array))
        
        # 如果是第一段，从0时刻开始
        if seg == 1
            seg_times = time_array[start_idx:end_idx]
        else
            # 否则从上一段的最后时刻接续
            seg_times = time_array[start_idx:end_idx]
        end
        
        # 设置仿真选项
        opt = JuBat.Option()
        opt.model = "sP2D"
        opt.mechanicalmodel = "full"
        opt.time = seg_times
        
        # 使用带安全检查的电流函数
        local_voltage = last_voltage
        opt.Current = t -> begin
            c = safe_current(t, local_voltage)
            return c
        end
        
        # 创建并运行当前段的仿真
        try
            case = JuBat.SetCase(param_dim, opt)
            result = JuBat.Solve(case)
            
            # 提取结果
            seg_times = result["time [s]"]
            seg_voltages = result["cell voltage [V]"]
            seg_concentrations = result["negative particle surface lithium concentration [mol/m^3]"]
            
            # 提取应力数据（如果有）
            if haskey(result, "negative particle surface tangential stress[Pa]")
                seg_stresses = result["negative particle surface tangential stress[Pa]"]
            else
                seg_stresses = zeros(size(seg_concentrations))
                println("警告: 段 $seg 未找到应力数据键")
            end
            
            # 处理维度
            if ndims(seg_concentrations) > 1
                seg_concentrations_plot = seg_concentrations[1, :]
                seg_stresses_plot = seg_stresses[1, :]
            else
                seg_concentrations_plot = seg_concentrations
                seg_stresses_plot = seg_stresses
            end
            
            # 点级别的错误处理：检查每个时间点的电压，过滤掉无效值
            valid_indices = []
            for i in 1:length(seg_voltages)
                if v_min <= seg_voltages[i] <= v_max
                    push!(valid_indices, i)
                else
                    println("发现无效电压 $(seg_voltages[i])V 在时间 $(seg_times[i])s，跳过该点")
                end
            end
            
            # 如果没有有效点，跳过整个段
            if isempty(valid_indices)
                println("段 $seg 没有有效的电压点，跳过整个段")
                continue
            end
            
            # 只保留有效数据点
            valid_times = seg_times[valid_indices]
            valid_voltages = seg_voltages[valid_indices]
            valid_concentrations = seg_concentrations_plot[valid_indices]
            valid_stresses = seg_stresses_plot[valid_indices]
            
            # 更新最后的有效状态
            last_valid_voltage = valid_voltages[end]
            last_valid_time = valid_times[end]

            
            
            # 合并结果
            append!(all_times, valid_times)
            append!(all_voltages, valid_voltages)
            
            if isempty(all_concentrations)
                all_concentrations = valid_concentrations
                all_stresses = valid_stresses
            else
                all_concentrations = vcat(all_concentrations, valid_concentrations)
                all_stresses = vcat(all_stresses, valid_stresses)
            end

            function smooth_data(data, window_size=5)
                n = length(data)
                smoothed = copy(data)
                for i in 1:n
                    start_idx = max(1, i - window_size ÷ 2)
                    end_idx = min(n, i + window_size ÷ 2)
                    smoothed[i] = mean(data[start_idx:end_idx])
                end
                return smoothed
            end
            all_voltages = smooth_data(all_voltages)
            all_concentrations = smooth_data(all_concentrations)
            all_stresses = smooth_data(all_stresses)
            
            println("段 $seg 完成: 有效点数 = $(length(valid_indices))/$(length(seg_voltages))")
            println("最终有效电压 = $(last_valid_voltage)V, 时间 = $(last_valid_time)s")
            
        catch e
            println("段 $seg 发生错误: $e")
            println("跳过该段，继续下一段...")
            continue
        end
    end
    
    return all_times, all_voltages, all_concentrations, all_stresses
end

# 设置电池参数
param_dim = JuBat.ChooseCell("LG M50")

println("负极杨氏模量: ", param_dim.NE.E)
println("负极泊松比: ", param_dim.NE.nu)
println("负极体积膨胀系数: ", param_dim.NE.Omega)

# 创建电流插值函数
current_interp_udds = LinearInterpolation(time_udds, current_udds, extrapolation_bc=Flat())
current_interp_us06 = LinearInterpolation(time_us06, current_us06, extrapolation_bc=Flat())

# 设置分段大小
segment_size = 100  # 每段100秒

# UDDS分段模拟
println("开始UDDS工况分段模拟")
udds_time_array = collect(0:1:maximum(time_udds))
time_udds_sim, voltage_udds, concentration_udds_plot, stress_udds_plot = 
    run_segmented_simulation(param_dim, udds_time_array, current_interp_udds, segment_size)

# 检查UDDS结果是否为空
if length(time_udds_sim) == 0
    println("警告：UDDS模拟未产生任何结果！")
else
    # 绘制 UDDS 结果图
    p1_udds = plot(time_udds_sim, voltage_udds, label="voltage (V)", xlabel="time (s)", ylabel="voltage (V)", lw=2, title="UDDS 工况")
    p2_udds = plot(time_udds_sim, concentration_udds_plot, label="concentration (mol/m³)", xlabel="time (s)", ylabel="concentration (mol/m³)", lw=2)
    p3_udds = plot(time_udds_sim, stress_udds_plot, label="stress (Pa)", xlabel="time (s)", ylabel="stress (Pa)", lw=2)
    plot_udds = plot(p1_udds, p2_udds, p3_udds, layout=(3, 1), size=(800, 600))
    savefig(plot_udds, "udds_results_segmented_sp2d.pdf")
    
    # 保存UDDS数据到CSV
    udds_df = DataFrame(
        "Time (s)" => time_udds_sim,
        "Voltage (V)" => voltage_udds,
        "Concentration (mol/m³)" => concentration_udds_plot,
        "Stress (Pa)" => stress_udds_plot
    )
    CSV.write("udds_results_segmented_sp2d.csv", udds_df)
    println("UDDS结果已保存")
end

# US06分段模拟
println("开始US06工况分段模拟")
us06_time_array = collect(0:1:maximum(time_us06))
time_us06_sim, voltage_us06, concentration_us06_plot, stress_us06_plot = 
    run_segmented_simulation(param_dim, us06_time_array, current_interp_us06, segment_size)

# 检查US06结果是否为空
if length(time_us06_sim) == 0
    println("警告：US06模拟未产生任何结果！")
else
    # 绘制 US06 结果图
    p1_us06 = plot(time_us06_sim, voltage_us06, label="voltage (V)", xlabel="time (s)", ylabel="voltage (V)", lw=2, title="US06 工况")
    p2_us06 = plot(time_us06_sim, concentration_us06_plot, label="concentration (mol/m³)", xlabel="time (s)", ylabel="concentration (mol/m³)", lw=2)
    p3_us06 = plot(time_us06_sim, stress_us06_plot, label="stress (Pa)", xlabel="time (s)", ylabel="stress (Pa)", lw=2)
    plot_us06 = plot(p1_us06, p2_us06, p3_us06, layout=(3, 1), size=(800, 600))
    savefig(plot_us06, "us06_results_segmented_sp2d.pdf")
    
    # 保存US06数据到CSV
    us06_df = DataFrame(
        "Time (s)" => time_us06_sim,
        "Voltage (V)" => voltage_us06,
        "Concentration (mol/m³)" => concentration_us06_plot,
        "Stress (Pa)" => stress_us06_plot
    )
    CSV.write("us06_results_segmented_sp2d.csv", us06_df)
    println("US06结果已保存")
end

println("所有模拟完成")