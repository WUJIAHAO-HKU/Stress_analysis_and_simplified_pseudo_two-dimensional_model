using Plots, CSV, DataFrames, FFTW, Statistics
include("../src/JuBat.jl")

# 频域分析与传递函数验证
function frequency_domain_analysis()
    println("开始频域分析...")
    
    # 定义频率范围 (Hz)
    frequencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    models = ["P2D", "sP2D"]
    
    # 存储频域响应数据
    frequency_response = Dict()
    
    # 设置电池参数
    param_dim = JuBat.ChooseCell("LG M50")
    param_dim.cell.v_h = 4.3
    
    for model in models
        model_response = Dict()
        
        for freq in frequencies
            println("分析 $model 模型，频率: $freq Hz")
            
            # 设置仿真参数
            opt = JuBat.Option()
            opt.mechanicalmodel = "full"
            opt.model = model
            opt.dtType = "fixed"
            
            # 计算仿真时间和时间步长
            period = 1.0 / freq
            total_time = max(10 * period, 1000)  # 至少10个周期或1000秒
            time_step = period / 100  # 每个周期100个点
            
            opt.dt = time_step
            opt.time = [0, total_time]
            
            # 设置正弦电流输入
            amplitude = 2.0  # 电流幅值 (A)
            bias_current = 3.0  # 偏置电流 (A)
            opt.Current = t -> bias_current + amplitude * sin(2 * π * freq * t)
            
            # 运行仿真
            calc_time = @elapsed begin
                case = JuBat.SetCase(param_dim, opt)
                result = JuBat.Solve(case)
            end
            
            # 提取结果
            time_result = result["time [s]"]
            voltage_result = result["cell voltage [V]"]
            current_input = [bias_current + amplitude * sin(2 * π * freq * t) for t in time_result]
            
            # 只分析稳态响应 (后半段数据)
            start_idx = div(length(time_result), 2)
            steady_time = time_result[start_idx:end] .- time_result[start_idx]
            steady_voltage = voltage_result[start_idx:end]
            steady_current = current_input[start_idx:end]
            
            # 计算幅值和相位响应
            voltage_fft = fft(steady_voltage)
            current_fft = fft(steady_current)
            
            # 找到主频率对应的频率bin
            n_samples = length(steady_voltage)
            freq_bins = fftfreq(n_samples, 1.0 / time_step)
            main_freq_idx = argmin(abs.(freq_bins .- freq))
            
            # 计算传递函数
            if abs(current_fft[main_freq_idx]) > 1e-10
                transfer_function = voltage_fft[main_freq_idx] / current_fft[main_freq_idx]
                magnitude = abs(transfer_function)
                phase = angle(transfer_function) * 180 / π
            else
                magnitude = 0.0
                phase = 0.0
            end
            
            # 计算电压交流成分的幅值
            voltage_ac_amplitude = 2 * abs(voltage_fft[main_freq_idx]) / n_samples
            
            model_response[freq] = Dict(
                "magnitude" => magnitude,
                "phase" => phase,
                "voltage_amplitude" => voltage_ac_amplitude,
                "calc_time" => calc_time,
                "steady_voltage" => steady_voltage,
                "steady_current" => steady_current,
                "steady_time" => steady_time
            )
        end
        
        frequency_response[model] = model_response
    end
    
    return frequency_response
end

# 创建频域可视化图
function create_frequency_plots(frequency_response)
    frequencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    
    # 1. Bode图 - 幅值
    p1 = plot(title="Bode Plot - Magnitude Response",
              xlabel="Frequency [Hz]",
              ylabel="Magnitude [V/A]",
              xscale=:log10,
              yscale=:log10,
              legend=:topright,
              grid=true)
    
    for model in ["P2D", "sP2D"]
        magnitudes = [frequency_response[model][f]["magnitude"] for f in frequencies]
        
        plot!(p1, frequencies, magnitudes,
              label="$model Model",
              marker=:circle,
              linewidth=2,
              markersize=6)
    end
    
    # 2. Bode图 - 相位
    p2 = plot(title="Bode Plot - Phase Response",
              xlabel="Frequency [Hz]",
              ylabel="Phase [degrees]",
              xscale=:log10,
              legend=:bottomleft,
              grid=true)
    
    for model in ["P2D", "sP2D"]
        phases = [frequency_response[model][f]["phase"] for f in frequencies]
        
        plot!(p2, frequencies, phases,
              label="$model Model",
              marker=:diamond,
              linewidth=2,
              markersize=6)
    end
    
    # 3. 电压响应幅值对比
    p3 = plot(title="Voltage Response Amplitude",
              xlabel="Frequency [Hz]",
              ylabel="Voltage Amplitude [V]",
              xscale=:log10,
              legend=:topright,
              grid=true)
    
    for model in ["P2D", "sP2D"]
        v_amplitudes = [frequency_response[model][f]["voltage_amplitude"] for f in frequencies]
        
        plot!(p3, frequencies, v_amplitudes,
              label="$model Model",
              marker=:star,
              linewidth=2,
              markersize=6)
    end
    
    # 4. 相对误差分析
    p4 = plot(title="Frequency Response Error Analysis",
              xlabel="Frequency [Hz]",
              ylabel="Relative Error [%]",
              xscale=:log10,
              legend=:topright,
              grid=true)
    
    # 计算sP2D相对于P2D的误差
    magnitude_errors = []
    phase_errors = []
    
    for freq in frequencies
        p2d_mag = frequency_response["P2D"][freq]["magnitude"]
        sp2d_mag = frequency_response["sP2D"][freq]["magnitude"]
        p2d_phase = frequency_response["P2D"][freq]["phase"]
        sp2d_phase = frequency_response["sP2D"][freq]["phase"]
        
        mag_error = abs(sp2d_mag - p2d_mag) / p2d_mag * 100
        phase_error = abs(sp2d_phase - p2d_phase)
        
        push!(magnitude_errors, mag_error)
        push!(phase_errors, phase_error)
    end
    
    plot!(p4, frequencies, magnitude_errors,
          label="Magnitude Error (%)",
          marker=:circle,
          linewidth=2,
          markersize=6)
    
    plot!(p4, frequencies, phase_errors,
          label="Phase Error (degrees)",
          marker=:triangle,
          linewidth=2,
          markersize=6)
    
    # 组合Bode图
    bode_plot = plot(p1, p2, p3, p4,
                    layout=(2, 2),
                    size=(1200, 800),
                    dpi=300)
    
    savefig(bode_plot, "frequency_domain_analysis.pdf")
    savefig(bode_plot, "frequency_domain_analysis.png")
    
    return bode_plot
end

# 创建时域响应示例图
function create_time_domain_examples(frequency_response)
    # 选择几个代表性频率显示时域响应
    example_freqs = [0.01, 0.1, 1.0, 5.0]
    
    plots_array = []
    
    for freq in example_freqs
        p = plot(title="Time Domain Response at $freq Hz",
                xlabel="Time [s]",
                ylabel="Voltage [V]",
                legend=:topright)
        
        for model in ["P2D", "sP2D"]
            time_data = frequency_response[model][freq]["steady_time"]
            voltage_data = frequency_response[model][freq]["steady_voltage"]
            
            # 只显示前几个周期
            period = 1.0 / freq
            max_time = min(5 * period, maximum(time_data))
            indices = time_data .<= max_time
            
            plot!(p, time_data[indices], voltage_data[indices],
                  label="$model Model",
                  linewidth=2)
        end
        
        push!(plots_array, p)
    end
    
    time_domain_plot = plot(plots_array...,
                           layout=(2, 2),
                           size=(1200, 800),
                           dpi=300)
    
    savefig(time_domain_plot, "time_domain_frequency_examples.pdf")
    savefig(time_domain_plot, "time_domain_frequency_examples.png")
    
    return time_domain_plot
end

# Nyquist图分析
function create_nyquist_plot(frequency_response)
    frequencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    
    p = plot(title="Nyquist Plot",
             xlabel="Real Part [V/A]",
             ylabel="Imaginary Part [V/A]",
             aspect_ratio=:equal,
             legend=:topright,
             grid=true)
    
    for model in ["P2D", "sP2D"]
        real_parts = []
        imag_parts = []
        
        for freq in frequencies
            magnitude = frequency_response[model][freq]["magnitude"]
            phase_rad = frequency_response[model][freq]["phase"] * π / 180
            
            real_part = magnitude * cos(phase_rad)
            imag_part = magnitude * sin(phase_rad)
            
            push!(real_parts, real_part)
            push!(imag_parts, imag_part)
        end
        
        plot!(p, real_parts, imag_parts,
              label="$model Model",
              marker=:circle,
              linewidth=2,
              markersize=6)
        
        # 标记特定频率点
        for (i, freq) in enumerate([0.001, 0.1, 1.0, 10.0])
            if freq in frequencies
                idx = findfirst(x -> x == freq, frequencies)
                annotate!(p, real_parts[idx], imag_parts[idx], 
                         text("$(freq)Hz", 8, :bottom))
            end
        end
    end
    
    savefig(p, "nyquist_plot.pdf")
    savefig(p, "nyquist_plot.png")
    
    return p
end

# 执行频域分析
println("开始频域分析与传递函数验证...")
frequency_response = frequency_domain_analysis()

# 生成可视化
bode_plot = create_frequency_plots(frequency_response)
time_examples = create_time_domain_examples(frequency_response)
nyquist_plot = create_nyquist_plot(frequency_response)

# 保存数据
freq_data = []
for model in ["P2D", "sP2D"]
    for freq in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
        push!(freq_data, Dict(
            "model" => model,
            "frequency" => freq,
            "magnitude" => frequency_response[model][freq]["magnitude"],
            "phase" => frequency_response[model][freq]["phase"],
            "voltage_amplitude" => frequency_response[model][freq]["voltage_amplitude"],
            "calc_time" => frequency_response[model][freq]["calc_time"]
        ))
    end
end

freq_df = DataFrame(freq_data)
CSV.write("frequency_domain_results.csv", freq_df)

# 打印分析结果
println("=== 频域分析结果 ===")
frequencies = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
for freq in [0.01, 0.1, 1.0, 10.0]  # 显示关键频率
    p2d_mag = frequency_response["P2D"][freq]["magnitude"]
    sp2d_mag = frequency_response["sP2D"][freq]["magnitude"]
    p2d_phase = frequency_response["P2D"][freq]["phase"]
    sp2d_phase = frequency_response["sP2D"][freq]["phase"]
    
    mag_error = abs(sp2d_mag - p2d_mag) / p2d_mag * 100
    phase_error = abs(sp2d_phase - p2d_phase)
    
    println("频率 $freq Hz:")
    println("  P2D:  幅值=$(round(p2d_mag, digits=4)), 相位=$(round(p2d_phase, digits=2))°")
    println("  sP2D: 幅值=$(round(sp2d_mag, digits=4)), 相位=$(round(sp2d_phase, digits=2))°")
    println("  误差: 幅值$(round(mag_error, digits=2))%, 相位$(round(phase_error, digits=2))°")
end

println("频域分析完成!") 