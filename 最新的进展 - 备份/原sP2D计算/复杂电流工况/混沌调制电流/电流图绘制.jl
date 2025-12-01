using DifferentialEquations
using Plots

# ======================
# 1. 定义Lorenz系统微分方程
# ======================
function lorenz_system(u, p, t)
    sigma, rho, beta = p
    x, y, z = u
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]
end

# ======================
# 2. 参数设置与求解
# ======================
# Lorenz系统参数 (经典混沌参数)
sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0

# 初始条件和时间设置
u0 = [1.0, 1.0, 1.0]  # 初始状态[x, y, z]
t_span = (0.0, 20.0)  # 时间范围
dt = 0.01             # 时间步长

# 求解Lorenz系统
prob = ODEProblem(lorenz_system, u0, t_span, [sigma, rho, beta])
sol = solve(prob, RK4(), dt=dt)

# 提取解
t = sol.t
x, y, z = sol.u[1,:], sol.u[2,:], sol.u[3,:]

# ======================
# 3. 混沌调制电流生成
# ======================
# 载波信号参数
f_carrier = 50          # 载波频率 (Hz)
carrier_amp = 1.0       # 载波幅度

# 生成载波信号
carrier = carrier_amp * sin.(2 * pi * f_carrier * t)

# 调制过程：用Lorenz系统的x分量调制载波幅度
modulated_current = x .* carrier

# ======================
# 4. 可视化
# ======================
p1 = plot(layout=(3, 1), size=(1500, 1200))

# (1) Lorenz系统x分量（混沌信号）
plot!(t, x, linewidth=0.8, color=:blue, title="Lorenz System Chaotic Signal (x-component)", xlabel="Time (s)", ylabel="Amplitude", subplot=1)

# (2) 载波信号
plot!(t, carrier, linewidth=0.8, color=:green, title="Carrier Signal (50 Hz Sine Wave)", xlabel="Time (s)", ylabel="Amplitude", subplot=2)

# (3) 混沌调制电流
plot!(t, modulated_current, linewidth=0.8, color=:red, title="Lorenz-Modulated Current", xlabel="Time (s)", ylabel="Current (A)", subplot=3)

savefig(p1, "混沌调制电流图.pdf")