using Plots, CSV, DataFrames
include("../src/JuBat.jl") 
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
path = "D:\\竞赛和课程文件\\课程文件\\毕业设计\\src\\data\\"
Crate = 1
i= 5*Crate
opt.Current = x-> i
opt.time = [0 3600]
opt.model = "sP2D" # choose model, other option is "sP2D"

opt.mechanicalmodel = "full"
case1 = JuBat.SetCase(param_dim, opt)
result1 = JuBat.Solve(case1)
        # 转置数据为列向量
        stress_t_vector = vec(result1["negative particle surface tangential stress[Pa]"])
        maxt_df = maximum(stress_t_vector)
        stress_t_df = DataFrame(stress_t = stress_t_vector)

        # 打印最大值
        println("Max tangential stress: ", maxt_df)

        # 保存 sig_tsurf.csv 文件
        file_path = "D:/竞赛和课程文件/课程文件/毕业设计/src/data/sig_tsurf_test1111.csv"
        CSV.write(file_path, stress_t_df)