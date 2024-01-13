using  JLD, GlobalSensitivity

#### Homma ESTIMATOR ####
# Sobol # 
SH =  load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Sobol/Sobol_HD.jld")["data"]
SH_E =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Sobol/Sobol_HD_conf.jld")["data"]
# Lattice Rule # 
LRH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Lattice Rule/LatticeRule_HD.jld")["data"]
LRH_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Lattice Rule/LatticeRule_HD_conf.jld")["data"]
# Golden #
GH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Golden/Golden_HD.jld")["data"]
GH_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Golden/Golden_HD_conf.jld")["data"]
# Uniform # 
UH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Uniform/Uniform_HD.jld")["data"]
UH_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Uniform/Uniform_HD_conf.jld")["data"]
# Latin Hypercube #
LHH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Latin hypercube/LatinHypercube_HD.jld")["data"]
LHH_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Homma/Latin hypercube/LatinHypercube_HD_conf.jld")["data"]

#### Sobol ESTIMATOR ####
# Sobol # 
SS =  load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Sobol/Sobol_SD.jld")["data"]
SS_E =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Sobol/Sobol_SD_conf.jld")["data"]
# Lattice Rule # 
LRS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Lattice Rule/LatticeRule_SD.jld")["data"]
LRS_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Lattice Rule/LatticeRule_SD_conf.jld")["data"]
# Golden #
GS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Golden/Golden_SD.jld")["data"]
GS_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Golden/Golden_SD_conf.jld")["data"]
# Uniform # 
US = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Uniform/Uniform_SD.jld")["data"]
US_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Uniform/Uniform_SD_conf.jld")["data"]
# Latin Hypercube #
LHS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Latin hypercube/LatinHypercube_SD.jld")["data"]
LHS_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Sobol/Latin hypercube/LatinHypercube_SD_conf.jld")["data"]

#### JANSEN ESTIMATOR ####
# Sobol # 
SJ =  load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Sobol/Sobol_JD.jld")["data"]
SJ_E =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Sobol/Sobol_JD_conf.jld")["data"]
# Lattice Rule # 
LRJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Lattice Rule/LatticeRule_JD.jld")["data"]
LRJ_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Lattice Rule/LatticeRule_JD_conf.jld")["data"]
# Golden #
GJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Golden/Golden_JD.jld")["data"]
GJ_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Golden/Golden_JD_conf.jld")["data"]
# Uniform # 
UJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Uniform/Uniform_JD.jld")["data"]
UJ_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Uniform/Uniform_JD_conf.jld")["data"]
# Latin Hypercube #
LHJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Latin hypercube/LatinHypercube_JD.jld")["data"]
LHJ_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Jansen/Latin hypercube/LatinHypercube_JD_conf.jld")["data"]

#### JANON ESTIMATOR ####
# Sobol # 
SJA = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Sobol/Sobol_JAD.jld")["data"]
SJA_E =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Sobol/Sobol_JAD_conf.jld")["data"]
# Lattice Rule # 
LRJA =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Lattice Rule/LatticeRule_JAD.jld")["data"]
LRJA_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Lattice Rule/LatticeRule_JAD_conf.jld")["data"]
# Golden #
GJA =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Golden/Golden_JAD.jld")["data"]
GJA_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Golden/Golden_JAD_conf.jld")["data"]
# Uniform # 
UJA =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Uniform/Uniform_JAD.jld")["data"]
UJA_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Uniform/Uniform_JAD_conf.jld")["data"]
# Latin Hypercube #
LHJA =load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Latin hypercube/LatinHypercube_JAD.jld")["data"]
LHJA_E=load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/Discrete data/Janon/Latin hypercube/LatinHypercube_JAD_conf.jld")["data"]


## 100k Sobol estimator data ##
test = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi BS and convergence/100k_sobol_dis.jld")["data"]




using CairoMakie
CairoMakie.activate!(type = "svg")
x = LinRange(10000, 30000, 5)
## ST - Csvn ##
begin
    f = Figure(resolution = (1400, 1200), backgroundcolor = RGBf(0.98, 0.98, 0.98))
    #f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98));
    gabcde = f[1:5,1]
    ga = gabcde[1, 1] = GridLayout()
    gb = gabcde[2, 1] = GridLayout()
    gc = gabcde[3, 1] = GridLayout()
    gd = gabcde[4, 1] = GridLayout()
    ge = gabcde[5, 1] = GridLayout()

    gfghij = f[1:5, 2] = GridLayout()
    gf = gfghij[1, 1] = GridLayout()
    gg = gfghij[2, 1] = GridLayout()
    gh = gfghij[3, 1] = GridLayout()
    gi = gfghij[4, 1] = GridLayout()
    gj = gfghij[5, 1] = GridLayout()

    gklmno = f[1:5, 3] = GridLayout()
    gk = gklmno[1, 1] = GridLayout()
    gl = gklmno[2, 1] = GridLayout()
    gm = gklmno[3, 1] = GridLayout()
    gn = gklmno[4, 1] = GridLayout()
    go = gklmno[5, 1] = GridLayout()

    gpqrst = f[1:5, 4] = GridLayout()
    gp = gpqrst[1, 1] = GridLayout()
    gq = gpqrst[2, 1] = GridLayout()
    gr = gpqrst[3, 1] = GridLayout()
    gs = gpqrst[4, 1] = GridLayout()
    gt = gpqrst[5, 1] = GridLayout()


    ax = Axis(ga[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{svn}", limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    Label(ga[1, 1, Top()], "Homma Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    Label(ga[1, 1, Left()], "Sobol Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = SH[1,31,:] .- 1.96 *SH_E[1,31,:]
    higherrors_LVP  = SH[1,31,:] .+ 1.96 *SH_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SH[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SH[2,31,:] .- 1.96 *SH_E[2,31,:]
    higherrors_SAP  = SH[2,31,:] .+ 1.96 *SH_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SH[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SH[3,31,:] .- 1.96 *SH_E[3,31,:]
    higherrors_LVV  = SH[3,31,:] .+ 1.96 *SH_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SH[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gb[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{svn}",  limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    Label(gb[1, 1, Left()], "LR Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = LRH[1,31,:] .- 1.96 *LRH_E[1,31,:]
    higherrors_LVP  = LRH[1,31,:] .+ 1.96 *LRH_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRH[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRH[2,31,:] .- 1.96 *LRH_E[2,31,:]
    higherrors_SAP  = LRH[2,31,:] .+ 1.96 *LRH_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRH[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRH[3,31,:] .- 1.96 *LRH_E[3,31,:]
    higherrors_LVV  = LRH[3,31,:] .+ 1.96 *LRH_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRH[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gc[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{svn}",  limits = (10000, 30000, 0.4, 0.8), xlabel = "Sample Size", xticks = [10000, 20000, 30000])
    Label(gc[1, 1, Left()], "GR Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = GH[1,31,:] .- 1.96 *GH_E[1,31,:]
    higherrors_LVP  = GH[1,31,:] .+ 1.96 *GH_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GH[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GH[2,31,:] .- 1.96 *GH_E[2,31,:]
    higherrors_SAP  = GH[2,31,:] .+ 1.96 *GH_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GH[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GH[3,31,:] .- 1.96 *GH_E[3,31,:]
    higherrors_LVV  = GH[3,31,:] .+ 1.96 *GH_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GH[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gd[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,   yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{svn}",  limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    Label(gd[1, 1, Left()], "Uniform Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = UH[1,31,:] .- 1.96 *UH_E[1,31,:]
    higherrors_LVP  = UH[1,31,:] .+ 1.96 *UH_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], UH[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = UH[2,31,:] .- 1.96 *UH_E[2,31,:]
    higherrors_SAP  = UH[2,31,:] .+ 1.96 *UH_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], UH[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = UH[3,31,:] .- 1.96 *UH_E[3,31,:]
    higherrors_LVV  = UH[3,31,:] .+ 1.96 *UH_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], UH[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(ge[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{svn}", xlabel = "Sample Size",  limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    Label(ge[1, 1, Left()], "LH Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = LHH[1,31,:] .- 1.96 *LHH_E[1,31,:]
    higherrors_LVP  = LHH[1,31,:] .+ 1.96 *LHH_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHH[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHH[2,31,:] .- 1.96 *LHH_E[2,31,:]
    higherrors_SAP  = LHH[2,31,:] .+ 1.96 *LHH_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHH[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHH[3,31,:] .- 1.96 *LHH_E[3,31,:]
    higherrors_LVV  = LHH[3,31,:] .+ 1.96 *LHH_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHH[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)


    ax = Axis(gf[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    Label(gf[1, 1, Top()], "Sobol Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    #LV.P
    lowerrors_LVP  = SS[1,31,:] .- 1.96 *SS_E[1,31,:]
    higherrors_LVP  = SS[1,31,:] .+ 1.96 *SS_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SS[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SS[2,31,:] .- 1.96 *SS_E[2,31,:]
    higherrors_SAP  = SS[2,31,:] .+ 1.96 *SS_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SS[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SS[3,31,:] .- 1.96 *SS_E[3,31,:]
    higherrors_LVV  = SS[3,31,:] .+ 1.96 *SS_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SS[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gg[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,   yticksize = 5, limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    #LV.P
    lowerrors_LVP  = LRS[1,31,:] .- 1.96 *LRS_E[1,31,:]
    higherrors_LVP  = LRS[1,31,:] .+ 1.96 *LRS_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRS[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRS[2,31,:] .- 1.96 *LRS_E[2,31,:]
    higherrors_SAP  = LRS[2,31,:] .+ 1.96 *LRS_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRS[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRS[3,31,:] .- 1.96 *LRS_E[3,31,:]
    higherrors_LVV  = LRS[3,31,:] .+ 1.96 *LRS_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRS[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gh[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    #LV.P
    lowerrors_LVP  = GS[1,31,:] .- 1.96 *GS_E[1,31,:]
    higherrors_LVP  = GS[1,31,:] .+ 1.96 *GS_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GS[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GS[2,31,:] .- 1.96 *GS_E[2,31,:]
    higherrors_SAP  = GS[2,31,:] .+ 1.96 *GS_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GS[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GS[3,31,:] .- 1.96 *GS_E[3,31,:]
    higherrors_LVV  = GS[3,31,:] .+ 1.96 *GS_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GS[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gi[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticksize = 5, limits = (10000, 30000, 0.4, 0.8), xticks = [10000, 20000, 30000])
    #LV.P
    lowerrors_LVP  = US[1,31,:] .- 1.96 *US_E[1,31,:]
    higherrors_LVP  = US[1,31,:] .+ 1.96 *US_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], US[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = US[2,31,:] .- 1.96 *US_E[2,31,:]
    higherrors_SAP  = US[2,31,:] .+ 1.96 *US_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], US[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = US[3,31,:] .- 1.96 *US_E[3,31,:]
    higherrors_LVV  = US[3,31,:] .+ 1.96 *US_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], US[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gj[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xticks = [10000, 20000, 30000] , limits = (10000, 30000, 0.4, 0.8), xlabel = "Sample Size")
    #LV.P
    lowerrors_LVP  = LHS[1,31,:] .- 1.96 *LHS_E[1,31,:]
    higherrors_LVP  = LHS[1,31,:] .+ 1.96 *LHS_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHS[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHS[2,31,:] .- 1.96 *LHS_E[2,31,:]
    higherrorss_SAP  = LHS[2,31,:] .+ 1.96 *LHS_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrorss_SAP), max.(lowerrors_SAP, higherrorss_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHS[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHS[3,31,:] .- 1.96 *LHS_E[3,31,:]
    higherrors_LVV  = LHS[3,31,:] .+ 1.96 *LHS_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV , higherrors_LVV ), max.(lowerrors_LVV , higherrors_LVV ), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHS[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)

    ax = Axis(gk[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    Label(gk[1, 1, Top()], "Jansen Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    #LV.P
    lowerrors_LVP  = SJ[1,31,:] .- 1.96 *SJ_E[1,31,:]
    higherrors_LVP  = SJ[1,31,:] .+ 1.96 *SJ_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SJ[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SJ[2,31,:] .- 1.96 *SJ_E[2,31,:]
    higherrors_SAP  = SJ[2,31,:] .+ 1.96 *SJ_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SJ[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SJ[3,31,:] .- 1.96 *SJ_E[3,31,:]
    higherrors_LVV  = SJ[3,31,:] .+ 1.96 *SJ_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SJ[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gl[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = LRJ[1,31,:] .- 1.96 *LRJ_E[1,31,:]
    higherrors_LVP  = LRJ[1,31,:] .+ 1.96 *LRJ_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRJ[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRJ[2,31,:] .- 1.96 *LRJ_E[2,31,:]
    higherrors_SAP  = LRJ[2,31,:] .+ 1.96 *LRJ_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRJ[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRJ[3,31,:] .- 1.96 *LRJ_E[3,31,:]
    higherrors_LVV  = LRJ[3,31,:] .+ 1.96 *LRJ_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRJ[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gm[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xlabel = "Sample Size", xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = GJ[1,31,:] .- 1.96 *GJ_E[1,31,:]
    higherrors_LVP  = GJ[1,31,:] .+ 1.96 *GJ_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GJ[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GJ[2,31,:] .- 1.96 *GJ_E[2,31,:]
    higherrors_SAP  = GJ[2,31,:] .+ 1.96 *GJ_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GJ[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GJ[3,31,:] .- 1.96 *GJ_E[3,31,:]
    higherrors_LVV  = GJ[3,31,:] .+ 1.96 *GJ_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GJ[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gn[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,   yticksize = 5, xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = UJ[1,31,:] .- 1.96 *UJ_E[1,31,:]
    higherrors_LVP  = UJ[1,31,:] .+ 1.96 *UJ_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], UJ[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = UJ[2,31,:] .- 1.96 *UJ_E[2,31,:]
    higherrors_SAP  = UJ[2,31,:] .+ 1.96 *UJ_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], UJ[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = UJ[3,31,:] .- 1.96 *UJ_E[3,31,:]
    higherrors_LVV  = UJ[3,31,:] .+ 1.96 *UJ_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], UJ[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(go[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,   yticksize = 5, xlabel = "Sample Size", xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = LHJ[1,31,:] .- 1.96 *LHJ_E[1,31,:]
    higherrors_LVP  = LHJ[1,31,:] .+ 1.96 *LHJ_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHJ[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHJ[2,31,:] .- 1.96 *LHJ_E[2,31,:]
    higherrors_SAP  = LHJ[2,31,:] .+ 1.96 *LHJ_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHJ[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHJ[3,31,:] .- 1.96 *LHJ_E[3,31,:]
    higherrors_LVV  = LHJ[3,31,:] .+ 1.96 *LHJ_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHJ[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)

    ax = Axis(gp[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    Label(gp[1, 1, Top()], "Janon Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    #LV.P
    lowerrors_LVP  = SJA[1,31,:] .- 1.96 *SJA_E[1,31,:]
    higherrors_LVP  = SJA[1,31,:] .+ 1.96 *SJA_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SJA[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SJA[2,31,:] .- 1.96 *SJA_E[2,31,:]
    higherrors_SAP  = SJA[2,31,:] .+ 1.96 *SJA_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SJA[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SJA[3,31,:] .- 1.96 *SJA_E[3,31,:]
    higherrors_LVV  = SJA[3,31,:] .+ 1.96 *SJA_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SJA[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gq[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = LRJA[1,31,:] .- 1.96 *LRJA_E[1,31,:]
    higherrors_LVP  = LRJA[1,31,:] .+ 1.96 *LRJA_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRJA[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRJA[2,31,:] .- 1.96 *LRJA_E[2,31,:]
    higherrors_SAP  = LRJA[2,31,:] .+ 1.96 *LRJA_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRJA[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRJA[3,31,:] .- 1.96 *LRJA_E[3,31,:]
    higherrors_LVV  = LRJA[3,31,:] .+ 1.96 *LRJA_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRJA[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gr[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,  yticksize = 5, xlabel = "Sample Size", xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = GJA[1,31,:] .- 1.96 *GJA_E[1,31,:]
    higherrors_LVP  = GJA[1,31,:] .+ 1.96 *GJA_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GJA[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GJA[2,31,:] .- 1.96 *GJA_E[2,31,:]
    higherrors_SAP  = GJA[2,31,:] .+ 1.96 *GJA_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GJA[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GJA[3,31,:] .- 1.96 *GJA_E[3,31,:]
    higherrors_LVV  = GJA[3,31,:] .+ 1.96 *GJA_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GJA[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gs[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,   yticksize = 5, xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = UJA[1,31,:] .- 1.96 *UJA_E[1,31,:]
    higherrors_LVP  = UJA[1,31,:] .+ 1.96 *UJA_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], UJA[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = UJA[2,31,:] .- 1.96 *UJA_E[2,31,:]
    higherrors_SAP  = UJA[2,31,:] .+ 1.96 *UJA_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], UJA[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = UJA[3,31,:] .- 1.96 *UJA_E[3,31,:]
    higherrors_LVV  = UJA[3,31,:] .+ 1.96 *UJA_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], UJA[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)
    hidexdecorations!(ax, grid = false)

    ax = Axis(gt[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5,   yticksize = 5, xlabel = "Sample Size", xticks = [10000, 20000, 30000], limits = (10000, 30000, 0.4, 0.8))
    #LV.P
    lowerrors_LVP  = LHJA[1,31,:] .- 1.96 *LHJA_E[1,31,:]
    higherrors_LVP  = LHJA[1,31,:] .+ 1.96 *LHJA_E[1,31,:]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHJA[1,31,:], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHJA[2,31,:] .- 1.96 *LHJA_E[2,31,:]
    higherrors_SAP  = LHJA[2,31,:] .+ 1.96 *LHJA_E[2,31,:]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHJA[2,31,:], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHJA[3,31,:] .- 1.96 *LHJA_E[3,31,:]
    higherrors_LVV  = LHJA[3,31,:] .+ 1.96 *LHJA_E[3,31,:]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHJA[3,31,:], label = "LV.V")
    vlines!([20000], label = "Results Taken", color = :red, linestyle = :dashdot)

    f[6, 1:4] = Legend(f, ax, orientation = :horizontal)

    for (label, layout) in zip(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"], [ga, gb, gc, gd, ge, gf, gg, gh, gi, gj, gk, gl, gm, gn, go, gp, gq, gr, gs, gt])
        Label(layout[1, 1, TopRight()], label,fontsize = 18,font = :bold,halign = :right)
    end

    resize_to_layout!(f)

    f
end 





