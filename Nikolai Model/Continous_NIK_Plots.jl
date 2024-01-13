using  JLD, GlobalSensitivity
#### Homma ESTIMATOR ####
# Sobol # 
SH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Homma/Sobol_HC.jld")["data"]
# Lattice Rule # 
LRH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Homma/LatticeRule_HC.jld")["data"]
# Golden #
GH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Homma/Golden_HC.jld")["data"]
# Uniform # 
UH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Homma/Uniform_HC.jld")["data"]
# Latin Hypercube #
LHH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Homma/LatinHypercube_HC.jld")["data"]


#### Sobol ESTIMATOR ####
# Sobol # 
SS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Sobol/Sobol_SC.jld")["data"]
# Lattice Rule # 
LRS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Sobol/LatticeRule_SC.jld")["data"]
# Golden #
GS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Sobol/Golden_SC.jld")["data"]
# Uniform # 
US = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Sobol/Uniform_SC.jld")["data"]
# Latin Hypercube #
LHS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Sobol/LatinHypercube_SC.jld")["data"]


#### JANSEN ESTIMATOR ####
# Sobol # 
SJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Jansen/Sobol_JC.jld")["data"]
# Lattice Rule # 
LRJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Jansen/LatticeRule_JC.jld")["data"]
# Golden #
GJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Jansen/Golden_JC.jld")["data"]
# Uniform # 
UJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Jansen/Uniform_JC.jld")["data"]
# Latin Hypercube #
LHJ = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Jansen/LatinHypercube_JC.jld")["data"]


#### JANON ESTIMATOR ####
# Sobol # 
SJA = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Janon/Sobol_JAC.jld")["data"]
# Lattice Rule # 
LRJA = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Janon/LatticeRule_JAC.jld")["data"]
# Golden #
GJA = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Janon/Golden_JAC.jld")["data"]
# Uniform # 
UJA = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Janon/Uniform_JAC.jld")["data"]
# Latin Hypercube #
LHJA = load("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Continous Data Nik/Janon/LatinHypercube_JAC.jld")["data"]

using CairoMakie
CairoMakie.activate!(type = "svg")
x = LinRange(0,1,300)

## ST - Cao ##
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

    ax = Axis(ga[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{ao}",  limits = (0, 1, 0, 0.7))
    Label(ga[1, 1, Top()], "Homma Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    Label(ga[1, 1, Left()], "Sobol Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = SH.ST[1:300,6] .- 1.96 *SH.ST_Conf_Int[1:300,6]
    higherrors_LVP  = SH.ST[1:300,6] .+ 1.96 *SH.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SH.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SH.ST[301:600,6] .- 1.96 *SH.ST_Conf_Int[301:600,6]
    higherrors_SAP  = SH.ST[301:600,6] .+ 1.96 *SH.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SH.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SH.ST[601:900,6] .- 1.96 *SH.ST_Conf_Int[601:900,6]
    higherrors_LVV  = SH.ST[601:900,6] .+ 1.96 *SH.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SH.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gb[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{ao}",  limits = (0, 1, 0, 0.7))
    Label(gb[1, 1, Left()], "LR Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = LRH.ST[1:300,6] .- 1.96 *LRH.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LRH.ST[1:300,6] .+ 1.96 *LRH.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRH.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRH.ST[301:600,6] .- 1.96 *LRH.ST_Conf_Int[301:600,6]
    higherrors_SAP  = LRH.ST[301:600,6] .+ 1.96 *LRH.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRH.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRH.ST[601:900,6] .- 1.96 *LRH.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LRH.ST[601:900,6] .+ 1.96 *LRH.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRH.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gc[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{ao}",   limits = (0, 1, 0, 0.7))
    Label(gc[1, 1, Left()], "GR Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = GH.ST[1:300,6] .- 1.96 *GH.ST_Conf_Int[1:300,6]
    higherrors_LVP  = GH.ST[1:300,6] .+ 1.96 *GH.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GH.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GH.ST[301:600,6] .- 1.96 *GH.ST_Conf_Int[301:600,6]
    higherrors_SAP  = GH.ST[301:600,6] .+ 1.96 *GH.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GH.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GH.ST[601:900,6] .- 1.96 *GH.ST_Conf_Int[601:900,6]
    higherrors_LVV  = GH.ST[601:900,6] .+ 1.96 *GH.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GH.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gd[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{ao}",  limits = (0, 1, 0, 0.7))
    Label(gd[1, 1, Left()], "Uniform Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = UH.ST[1:300,6] .- 1.96 *UH.ST_Conf_Int[1:300,6]
    higherrors_LVP  = UH.ST[1:300,6] .+ 1.96 *UH.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], UH.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = UH.ST[301:600,6] .- 1.96 *UH.ST_Conf_Int[301:600,6]
    higherrors_SAP  = UH.ST[301:600,6] .+ 1.96 *UH.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], UH.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = UH.ST[601:900,6] .- 1.96 *UH.ST_Conf_Int[601:900,6]
    higherrors_LVV  = UH.ST[601:900,6] .+ 1.96 *UH.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], UH.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(ge[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, ylabel = L"S_{T}~\text{for}~C_{ao}", xlabel = "Cycle Time (s)",  limits = (0, 1, 0, 0.7))
    Label(ge[1, 1, Left()], "LH Sampling", valign = :center, font = :bold, rotation = pi/2, padding = (0, 60, 0, 0))
    #LV.P
    lowerrors_LVP  = LHH.ST[1:300,6] .- 1.96 *LHH.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LHH.ST[1:300,6] .+ 1.96 *LHH.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHH.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHH.ST[301:600,6] .- 1.96 *LHH.ST_Conf_Int[301:600,6]
    higherrors_SAP  = LHH.ST[301:600,6] .+ 1.96 *LHH.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHH.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHH.ST[601:900,6] .- 1.96 *LHH.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LHH.ST[601:900,6] .+ 1.96 *LHH.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHH.ST[601:900,6], label = "LV.V")


    ax = Axis(gf[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    Label(gf[1, 1, Top()], "Sobol Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    #LV.P
    lowerrors_LVP  = SS.ST[1:300,6] .- 1.96 *SS.ST_Conf_Int[1:300,6]
    higherrors_LVP  = SS.ST[1:300,6] .+ 1.96 *SS.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SS.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SS.ST[301:600,6] .- 1.96 *SS.ST_Conf_Int[301:600,6]
    higherrors_SAP  = SS.ST[301:600,6] .+ 1.96 *SS.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SS.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SS.ST[601:900,6] .- 1.96 *SS.ST_Conf_Int[601:900,6]
    higherrors_LVV  = SS.ST[601:900,6] .+ 1.96 *SS.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SS.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gg[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = LRS.ST[1:300,6] .- 1.96 *LRS.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LRS.ST[1:300,6] .+ 1.96 *LRS.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRS.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRS.ST[301:600,6] .- 1.96 *LRS.ST_Conf_Int[301:600,6]
    higherrors_SAP  = LRS.ST[301:600,6] .+ 1.96 *LRS.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRS.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRS.ST[601:900,6] .- 1.96 *LRS.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LRS.ST[601:900,6] .+ 1.96 *LRS.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRS.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gh[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,   limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = GS.ST[1:300,6] .- 1.96 *GS.ST_Conf_Int[1:300,6]
    higherrors_LVP  = GS.ST[1:300,6] .+ 1.96 *GS.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GS.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GS.ST[301:600,6] .- 1.96 *GS.ST_Conf_Int[301:600,6]
    higherrors_SAP  = GS.ST[301:600,6] .+ 1.96 *GS.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GS.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GS.ST[601:900,6] .- 1.96 *GS.ST_Conf_Int[601:900,6]
    higherrors_LVV  = GS.ST[601:900,6] .+ 1.96 *GS.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GS.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gi[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = US.ST[1:300,6] .- 1.96 *US.ST_Conf_Int[1:300,6]
    higherrors_LVP  = US.ST[1:300,6] .+ 1.96 *US.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], US.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = US.ST[301:600,6] .- 1.96 *US.ST_Conf_Int[301:600,6]
    higherrors_SAP  = US.ST[301:600,6] .+ 1.96 *US.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], US.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = US.ST[601:900,6] .- 1.96 *US.ST_Conf_Int[601:900,6]
    higherrors_LVV  = US.ST[601:900,6] .+ 1.96 *US.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], US.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gj[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, xlabel = "Cycle Time (s)",  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = LHS.ST[1:300,6] .- 1.96 *LHS.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LHS.ST[1:300,6] .+ 1.96 *LHS.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHS.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHS.ST[301:600,6] .- 1.96 *LHS.ST_Conf_Int[301:600,6]
    higherrors_SAP  = LHS.ST[301:600,6] .+ 1.96 *LHS.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHS.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHS.ST[601:900,6] .- 1.96 *LHS.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LHS.ST[601:900,6] .+ 1.96 *LHS.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHS.ST[601:900,6], label = "LV.V")

    ax = Axis(gk[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    Label(gk[1, 1, Top()], "Jansen Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    #LV.P
    lowerrors_LVP  = SJ.ST[1:300,6] .- 1.96 *SJ.ST_Conf_Int[1:300,6]
    higherrors_LVP  = SJ.ST[1:300,6] .+ 1.96 *SJ.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SJ.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SJ.ST[301:600,6] .- 1.96 *SJ.ST_Conf_Int[301:600,6]
    higherrors_SAP  = SJ.ST[301:600,6] .+ 1.96 *SJ.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SJ.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SJ.ST[601:900,6] .- 1.96 *SJ.ST_Conf_Int[601:900,6]
    higherrors_LVV  = SJ.ST[601:900,6] .+ 1.96 *SJ.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SJ.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gl[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = LRJ.ST[1:300,6] .- 1.96 *LRJ.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LRJ.ST[1:300,6] .+ 1.96 *LRJ.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRJ.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRJ.ST[301:600,6] .- 1.96 *LRJ.ST_Conf_Int[301:600,6]
    higherrors_SAP  = LRJ.ST[301:600,6] .+ 1.96 *LRJ.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRJ.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRJ.ST[601:900,6] .- 1.96 *LRJ.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LRJ.ST[601:900,6] .+ 1.96 *LRJ.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRJ.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gm[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = GJ.ST[1:300,6] .- 1.96 *GJ.ST_Conf_Int[1:300,6]
    higherrors_LVP  = GJ.ST[1:300,6] .+ 1.96 *GJ.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GJ.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GJ.ST[301:600,6] .- 1.96 *GJ.ST_Conf_Int[301:600,6]
    higherrors_SAP  = GJ.ST[301:600,6] .+ 1.96 *GJ.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GJ.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GJ.ST[601:900,6] .- 1.96 *GJ.ST_Conf_Int[601:900,6]
    higherrors_LVV  = GJ.ST[601:900,6] .+ 1.96 *GJ.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GJ.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gn[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = UJ.ST[1:300,6] .- 1.96 *UJ.ST_Conf_Int[1:300,6]
    higherrors_LVP  = UJ.ST[1:300,6] .+ 1.96 *UJ.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], UJ.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = UJ.ST[301:600,6] .- 1.96 *UJ.ST_Conf_Int[301:600,6]
    higherrors_SAP  = UJ.ST[301:600,6] .+ 1.96 *UJ.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], UJ.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = UJ.ST[601:900,6] .- 1.96 *UJ.ST_Conf_Int[601:900,6]
    higherrors_LVV  = UJ.ST[601:900,6] .+ 1.96 *UJ.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], UJ.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(go[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, xlabel = "Cycle Time (s)", limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = LHJ.ST[1:300,6] .- 1.96 *LHJ.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LHJ.ST[1:300,6] .+ 1.96 *LHJ.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHJ.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHJ.ST[301:600,6] .- 1.96 *LHJ.ST_Conf_Int[301:600,6]
    higherrorss_SAP  = LHJ.ST[301:600,6] .+ 1.96 *LHJ.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrorss_SAP), max.(lowerrors_SAP, higherrorss_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHJ.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHJ.ST[601:900,6] .- 1.96 *LHJ.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LHJ.ST[601:900,6] .+ 1.96 *LHJ.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV , higherrors_LVV ), max.(lowerrors_LVV , higherrors_LVV ), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHJ.ST[601:900,6], label = "LV.V")


    ax = Axis(gp[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    Label(gp[1, 1, Top()], "Janon Estimator", valign = :center, font = :bold, padding = (0, 0, 0, 0))
    #LV.P
    lowerrors_LVP  = SJA.ST[1:300,6] .- 1.96 *SJA.ST_Conf_Int[1:300,6]
    higherrors_LVP  = SJA.ST[1:300,6] .+ 1.96 *SJA.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], SJA.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = SJA.ST[301:600,6] .- 1.96 *SJA.ST_Conf_Int[301:600,6]
    higherrors_SAP  = SJA.ST[301:600,6] .+ 1.96 *SJA.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], SJA.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = SJA.ST[601:900,6] .- 1.96 *SJA.ST_Conf_Int[601:900,6]
    higherrors_LVV  = SJA.ST[601:900,6] .+ 1.96 *SJA.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], SJA.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gq[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = LRJA.ST[1:300,6] .- 1.96 *LRJA.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LRJA.ST[1:300,6] .+ 1.96 *LRJA.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LRJA.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LRJA.ST[301:600,6] .- 1.96 *LRJA.ST_Conf_Int[301:600,6]
    higherrors_SAP  = LRJA.ST[301:600,6] .+ 1.96 *LRJA.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LRJA.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LRJA.ST[601:900,6] .- 1.96 *LRJA.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LRJA.ST[601:900,6] .+ 1.96 *LRJA.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LRJA.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gr[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = GJA.ST[1:300,6] .- 1.96 *GJA.ST_Conf_Int[1:300,6]
    higherrors_LVP  = GJA.ST[1:300,6] .+ 1.96 *GJA.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], GJA.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = GJA.ST[301:600,6] .- 1.96 *GJA.ST_Conf_Int[301:600,6]
    higherrors_SAP  = GJA.ST[301:600,6] .+ 1.96 *GJA.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], GJA.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = GJA.ST[601:900,6] .- 1.96 *GJA.ST_Conf_Int[601:900,6]
    higherrors_LVV  = GJA.ST[601:900,6] .+ 1.96 *GJA.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], GJA.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gs[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5,  limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = UJA.ST[1:300,6] .- 1.96 *UJA.ST_Conf_Int[1:300,6]
    higherrors_LVP  = UJA.ST[1:300,6] .+ 1.96 *UJA.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], UJA.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = UJA.ST[301:600,6] .- 1.96 *UJA.ST_Conf_Int[301:600,6]
    higherrors_SAP  = UJA.ST[301:600,6] .+ 1.96 *UJA.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrors_SAP), max.(lowerrors_SAP, higherrors_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], UJA.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = UJA.ST[601:900,6] .- 1.96 *UJA.ST_Conf_Int[601:900,6]
    higherrors_LVV  = UJA.ST[601:900,6] .+ 1.96 *UJA.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV, higherrors_LVV), max.(lowerrors_LVV, higherrors_LVV), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], UJA.ST[601:900,6], label = "LV.V")
    hidexdecorations!(ax, grid = false)

    ax = Axis(gt[1,1], xgridstyle = :dash, ygridstyle = :dash,xticksize = 0.5, yticks = [0.0, 0.2, 0.4, 0.6],  yticksize = 5, xlabel = "Cycle Time (s)", limits = (0, 1, 0, 0.7))
    #LV.P
    lowerrors_LVP  = LHJA.ST[1:300,6] .- 1.96 *LHJA.ST_Conf_Int[1:300,6]
    higherrors_LVP  = LHJA.ST[1:300,6] .+ 1.96 *LHJA.ST_Conf_Int[1:300,6]
    band!(x[1:end], min.(lowerrors_LVP, higherrors_LVP), max.(lowerrors_LVP, higherrors_LVP), color=(:blue, 0.5),  transparency=true)
    lines!(x[1:end], LHJA.ST[1:300,6], label = "LV.P")
    #SAT.P
    lowerrors_SAP  = LHJA.ST[301:600,6] .- 1.96 *LHJA.ST_Conf_Int[301:600,6]
    higherrorss_SAP  = LHJA.ST[301:600,6] .+ 1.96 *LHJA.ST_Conf_Int[301:600,6]
    band!(x[1:end], min.(lowerrors_SAP, higherrorss_SAP), max.(lowerrors_SAP, higherrorss_SAP), color=(:yellow, 0.5),  transparency=true)
    lines!(x[1:end], LHJA.ST[301:600,6], label = "SAT.P")
    #LV.V
    lowerrors_LVV  = LHJA.ST[601:900,6] .- 1.96 *LHJA.ST_Conf_Int[601:900,6]
    higherrors_LVV  = LHJA.ST[601:900,6] .+ 1.96 *LHJA.ST_Conf_Int[601:900,6]
    band!(x[1:end], min.(lowerrors_LVV , higherrors_LVV ), max.(lowerrors_LVV , higherrors_LVV ), color=(:green, 0.5),  transparency=true)
    lines!(x[1:end], LHJA.ST[601:900,6], label = "LV.V")


    f[6, 1:4] = Legend(f, ax, orientation = :horizontal)

    for (label, layout) in zip(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"], [ga, gb, gc, gd, ge, gf, gg, gh, gi, gj, gk, gl, gm, gn, go, gp, gq, gr, gs, gt])
        Label(layout[1, 1, TopRight()], label,fontsize = 18,font = :bold,halign = :right)
    end

    resize_to_layout!(f)

    f
end 