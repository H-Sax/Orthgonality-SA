using DifferentialEquations,LinearAlgebra, Distributions, Random, GlobalSensitivity, QuasiMonteCarlo, Statistics, LaTeXStrings, Plots, HDF5

function Valve(R, deltaP)
    q = 0.0
    if (-deltaP) < 0.0 
        q =  deltaP/R
    else
        q = 0.0
    end
    return q

end

function ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
    τₑₛ = τₑₛ*τ
    τₑₚ = τₑₚ*τ
    #τ = 4/3(τₑₛ+τₑₚ)
    tᵢ = rem(t + (1 - Eshift) * τ, τ)

    Eₚ = (tᵢ <= τₑₛ) * (1 - cos(tᵢ / τₑₛ * pi)) / 2 +
         (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * (1 + cos((tᵢ - τₑₛ) / (τₑₚ - τₑₛ) * pi)) / 2 +
         (tᵢ <= τₑₚ) * 0

    E = Eₘᵢₙ + (Eₘₐₓ - Eₘᵢₙ) * Eₚ

    return E
end

function DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)

    τₑₛ = τₑₛ*τ
    τₑₚ = τₑₚ*τ
    #τ = 4/3(τₑₛ+τₑₚ)
    tᵢ = rem(t + (1 - Eshift) * τ, τ)

    DEₚ = (tᵢ <= τₑₛ) * pi / τₑₛ * sin(tᵢ / τₑₛ * pi) / 2 +
          (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * pi / (τₑₚ - τₑₛ) * sin((τₑₛ - tᵢ) / (τₑₚ - τₑₛ) * pi) / 2
    (tᵢ <= τₑₚ) * 0
    DE = (Eₘₐₓ - Eₘᵢₙ) * DEₚ

    return DE
end

#Shi timing parameters
Eshift = 0.0
Eₘᵢₙ = 0.03
τₑₛ = 0.3
τₑₚ = 0.45 
Eₘₐₓ = 1.5
Rmv = 0.06
τ = 1.0

function NIK!(du, u, p, t)
    pLV, psa, psv, Vlv, Qav, Qmv, Qs = u 
    τₑₛ, τₑₚ, Rmv, Zao, Rs, Csa, Csv, Eₘₐₓ, Eₘᵢₙ = p
    # pressures (more readable names)
    # the differential equations
    du[1] = (Qmv - Qav) * ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift) + pLV / ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift) * DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
    # 1 Left Ventricle
    du[2] = (Qav - Qs ) / Csa #Systemic arteries     
    du[3] = (Qs - Qmv) / Csv # Venous
    du[4] = Qmv - Qav # volume
    du[5]    = Valve.(Zao, (pLV - psa)) - Qav # AV 
    du[6]   = Valve(Rmv, (psv - pLV)) - Qmv # MV
    du[7]     = (du[2] - du[3]) / Rs # Systemic flow
    nothing 
end
##
M = [1.  0  0  0  0  0  0
     0  1.  0  0  0  0  0
     0  0  1.  0  0  0  0
     0  0  0  1.  0  0  0
     0  0  0  0  0  0  0
     0  0  0  0  0  0  0 
     0  0  0  0  0  0  1. ]
     
Nik_ODE = ODEFunction(NIK!,mass_matrix=M)

u0 = [8.0, 8.0, 8.0, 265.0, 0.0, 0.0, 0.0]

p = [0.3, 0.45, 0.06, 0.033, 1.11, 1.13, 11.0, 1.5, 0.03]

tspan = (0, 20)

prob = ODEProblem(Nik_ODE, u0, tspan, p)
 
saveat = LinRange(15, 16, 300)
@time sol = solve(prob, Rodas4(autodiff = false), reltol = 1e-8, abstol = 1e-8, saveat = saveat)


###### Global Sensitiivty analysis on Discrete outputs ####
# mean(LV pressure), Maximum(SA pressure), Maximum(LV volume) - Minimum(LV Volume) 

circ_time = function (p)
    prob_func(prob,i,repeat) = remake(prob;u0 = [8.0, 8.0, 8.0, 265.0, 0.0, 0.0, 0.0], p=p[:,i])
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
    sol = solve(ensemble_prob,Rodas4(autodiff = false), reltol = 1e-8, abstol = 1e-8,EnsembleThreads();saveat=saveat,trajectories=size(p,2))
    out = zeros(3,size(p,2))
    for i in 1:size(p,2)
      out[1,i] = mean(sol[i][1,:])
      out[2,i] = maximum(sol[i][2,:])
      out[3,i] = maximum(sol[i][4,:])
    end
    out
end

### Sobol Sampler ###

N = 6600
lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
bounds = tuple.(lb,ub)
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)
println("Sobol Jansen")
@time Sobol_JS = gsa(circ_time,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
save("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JS_ST.jld", "data", Sobol_JS.ST)
### Lattice Rule Sampler ###

N = 6600
lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
bounds = tuple.(lb,ub)
sampler = LatticeRuleSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)
println("Lattice Rule Jansen")
@time Sobol_JLR = gsa(circ_time,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
save("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JLR_ST.jld", "data", Sobol_JLR.ST)
### Golden Sampler ###

N = 6600
lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
bounds = tuple.(lb,ub)
sampler = GoldenSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)
println("Golden Jansen")
@time Sobol_JG = gsa(circ_time,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
save("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JG_ST.jld", "data", Sobol_JG.ST)

### Uniform Sampler ###

N = 6600
lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
bounds = tuple.(lb,ub)
sampler = UniformSample()
#A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)
A = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/A_U_Sample.h5", "r") do file
    read(file, "A")
end
B = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/B_U_Sample.h5", "r") do file
    read(file, "A")
end
println("Uniform Jansen")
@time Sobol_JU = gsa(circ_time,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
save("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JU_ST.jld", "data", Sobol_JU.ST)
### Latin Hypercube Sampler ###

N = 6600
lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
bounds = tuple.(lb,ub)
sampler = LatinHypercubeSample()
#A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)
A = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/A_LH_Sample.h5", "r") do file
    read(file, "A")
end
B = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/B_LH_Sample.h5", "r") do file
    read(file, "A")
end
println("Latin Hypercube Jansen")
@time Sobol_JLH = gsa(circ_time,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
save("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JLH_ST.jld", "data", Sobol_JLH.ST)