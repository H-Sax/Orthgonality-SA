using Pkg
Pkg.add("Dates")
using Dates
println("Running on ", Threads.nthreads(), " threads.")
println(now())

Pkg.add("DifferentialEquations")
Pkg.add("Statistics")
Pkg.add("Distributions")
Pkg.add("QuasiMonteCarlo")
Pkg.add("JLD")
@time "Setup modules" using DifferentialEquations, Statistics, Distributions, QuasiMonteCarlo, JLD, GlobalSensitivity

function Valve(R, deltaP, open)
    dq = 0.0
    if (-open) < 0.0 
        dq =  deltaP/R
    else
        dq = 0.0
    end
    return dq

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
    du[5]    = Valve(Zao, (du[1] - du[2]), u[1] - u[2])  # AV 
    du[6]   = Valve(Rmv, (du[3] - du[1]), u[3] - u[1])  # MV
    du[7]     = (du[2] - du[3]) / Rs # Systemic flow
end
##

u0 = [8.0, 8.0, 8.0, 265.0, 0.0, 0.0, 0.0]

p = [0.3, 0.45, 0.06, 0.033, 1.11, 1.13, 11.0, 1.5, 0.03]

tspan = (0, 10)

prob = ODEProblem(NIK!, u0, tspan, p)

x = LinRange(7,8,100)
@time sol = solve(prob, Vern7(),  reltol = 1e-5, abstol = 1e-5, saveat = x)

circ_time_post = function (p)
    out_global
end

###### Global Sensitiivty analysis on continous outputs ####
# LV pressure, SA pressure, LV volume 

circ_time_idxs_chunked = function (p)

    global parameterset = p

    Np::Int64 = size(p,2)

    println("Running ", Np, " cases.")

    chunksize::Int64 = Np/chunks

    println("Using ", chunks, " chunk(s) of size ", chunksize, ".")


    out = zeros(3,Np)

    for k in 1:chunks
        offset = Int((k - 1) * Np/chunks)
        startindx = Int(offset + 1)
        endindx = Int(k * Np/chunks)

        println("Starting chunk ", k, ", from ", startindx, " to ", endindx, " with offset ", offset, ".")

        pchunk = p[:,startindx:endindx]

        prob_func(prob,i,repeat) = remake(prob; u0 = [8.0, 8.0, 8.0, 265.0, 0.0, 0.0, 0.0], p=pchunk[:,i])

        ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)

        @time "    ODE Ensemble solved in " sol = solve(ensemble_prob, Vern7(),  reltol = 1e-5, abstol = 1e-5, EnsembleThreads();saveat = x,trajectories=chunksize)
        @time "    Results processed in" Threads.@threads for i in 1:chunksize
          # println(i+offset)
          out[1,i+offset] = mean(sol[i][1,:]) #LV.P
          out[2,i+offset] = maximum(sol[i][2,:]) #SA.P
          out[3,i+offset] = maximum(sol[i][4,:]) #LV.V
        end
    end
    global out_global = out
    out
end

## 100k test ##
println("Sobol Estimator - Sobol")
#save("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Discrete/Discrete Data/Sobol_HS.jld", "data", Sobol_HS)
chunks::Int64 = 10000
samples = 100000
lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
@time Sobol_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
save("Nik_Sobol_S.jld", "data", Sobol_S)


###### Sobol Sampler #####
#
#println("Jansen Estimator - Sobol")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = Sobol_J.ST
#    b2[:,:,i] = Sobol_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Sobol_J_ST.jld", "data", b1)
#save("Sobol_J_ST_conf.jld", "data", b2)
#
#println("Sobol Estimator - Sobol")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = Sobol_S.ST
#    b2[:,:,i] = Sobol_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Sobol_S_ST.jld", "data", b1)
#save("Sobol_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - Sobol")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = Sobol_H.ST
#    b2[:,:,i] = Sobol_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Sobol_H_ST.jld", "data", b1)
#save("DN_Sobol_H_ST_conf.jld", "data", b2)
#
#println("Janon Estimator - Sobol")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = Sobol_JA.ST
#    b2[:,:,i] = Sobol_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Sobol_JA_ST.jld", "data", b1)
#save("DN_Sobol_JA_ST_conf.jld", "data", b2)

###### Lattice Rule Sampling ########
#
#println("Jansen Estimator - Lattice Rule")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = lattice_J.ST
#    b2[:,:,i] = lattice_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Lattice_J_ST.jld", "data", b1)
#save("Lattice_J_ST_conf.jld", "data", b2)
#
#println("Sobol Estimator - Lattice Rule")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = lattice_S.ST
#    b2[:,:,i] = lattice_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Lattice_S_ST.jld", "data", b1)
#save("Lattice_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - Lattice Rule")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = lattice_H.ST
#    b2[:,:,i] = lattice_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Lattice_H_ST.jld", "data", b1)
#save("DN_Lattice_H_ST_conf.jld", "data", b2)
#
#
#println("Janon Estimator - Lattice Rule")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = lattice_JA.ST
#    b2[:,:,i] = lattice_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Lattice_JA_ST.jld", "data", b1)
#save("DN_Lattice_JA_ST_conf.jld", "data", b2)

#### Golden Sampling ###
#
#println("Jansen Estimator - Golden")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = Golden_J.ST
#    b2[:,:,i] = Golden_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Golden_J_ST.jld", "data", b1)
#save("Golden_J_ST_conf.jld", "data", b2)
#
#println("Sobol Estimator - Golden")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = Golden_S.ST
#    b2[:,:,i] = Golden_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Golden_S_ST.jld", "data", b1)
#save("Golden_S_ST_conf.jld", "data", b2)


#println("Homma Estimator - Golden")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = Golden_H.ST
#    b2[:,:,i] = Golden_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Golden_H_ST.jld", "data", b1)
#save("DN_Golden_H_ST_conf.jld", "data", b2)
#
#println("Janon Estimator - Golden")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = Golden_JA.ST
#    b2[:,:,i] = Golden_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Golden_JA_ST.jld", "data", b1)
#save("DN_Golden_JA_ST_conf.jld", "data", b2)

## Uniform Sampling 

#println("Jansen Estimator - Uniform")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = Uniform_J.ST
#    b2[:,:,i] = Uniform_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Uniform_J_ST.jld", "data", b1)
#save("Uniform_J_ST_conf.jld", "data", b2)
#
#
#println("Sobol Estimator - Uniform")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = Uniform_S.ST
#    b2[:,:,i] = Uniform_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Uniform_S_ST.jld", "data", b1)
#save("Uniform_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - Uniform")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = Uniform_H.ST
#    b2[:,:,i] = Uniform_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Uniform_H_ST.jld", "data", b1)
#save("DN_Uniform_H_ST_conf.jld", "data", b2)
#
#
#println("Janon Estimator - Uniform")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = Uniform_JA.ST
#    b2[:,:,i] = Uniform_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_Uniform_JA_ST.jld", "data", b1)
#save("DN_Uniform_JA_ST_conf.jld", "data", b2)



### Latin HyperCube ###

#println("Jansen Estimator - LH")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = LH_J.ST
#    b2[:,:,i] = LH_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("LH_J_ST.jld", "data", b1)
#save("LH_J_ST_conf.jld", "data", b2)


#println("Sobol Estimator - LH")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = LH_S.ST
#    b2[:,:,i] = LH_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("LH_S_ST.jld", "data", b1)
#save("LH_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - LH")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = LH_H.ST
#    b2[:,:,i] = LH_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_LH_H_ST.jld", "data", b1)
#save("DN_LH_H_ST_conf.jld", "data", b2)
#
#
#println("Janon Estimator - LH")
#b1 = zeros(3,9,20)
#b2 = zeros(3,9,20)
#N = range(start = 2000, stop = 40000, step = 2000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [0.21, 0.36, 0.042, 0.0231, 0.777, 0.791, 7.7, 1.05, 0.021]
#    ub = [0.34, 0.585, 0.078, 0.0429, 1.443, 1.469, 14.3, 1.95, 0.039]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = LH_JA.ST
#    b2[:,:,i] = LH_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("DN_LH_JA_ST.jld", "data", b1)
#save("DN_LH_JA_ST_conf.jld", "data", b2)