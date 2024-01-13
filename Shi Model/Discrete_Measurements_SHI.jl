using Pkg
Pkg.add("Dates")
using Dates
println("Running on ", Threads.nthreads(), " threads.")
println(now())

Pkg.add("DifferentialEquations")
Pkg.add("Statistics")
Pkg.add("CirculatorySystemModels")
Pkg.add("Distributions")
Pkg.add("QuasiMonteCarlo")
Pkg.add("JLD")
Pkg.add("ModelingToolkit")
Pkg.add("GlobalSensitivity")
@time "Setup modules" using ModelingToolkit, DifferentialEquations, Statistics, CirculatorySystemModels,  Distributions, QuasiMonteCarlo, JLD, GlobalSensitivity
##
# Shi parameters for chambers and circulation that will not change 
begin 
  τ::Float64 = 1.0
  v0_lv::Float64 = 5.0
  p0_lv::Float64 = 1.0
  Emin_lv::Float64 = 0.1
  Emax_lv::Float64 = 2.5
  τes_lv::Float64 = 0.3
  τed_lv::Float64 = 0.45
  Eshift_lv::Float64 = 0.0
  LV_Vt0::Float64 = 500
  ### LA Atrium Parameters #### Checked
  v0_la::Float64 = 4.0
  p0_la::Float64 = 1.0
  Emin_la::Float64 = 0.15
  Emax_la::Float64 = 0.25
  τpwb_la::Float64 = 0.92
  τpww_la::Float64 = 0.09
  τes_la::Float64 = τpww_la/2
  τed_la::Float64 = τpww_la
  Eshift_la::Float64 = τpwb_la
  LA_Vt0::Float64 = 20
  #####
  Csas::Float64 = 0.08
  Rsas::Float64 = 0.003
  Lsas::Float64 = 6.2e-5
  pt0sas::Float64 = 100.0
  qt0sas::Float64 = 0.0
  ## Systemic Artery #### Checked
  Csat::Float64 = 1.6
  Rsat::Float64 = 0.05
  Lsat::Float64 = 0.0017
  pt0sat::Float64 = 100.0
  qt0sat::Float64 = 0.0
  ## Systemic Arteriole #### Checked
  Rsar::Float64 = 0.5
  ## Systemic Capillary #### Checked 
  Rscp::Float64 = 0.52
  ## Systemic Vein #### Checked
  Csvn::Float64 = 20.5
  Rsvn::Float64 = 0.075
  pt0svn::Float64 = 0.0
  qt0svn::Float64 = 0.0
  ## Valve Parameters 
  Zao::Float64 = 0.033
  Rmv::Float64 = 0.06
end

@time "Setup Model" begin
@parameters t
@named LV = ShiChamber(V₀=v0_lv, p₀ = p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
@named LA = ShiChamber(V₀=v0_la, p₀ = p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la/2, τₑₚ=τpww_la, Eshift=τpwb_la)

@named AV = ResistorDiode(R = Zao)
@named MV = ResistorDiode(R = Rmv)

@named Sys_loop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas, SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat, SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn)

circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Sys_loop.in)
    connect(Sys_loop.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

# Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)


@named circ_model = compose(_circ_model,[LV, LA, AV, MV, Sys_loop ])
end

## And simplify it
@time "Simplify model"  circ_sys = structural_simplify(circ_model)

u0 =  [LV_Vt0, LA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn]
@time "Define problem" prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

saveat = LinRange(15, 16, 150)
@time "Solve first time" sol = solve(prob, Tsit5(),  reltol=1e-4, abstol=1e-4, saveat = saveat)
@time "Solve second time" sol = solve(prob, Tsit5(),  reltol=1e-4, abstol=1e-4, saveat = saveat)



LV_p = function (LV_V,tsol,p)
    v0_lv, p0_lv, Emin_lv, Emax_lv, τ, τes_lv, τed_lv, Eshift_lv = p[1:8]
    @. p0_lv + (Emin_lv + (Emax_lv - Emin_lv)*((1//2)*(1 - cos((π*rem(tsol + τ*(1 - Eshift_lv), τ)) / τes_lv)) *(rem(tsol + τ*(1 - Eshift_lv), τ) <= τes_lv) + (1//2)*(1 + cos((π*(rem(tsol + τ*(1 - Eshift_lv), τ) - τes_lv)) / (τed_lv - τes_lv)))*(rem(tsol + τ*(1 - Eshift_lv), τ) <= τed_lv)*(rem(tsol + τ*(1 - Eshift_lv), τ) > τes_lv)))*(LV_V - v0_lv)
end


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

        prob_func(prob,i,repeat) = remake(prob; u0 = [LV_Vt0, LA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn], p=pchunk[:,i])

        ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)

        @time "    ODE Ensemble solved in " sol = solve(ensemble_prob, Tsit5(),  reltol=1e-4, abstol=1e-4, EnsembleThreads();saveat = saveat,trajectories=chunksize)
        @time "    Results processed in" Threads.@threads for i in 1:chunksize
          # println(i+offset)
          tsol = sol[i].t
          LV_V = sol[i][1,:]
          out[1,i+offset] = mean(LV_p(LV_V,tsol, pchunk[:,i]))
          out[2,i+offset] = maximum(sol[i][5,:])
          out[3,i+offset] = maximum(LV_V)
        end
    end
    global out_global = out
    out
end

## 100K sample test Sobol estimator ##

println("Sobol Estimator - Sobol")
#save("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Discrete/Discrete Data/Sobol_HS.jld", "data", Sobol_HS)
chunks::Int64 = 10000
samples = 100000
lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
@time Sobol_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
save("Shi_Sobol_S.jld", "data", Sobol_S)
#save("Shi_Sobol_S_ST_conf.jld", "data", b2)


###### Sobol Sampler #####
#
#println("Jansen Estimator - Sobol")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = Sobol_J.ST
#    b2[:,:,i] = Sobol_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Sobol_J_ST.jld", "data", b1)
#save("Shi_Sobol_J_ST_conf.jld", "data", b2)
#
#println("Sobol Estimator - Sobol")
##save("/home/harry/Desktop/PhD/Year 2/GSA project/Nik BS and convergence/Discrete/Discrete Data/Sobol_HS.jld", "data", Sobol_HS)
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = Sobol_S.ST
#    b2[:,:,i] = Sobol_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Sobol_S_ST.jld", "data", b1)
#save("Shi_Sobol_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - Sobol")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = Sobol_H.ST
#    b2[:,:,i] = Sobol_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Sobol_H_ST.jld", "data", b1)
#save("ShiD_Sobol_H_ST_conf.jld", "data", b2)
#
#println("Janon Estimator - Sobol")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    sampler = SobolSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
#    @time Sobol_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = Sobol_JA.ST
#    b2[:,:,i] = Sobol_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Sobol_JA_ST.jld", "data", b1)
#save("ShiD_Sobol_JA_ST_conf.jld", "data", b2)

###### Lattice Rule Sampling ########
#
#println("Jansen Estimator - Lattice Rule")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = lattice_J.ST
#    b2[:,:,i] = lattice_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Lattice_J_ST.jld", "data", b1)
#save("Shi_Lattice_J_ST_conf.jld", "data", b2)
#
#println("Sobol Estimator - Lattice Rule")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = lattice_S.ST
#    b2[:,:,i] = lattice_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Lattice_S_ST.jld", "data", b1)
#save("Shi_Lattice_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - Lattice Rule")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = lattice_H.ST
#    b2[:,:,i] = lattice_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Lattice_H_ST.jld", "data", b1)
#save("ShiD_Lattice_H_ST_conf.jld", "data", b2)

#println("Janon Estimator - Lattice Rule")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatticeRuleSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time lattice_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = lattice_JA.ST
#    b2[:,:,i] = lattice_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Lattice_JA_ST.jld", "data", b1)
#save("ShiD_Lattice_JA_ST_conf.jld", "data", b2)

#### Golden Sampling ###
#
#println("Jansen Estimator - Golden")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = Golden_J.ST
#    b2[:,:,i] = Golden_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Golden_J_ST.jld", "data", b1)
#save("Shi_Golden_J_ST_conf.jld", "data", b2)
#
#println("Sobol Estimator - Golden")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = Golden_S.ST
#    b2[:,:,i] = Golden_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Golden_S_ST.jld", "data", b1)
#save("Shi_Golden_S_ST_conf.jld", "data", b2)


#println("Homma Estimator - Golden")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = Golden_H.ST
#    b2[:,:,i] = Golden_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Golden_H_ST.jld", "data", b1)
#save("ShiD_Golden_H_ST_conf.jld", "data", b2)

#println("Janon Estimator - Golden")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = GoldenSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time Golden_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = Golden_JA.ST
#    b2[:,:,i] = Golden_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Golden_JA_ST.jld", "data", b1)
#save("ShiD_Golden_JA_ST_conf.jld", "data", b2)

## Uniform Sampling 

#println("Jansen Estimator - Uniform")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = Uniform_J.ST
#    b2[:,:,i] = Uniform_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Uniform_J_ST.jld", "data", b1)
#save("Shi_Uniform_J_ST_conf.jld", "data", b2)
#
#
#println("Sobol Estimator - Uniform")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = Uniform_S.ST
#    b2[:,:,i] = Uniform_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_Uniform_S_ST.jld", "data", b1)
#save("Shi_Uniform_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - Uniform")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = Uniform_H.ST
#    b2[:,:,i] = Uniform_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Uniform_H_ST.jld", "data", b1)
#save("ShiD_Uniform_H_ST_conf.jld", "data", b2)
#
#println("Janon Estimator - Uniform")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    A = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    B = QuasiMonteCarlo.sample(samples, lb, ub, Uniform())
#    @time Uniform_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#    b1[:,:,i] = Uniform_JA.ST
#    b2[:,:,i] = Uniform_JA.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_Uniform_JA_ST.jld", "data", b1)
#save("ShiD_Uniform_JA_ST_conf.jld", "data", b2)

### Latin HyperCube ###

#println("Jansen Estimator - LH")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_J = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)
#    b1[:,:,i] = LH_J.ST
#    b2[:,:,i] = LH_J.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_LH_J_ST.jld", "data", b1)
#save("Shi_LH_J_ST_conf.jld", "data", b2)
#
#
#println("Sobol Estimator - LH")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_S = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Sobol2007)
#    b1[:,:,i] = LH_S.ST
#    b2[:,:,i] = LH_S.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("Shi_LH_S_ST.jld", "data", b1)
#save("Shi_LH_S_ST_conf.jld", "data", b2)

#println("Homma Estimator - LH")
#b1 = zeros(3,31,5)
#b2 = zeros(3,31,5)
#N = range(start = 10000, stop = 30000, step = 5000) #± 10%
#chunks::Int64 = 1000
#for i in 1:length(N)
#    samples = N[i]
#    lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#    ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#    bounds = tuple.(lb,ub)
#    sampler = LatinHypercubeSample()
#    A,B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
#    @time LH_H = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Homma1996)
#    b1[:,:,i] = LH_H.ST
#    b2[:,:,i] = LH_H.ST_Conf_Int
#    println("Iteration at ", i)
#end 
#save("ShiD_LH_H_ST.jld", "data", b1)
#save("ShiD_LH_H_ST_conf.jld", "data", b2)

#println("Janon Estimator - LH")
#chunks::Int64 = 1000
#
#samples = 10000
#lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#A = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#B = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#@time LH_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#save("ShiD_LH_JA_10k.jld", "data", LH_JA)
#
#samples = 15000
#lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#A = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#B = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#@time LH_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#save("ShiD_LH_JA_15k.jld", "data", LH_JA)
#
#samples = 20000
#chunks::Int64 = 2000
#lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#A = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#B = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#@time LH_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#save("ShiD_LH_JA_20k.jld", "data", LH_JA)
#
#samples = 25000
#chunks::Int64 = 5000
#lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#A = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#B = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#@time LH_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#save("ShiD_LH_JA_25k.jld", "data", LH_JA)
#
#samples = 30000
#lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
#ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
#A = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#B = QuasiMonteCarlo.sample(samples, lb, ub, LatinHypercubeSample())
#@time LH_JA = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Janon2014)
#save("ShiD_LH_JA_30k.jld", "data", LH_JA)

