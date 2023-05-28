using Dates
println("Running on ", Threads.nthreads(), " threads.")
println(now())

@time "Setup modules" using ModelingToolkit, DifferentialEquations, Statistics, CirculatorySystemModels,  Distributions, QuasiMonteCarlo, Serialization
include("GlobalSensitivity.jl")

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
``
## And simplify it
@time "Simplify model"  circ_sys = structural_simplify(circ_model)

u0 =  [LV_Vt0, LA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn]
@time "Define problem" prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))

saveat = LinRange(15, 16, 300)
@time "Solve first time" sol = solve(prob, Tsit5(),  reltol=1e-8, abstol=1e-8, saveat = saveat)
@time "Solve second time" sol = solve(prob, Tsit5(),  reltol=1e-8, abstol=1e-8, saveat = saveat)



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

        @time "    ODE Ensemble solved in " sol = solve(ensemble_prob,Vern7(), reltol=1e-8, abstol=1e-8, EnsembleThreads();saveat = saveat,trajectories=chunksize)
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


chunks::Int64 = 100

##### Sobol Sampler #####
N = 19800
# 30%
lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
bounds = tuple.(lb,ub)
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)

println("Homma Estimator - Sobol")
@time Sobol_HS = gsa(circ_time_idxs_chunked,Sobol(),A,B,batch = true, Ei_estimator = :Homma1996)
serialize("Sobol_HS.jls", Sobol_HS)


println("Sobol Estimator - Sobol")
@time Sobol_SS = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Sobol2007)
serialize("Sobol_SS.jls", Sobol_SS)


println("Jansen Estimator - Sobol")
@time Sobol_JS = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
serialize("Sobol_JS.jls", Sobol_JS)

println("Janon Estimator - Sobol")
@time Sobol_JaS = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Janon2014)
serialize("Sobol_JaS.jls", Sobol_JaS)

##### Lattice Rule Sampling ########
N = 19800
# 30%
lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
bounds = tuple.(lb,ub)
sampler = LatticeRuleSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)

println("Homma Estimator - Lattice Rule")
@time Sobol_HLR = gsa(circ_time_idxs_chunked,Sobol(),A,B,batch = true, Ei_estimator = :Homma1996)
serialize("Sobol_HLR.jls", Sobol_HLR)

println("Sobol Estimator - Lattice Rule")
@time Sobol_SLR = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Sobol2007)
serialize("Sobol_SLR.jls", Sobol_SLR)

println("Jansen Estimator - Lattice Rule")
@time Sobol_JLR = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
serialize("Sobol_JLR.jls", Sobol_JLR)

println("Janon Estimator - Lattice Rule")
@time Sobol_JaLR = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Janon2014)
serialize("Sobol_JaLR.jls", Sobol_JaLR)
### Golden Sampling ###
N = 19800
# 30%
lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
bounds = tuple.(lb,ub)
sampler = GoldenSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N, lb, ub, sampler)

println("Homma Estimator - Golden")
@time Sobol_HG = gsa(circ_time_idxs_chunked,Sobol(),A,B,batch = true, Ei_estimator = :Homma1996)
serialize("Sobol_HG.jls", Sobol_HG)

println("Sobol Estimator - Golden")
@time Sobol_SG = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Sobol2007)
serialize("Sobol_SG.jls", Sobol_SG)

println("Jansen Estimator - Golden")
@time Sobol_JG = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
serialize("Sobol_JG.jls", Sobol_JG)

println("Janon Estimator - Golden")
@time Sobol_JaG = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Janon2014)
serialize("Sobol_JaG.jls", Sobol_JaG)

## Uniform Sampling 
N = 19800
# 30%
lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
bounds = tuple.(lb,ub)
sampler = UniformSample()
A = deserialize("Uniform_Sampler_A.jls")
B = deserialize("Uniform_Sampler_B.jls")

println("Homma Estimator - Uniform")
@time Sobol_HU = gsa(circ_time_idxs_chunked,Sobol(),A,B,batch = true, Ei_estimator = :Homma1996)
serialize("Sobol_HU.jls", Sobol_HU)

println("Sobol Estimator - Uniform")
@time Sobol_SU = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Sobol2007)
serialize("Sobol_SU.jls", Sobol_SU)

println("Jansen Estimator - Uniform")
@time Sobol_JU = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
serialize("Sobol_JU.jls", Sobol_JU)

println("Janon Estimator - Uniform")
@time Sobol_JaU = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Janon2014)
serialize("Sobol_JaU.jls", Sobol_JaU)

### Latin HyperCube ###
N = 19800
# 30%
lb = [5.0, 1.0, 0.07, 1.75, 1.0, 0.21, 0.360, 0.0, 4.0, 1.0, 0.105, 0.175, 1.0, 0.0315, 0.063, 0.92, 0.0231, 0.042, 0.0, 0.056, 0.0021, 4.34e-5, 0.0, 1.12, 0.035, 0.00119, 0.35, 0.364, 0.0525, 0.0, 14.35]
ub = [5.0, 1.0, 0.13, 3.25, 1.0, 0.34, 0.585, 0.0, 4.0, 1.0, 0.165, 0.325, 1.0, 0.0585, 0.117, 0.92, 0.0429, 0.078, 0.0, 0.104, 0.0039, 8.06e-5, 0.0, 2.08, 0.065, 0.00221, 0.65, 0.676, 0.0975, 0.0, 26.65]
bounds = tuple.(lb,ub)
sampler = LatinHypercubeSample()
A = deserialize("Latin_Hypercube_Sampler_A.jls")
B = deserialize("Latin_Hypercube_Sampler_B.jls")

println("Homma Estimator - Uniform")
@time Sobol_HLH = gsa(circ_time_idxs_chunked,Sobol(),A,B,batch = true, Ei_estimator = :Homma1996)
serialize("Sobol_HLH.jls", Sobol_HLH)

println("Sobol Estimator - Uniform")
@time Sobol_SLH = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Sobol2007)
serialize("Sobol_SLH.jls", Sobol_SLH)

println("Jansen Estimator - Uniform")
@time Sobol_JLH = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Jansen1999)
serialize("Sobol_JLH.jls", Sobol_JLH)

println("Janon Estimator - Uniform")
@time Sobol_JaLH = gsa(circ_time_post,Sobol(),A,B,batch = true, Ei_estimator = :Janon2014)
serialize("Sobol_JaLH.jls", Sobol_JaLH)

Sobol_JaLH = deserialize("/home/harry/Downloads/Discrete_Measurements/Sobol_JaLH.jls")
save("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/DATA/Janon/Sobol_JaLH.jld", "data", Sobol_JaLH.ST)