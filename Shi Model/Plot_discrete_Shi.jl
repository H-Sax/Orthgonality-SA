using CairoMakie, HDF5, Statistics, LinearAlgebra, JLD
CairoMakie.activate!(type = "svg")
using ModelingToolkit, DifferentialEquations, Statistics, CirculatorySystemModels,  Distributions, QuasiMonteCarlo
# problem defintion
begin 
    
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

end

t = Array(LinRange(15, 16, 300))

ST = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon/Sobol_JaLH.jld")["data"]


## Model analysis
begin
    # Parameter importance PCA 

    F = transpose(ST)*ST

    e_deomp=eigen(F)

    λ = abs.(e_deomp.values)
    Q = abs.(e_deomp.vectors)

    e_value_sum = sum(λ)

    e = Vector{Float64}(undef,20)
    for i in 1:20
        for j in 1:20
        e[i] = sum(λ[j]*Q[i,j])/e_value_sum
        end 
    end 
    p=sortperm(e,rev=true)

    # orthogonality calculation

    Orth_heat1 = Matrix{Float64}(undef,20,20)
    for j in 1:20
        for i in 1:20
            if i==j
                Orth_heat1[i,j] = 0
            else 

            Orth_heat1[i,j] = sin(acos(((transpose(ST[:,i])*ST[:,j]))/(norm(ST[:,i])*norm(ST[:,j]))-1e-15))  #Slight numerical rounding error without the additional add on 
            end 
        end 
    end 
    data = Orth_heat1

    # histogram Plot 
    a1 = Orth_heat1[2:end,1]
    a2 = Orth_heat1[3:end,2]
    a3 = Orth_heat1[4:end,3]
    a4 = Orth_heat1[5:end,4]
    a5 = Orth_heat1[6:end,5]
    a6 = Orth_heat1[7:end,6]
    a7 = Orth_heat1[8:end,7]
    a8 = Orth_heat1[9:end,8]
    a9 = Orth_heat1[10:end,9]
    a10 = Orth_heat1[11:end,10]
    a11 = Orth_heat1[12:end,11]
    a12 = Orth_heat1[13:end,12]
    a13 = Orth_heat1[14:end,13]
    a14 = Orth_heat1[15:end,14]
    a15 = Orth_heat1[16:end,15]
    a16 = Orth_heat1[17:end,16]
    a17 = Orth_heat1[18:end,17]
    a18 = Orth_heat1[19:end,18]
    a19 = Orth_heat1[20:end,19]
    a = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19))
end


## plot the final figure 
begin
    f = Figure(resolution = (900, 600),backgroundcolor = RGBf(0.98, 0.98, 0.98));

    #Heatmap of measurements for Sobol analysis
    ax = Axis(f[1,1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:3, [L"LV.P", L"SA.P", L"LV.V"]), yticks = 1:20, title = L"Sobol - Total~Order", xlabel = L"Measurements", ylabel = L"Parameters")
    hm = CairoMakie.heatmap!(ax,ST, colormap=:plasma)
    for i in 1:3, j in 1:20
        txtcolor = ST[i, j] < -1000.0 ? :white : :black
        text!(ax, "$(round(ST[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    CairoMakie.Colorbar(f[1,2],hm);

    # PCA importance measure 
    ax = Axis(f[1,3],title=L"Parameter~Importance", xticks = 1:20, xlabel = L"Parameters", ylabel = L"Importance")
    CairoMakie.scatter!( e[p], label=L"Condition~Number")

    # orthogonality heatmap 
    ax1 = Axis(f[2,1], xticks = 1:20, yticks = 1:20, title = L"ST-Orthogonality~Matrix")
    hm1 = CairoMakie.heatmap!(ax1,data, colormap=:plasma)
    for i in 1:20, j in 1:20
        txtcolor = data[i, j] < -0.0 ? :white : :black
        text!(ax1, "$(round(data[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    CairoMakie.Colorbar(f[2,2],hm1, label = L"Orthogonality~Score", ticks = 0.0:0.2:1.0)

    # Histogram of orthogonality 
    ax = Axis(f[2,3], xticks = 0.0:0.1:1.0, title = L"Orthogonality~Spread", xlabel = L"Orthogonality~Score", ylabel = L"Density")
    hist!(ax, a, color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1,1], f[1,3], f[2,1], f[2,3]])
        Label(layout[1, 1, TopRight()], label,fontsize = 18, font = :bold,halign = :right)
    end

    f
end 


#save("/home/harry/Desktop/PhD/Year 2/GSA project/Plots/Janon/Discrete/Latin_Hypercube.svg",f)

#### plots for orthogonality ranking ### 
begin
    using LaTeXStrings, Plots
    data

    a = zeros(20)

    for i in 1:20
    a[i] = sum(data[:,i])/length(a)
    end

    p=sortperm(a,rev=true)

    param_label = [L"Emin_{lv}", L"Emax_{lv}", L"\tau es_{LV}", L"\tau ep_{LV}",L"Emin_{la}", L"Emax_{la}", L"\tau es_{La}", L"\tau ep_{La}", L"Zao", L"Rmv", L"C_{sas}", L"R_{sas}", L"L_{sas}", L"C_{sat}", L"R_{sat}", L"L_{sat}", L"R_{sar}", L"R_{scp}", L"R_{svn}", L"C_{svn}"]
    Plots.scatter(a[p], title=L"Orth~Ranking",legend=false, xticks = (1:20, param_label[p]))
    Plots.plot!(size = (900, 400))
end