using CairoMakie, HDF5, Statistics, LinearAlgebra, JLD
CairoMakie.activate!(type = "svg")
# problem defintion
begin 
    using DifferentialEquations
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

end

t = Array(LinRange(15, 16, 300))

ST = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Discrete/Sobol_JaU_ST.jld")["data"]


## Model analysis
begin
    # Parameter importance PCA 

    F = transpose(ST)*ST

    e_deomp=eigen(F)

    λ = abs.(e_deomp.values)
    Q = abs.(e_deomp.vectors)

    e_value_sum = sum(λ)

    e = Vector{Float64}(undef,9)
    for i in 1:9
        for j in 1:9
        e[i] = sum(λ[j]*Q[i,j])/e_value_sum
        end 
    end 
    p=sortperm(e,rev=true)

    # orthogonality calculation

    Orth_heat1 = Matrix{Float64}(undef,length(prob.p),length(prob.p))
    for j in 1:length(prob.p)
        for i in 1:length(prob.p)
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
    a = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))
end


## plot the final figure 
begin
    f = Figure(resolution = (900, 600),backgroundcolor = RGBf(0.98, 0.98, 0.98));

    #Heatmap of measurements for Sobol analysis
    ax = Axis(f[1,1], xticklabelrotation = π / 3, xticklabelalign = (:right, :center), xticks = (1:3, [L"LV.P", L"SA.P", L"LV.V"]), yticks = (1:9, [L"τ_{es}", L"τ_{ep}", L"Rmv", L"Zao", L"Rs", L"Csa", L"Csv", L"E_{max}", L"E_{min}"]), title = L"Sobol - Total~Order", xlabel = L"Measurements", ylabel = L"Parameters")
    hm = CairoMakie.heatmap!(ax,ST, colormap=:plasma)
    for i in 1:3, j in 1:9
        txtcolor = ST[i, j] < -1000.0 ? :white : :black
        text!(ax, "$(round(ST[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    CairoMakie.Colorbar(f[1,2],hm);

    # PCA importance measure 
    ax = Axis(f[1,3],title=L"Parameter~Importance", xticks = (1:9, [L"τ_{es}", L"τ_{ep}", L"Rmv", L"Zao", L"Rs", L"Csa", L"Csv", L"E_{max}", L"E_{min}"][p]), xlabel = L"Parameters", ylabel = L"Importance")
    CairoMakie.scatter!( e[p], label=L"Condition~Number")

    # orthogonality heatmap 
    ax1 = Axis(f[2,1], xticks = (1:9, [L"τ_{es}", L"τ_{ep}", L"Rmv", L"Zao", L"Rs", L"Csa", L"Csv", L"E_{max}", L"E_{min}"]), yticks = (1:9, [L"τ_{es}", L"τ_{ep}", L"Rmv", L"Zao", L"Rs", L"Csa", L"Csv", L"E_{max}", L"E_{min}"]), title = L"ST-Orthogonality~Matrix")
    hm1 = CairoMakie.heatmap!(ax1,data, colormap=:plasma)
    for i in 1:9, j in 1:9
        txtcolor = data[i, j] < -0.0 ? :white : :black
        text!(ax1, "$(round(data[i,j], digits = 2))", position = (i, j),
            color = txtcolor, align = (:center, :center), fontsize = 15)
    end
    CairoMakie.Colorbar(f[2,2],hm1, label = L"Orthogonality~Score", ticks = 0.0:0.2:1.0)

    # Histogram of orthogonality 
    ax = Axis(f[2,3], xticks = 0.0:0.1:1.0, title = L"Orthogonality~Spread", xlabel = L"Orthogonality~Score", ylabel = L"Density")
    hist!(ax, a, color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1,1], f[1,3], f[2,1], f[2,3]])
        Label(layout[1, 1, TopRight()], label,fontsize = 18,font = :bold,halign = :right)
    end

    f
end 


#save("/home/harry/Desktop/PhD/Year 2/GSA project/Plots/Janon/Discrete/Latin_Hypercube.svg",f)

#### plots for orthogonality ranking ### 
begin
    using LaTeXStrings, Plots
    data

    a = zeros(9)

    for i in 1:9
    a[i] = sum(data[:,i])/length(a)
    end

    p=sortperm(a,rev=true)

    param_label = [L"\tau_{es}", L"\tau_{ep}", L"Rmv", L"Zao", L"Rs", L"Csa", L"Csv",L"E_{max}", L"E_{min}"]
    Plots.scatter(a[p], title=L"Orth~Ranking",legend=false, xticks = (1:9, param_label[p]))
end