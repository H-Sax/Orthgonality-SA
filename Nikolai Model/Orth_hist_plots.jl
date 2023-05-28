using CairoMakie, HDF5, Statistics, LinearAlgebra, JLD
CairoMakie.activate!(type = "svg")

diff = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JS_ST.jld")["data"]
### Discrete
begin
    ## Homma ##
    ST_HS = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Discrete/Sobol_HS_ST.jld")["data"]
    ST_HLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Discrete/Sobol_HLR_ST.jld")["data"]
    ST_HG = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Discrete/Sobol_HG_ST.jld")["data"]
    ST_HU = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Discrete/Sobol_HU_ST.jld")["data"]
    ST_HLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Discrete/Sobol_HLH_ST.jld")["data"]

    ## Sobol ##
    ST_SS = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Discrete/Sobol_SS_ST.jld")["data"]
    ST_SLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Discrete/Sobol_SLR_ST.jld")["data"]
    ST_SG = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Discrete/Sobol_SG_ST.jld")["data"]
    ST_SU = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Discrete/Sobol_SU_ST.jld")["data"]
    ST_SLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Discrete/Sobol_SLH_ST.jld")["data"]

    ## Jansen ##
    ST_JS = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JS_ST.jld")["data"]
    ST_JLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JLR_ST.jld")["data"]
    ST_JG = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JG_ST.jld")["data"]
    ST_JU = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JU_ST.jld")["data"]
    ST_JLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Discrete/Sobol_JLH_ST.jld")["data"]

    ## Janon ##
    ST_JaS = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Discrete/Sobol_JaS_ST.jld")["data"]
    ST_JaLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Discrete/Sobol_JaLR_ST.jld")["data"]
    ST_JaG = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Discrete/Sobol_JaG_ST.jld")["data"]
    ST_JaU = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Discrete/Sobol_JaU_ST.jld")["data"]
    ST_JaLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Discrete/Sobol_JaLH_ST.jld")["data"]
end


## Continous 

begin
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
    
    ## Homma ##
    ST_HS =  h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Sobol_HS_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_HLR =  h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Sobol_HLR_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_HG =  h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Sobol_HG_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_HU =  h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Sobol_HU_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_HLH = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Homma/Sobol_HLH_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end

    ## Sobol ##
    ST_SS = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Sobol_SS_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_SLR = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Sobol_SLR_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_SG = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Sobol_SG_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_SU = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Sobol_SU_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_SLH = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/SOBOL/Sobol_SLH_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end

    ## Jansen ##
    ST_JS = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Sobol_JS_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_JLR = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Sobol_JLR_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_JG = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Sobol_JG_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_JU = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Sobol_JU_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end
    ST_JLH = h5open("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/JANSEN/Sobol_JLH_ST.h5", "r") do file
        read(file, "A")  # alternatively, say "@write file A"
    end

    ## Janon ##
    ST_JaS = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Sobol_JaS_ST.jld")["data"]
    ST_JaLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Sobol_JaLR_ST.jld")["data"]
    ST_JaG = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Sobol_JaG_ST.jld")["data"]
    ST_JaU = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Sobol_JaU_ST.jld")["data"]
    ST_JaLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/DATA/Janon/Sobol_JaLH_ST.jld")["data"]
    t = Array(LinRange(15, 16, 300))

    ## Breaking up into parts

    # Hommma
    ST_LVp_HS = ST_HS[1:300,:]
    ST_SAp_HS = ST_HS[301:600,:]
    ST_LVv_HS = ST_HS[601:900,:]

    ST_LVp_HLR = ST_HLR[1:300,:]
    ST_SAp_HLR = ST_HLR[301:600,:]
    ST_LVv_HLR = ST_HLR[601:900,:]

    ST_LVp_HG = ST_HG[1:300,:]
    ST_SAp_HG = ST_HG[301:600,:]
    ST_LVv_HG = ST_HG[601:900,:]

    ST_LVp_HU = ST_HU[1:300,:]
    ST_SAp_HU = ST_HU[301:600,:]
    ST_LVv_HU = ST_HU[601:900,:]

    ST_LVp_HLH = ST_HLH[1:300,:]
    ST_SAp_HLH = ST_HLH[301:600,:]
    ST_LVv_HLH = ST_HLH[601:900,:]

    # Sobol

    ST_LVp_SS = ST_SS[1:300,:]
    ST_SAp_SS = ST_SS[301:600,:]
    ST_LVv_SS = ST_SS[601:900,:]

    ST_LVp_SLR = ST_SLR[1:300,:]
    ST_SAp_SLR = ST_SLR[301:600,:]
    ST_LVv_SLR = ST_SLR[601:900,:]

    ST_LVp_SG = ST_SG[1:300,:]
    ST_SAp_SG = ST_SG[301:600,:]
    ST_LVv_SG = ST_SG[601:900,:]

    ST_LVp_SU = ST_SU[1:300,:]
    ST_SAp_SU = ST_SU[301:600,:]
    ST_LVv_SU = ST_SU[601:900,:]

    ST_LVp_SLH = ST_SLH[1:300,:]
    ST_SAp_SLH = ST_SLH[301:600,:]
    ST_LVv_SLH = ST_SLH[601:900,:]


    ## jansen 
    ST_LVp_JS = ST_JS[1:300,:]
    ST_SAp_JS = ST_JS[301:600,:]
    ST_LVv_JS = ST_JS[601:900,:]

    ST_LVp_JLR = ST_JLR[1:300,:]
    ST_SAp_JLR = ST_JLR[301:600,:]
    ST_LVv_JLR = ST_JLR[601:900,:]

    ST_LVp_JG = ST_JG[1:300,:]
    ST_SAp_JG = ST_JG[301:600,:]
    ST_LVv_JG = ST_JG[601:900,:]

    ST_LVp_JU = ST_JU[1:300,:]
    ST_SAp_JU = ST_JU[301:600,:]
    ST_LVv_JU = ST_JU[601:900,:]

    ST_LVp_JLH = ST_JLH[1:300,:]
    ST_SAp_JLH = ST_JLH[301:600,:]
    ST_LVv_JLH = ST_JLH[601:900,:]

    ## Janon

    ST_LVp_JaS = ST_JaS[1:300,:]
    ST_SAp_JaS = ST_JaS[301:600,:]
    ST_LVv_JaS = ST_JaS[601:900,:]

    ST_LVp_JaLR = ST_JaLR[1:300,:]
    ST_SAp_JaLR = ST_JaLR[301:600,:]
    ST_LVv_JaLR = ST_JaLR[601:900,:]

    ST_LVp_JaG = ST_JaG[1:300,:]
    ST_SAp_JaG = ST_JaG[301:600,:]
    ST_LVv_JaG = ST_JaG[601:900,:]

    ST_LVp_JaU = ST_JaU[1:300,:]
    ST_SAp_JaU = ST_JaU[301:600,:]
    ST_LVv_JaU = ST_JaU[601:900,:]

    ST_LVp_JaLH = ST_JLH[1:300,:]
    ST_SAp_JaLH = ST_JLH[301:600,:]
    ST_LVv_JaLH = ST_JLH[601:900,:]

    function dis_LVP(S)
        S_LVP = Vector{Float64}(undef,length(prob.p))
        for i in 1:length(prob.p)
            for j in 1:length(t)
            S_LVP[i] = (sum(S[1:j,i]*var(sol[1,1:j])))/(sum(var(sol[1,1:j])))/length(t)
            end 
        end 
        return S_LVP 
    end 

    function dis_SAP(S)
        S_SAP = Vector{Float64}(undef,length(prob.p))
        for i in 1:length(prob.p)
            for j in 1:length(t)
            S_SAP[i] = (sum(S[1:j,i]*var(sol[2,1:j])))/(sum(var(sol[2,1:j])))/length(t)
            end 
        end 
        return S_SAP
    end 

    function dis_LVV(S)
        S_LVV = Vector{Float64}(undef,length(prob.p))
        for i in 1:length(prob.p)
            for j in 1:length(t)
            S_LVV[i] = (sum(S[1:j,i]*var(sol[4,1:j])))/(sum(var(sol[4,1:j])))/length(t)
            end 
        end 
        return S_LVV
    end 

    ## Homma Sobol 

    SHS_LVP = dis_LVP(ST_LVp_HS)
    SHS_SAP = dis_SAP(ST_SAp_HS)
    SHS_LVV = dis_LVV(ST_LVv_HS)

    show(SHS_LVV)
    ST_HS = [SHS_LVP SHS_SAP SHS_LVV]'

    ## Homma Lattice Rule 

    SHLR_LVP = dis_LVP(ST_LVp_HLR)
    SHLR_SAP = dis_SAP(ST_SAp_HLR)
    SHLR_LVV = dis_LVV(ST_LVv_HLR)

    ST_HLR = [SHLR_LVP SHLR_SAP SHLR_LVV]'

    ## Homma Golden

    SHG_LVP = dis_LVP(ST_LVp_HG)
    SHG_SAP = dis_SAP(ST_SAp_HG)
    SHG_LVV = dis_LVV(ST_LVv_HG)

    ST_HG = [SHG_LVP SHG_SAP SHG_LVV]'
    
    ## Homma Golden

    SHU_LVP = dis_LVP(ST_LVp_HU)
    SHU_SAP = dis_SAP(ST_SAp_HU)
    SHU_LVV = dis_LVV(ST_LVv_HU)
    
    ST_HU = [SHU_LVP SHU_SAP SHU_LVV]'
    
    ## Homma Uniform
    SHU_LVP = dis_LVP(ST_LVp_HU)
    SHU_SAP = dis_SAP(ST_SAp_HU)
    SHU_LVV = dis_LVV(ST_LVv_HU)
    ST_HU = [SHU_LVP SHU_SAP SHU_LVV]'

    ## Homma Latin Hypercube
    SHLH_LVP = dis_LVP(ST_LVp_HLH)
    SHLH_SAP = dis_SAP(ST_SAp_HLH)
    SHLH_LVV = dis_LVV(ST_LVv_HLH)
    ST_HLH = [SHLH_LVP SHLH_SAP SHLH_LVV]'


    ## Sobol Sobol 

    SSS_LVP = dis_LVP(ST_LVp_SS)
    SSS_SAP = dis_SAP(ST_SAp_SS)
    SSS_LVV = dis_LVV(ST_LVv_SS)

    ST_SS = [SSS_LVP SSS_SAP SSS_LVV]'

    ## Sobol Lattice Rule 

    SSLR_LVP = dis_LVP(ST_LVp_SLR)
    SSLR_SAP = dis_SAP(ST_SAp_SLR)
    SSLR_LVV = dis_LVV(ST_LVv_SLR)

    ST_SLR = [SSLR_LVP SSLR_SAP SSLR_LVV]'

    ## Sobol Golden

    SSG_LVP = dis_LVP(ST_LVp_SG)
    SSG_SAP = dis_SAP(ST_SAp_SG)
    SSG_LVV = dis_LVV(ST_LVv_SG)

    ST_SG = [SSG_LVP SSG_SAP SSG_LVV]'
    
    ## Sobol Golden

    SSU_LVP = dis_LVP(ST_LVp_SU)
    SSU_SAP = dis_SAP(ST_SAp_SU)
    SSU_LVV = dis_LVV(ST_LVv_SU)
    
    ST_SU = [SSU_LVP SSU_SAP SSU_LVV]'
    
    ## Sobol Uniform
    SSU_LVP = dis_LVP(ST_LVp_SU)
    SSU_SAP = dis_SAP(ST_SAp_SU)
    SSU_LVV = dis_LVV(ST_LVv_SU)
    ST_SU = [SSU_LVP SSU_SAP SSU_LVV]'

    ## Sobol Latin Hypercube
    SSLH_LVP = dis_LVP(ST_LVp_SLH)
    SSLH_SAP = dis_SAP(ST_SAp_SLH)
    SSLH_LVV = dis_LVV(ST_LVv_SLH)
    ST_SLH = [SSLH_LVP SSLH_SAP SSLH_LVV]'

    ## Jansen Sobol 

    SJS_LVP = dis_LVP(ST_LVp_JS)
    SJS_SAP = dis_SAP(ST_SAp_JS)
    SJS_LVV = dis_LVV(ST_LVv_JS)

    ST_JS = [SJS_LVP SJS_SAP SJS_LVV]'

    ## Jansen Lattice Rule 

    SJLR_LVP = dis_LVP(ST_LVp_JLR)
    SJLR_SAP = dis_SAP(ST_SAp_JLR)
    SJLR_LVV = dis_LVV(ST_LVv_JLR)

    ST_JLR = [SJLR_LVP SJLR_SAP SJLR_LVV]'

    ## Jansen Golden

    SJG_LVP = dis_LVP(ST_LVp_JG)
    SJG_SAP = dis_SAP(ST_SAp_JG)
    SJG_LVV = dis_LVV(ST_LVv_JG)

    ST_JG = [SJG_LVP SJG_SAP SJG_LVV]'
    
    ## Jansen Golden

    SJU_LVP = dis_LVP(ST_LVp_JU)
    SJU_SAP = dis_SAP(ST_SAp_JU)
    SJU_LVV = dis_LVV(ST_LVv_JU)
    
    ST_JU = [SJU_LVP SJU_SAP SJU_LVV]'
    
    ## Jansen Uniform
    SJU_LVP = dis_LVP(ST_LVp_JU)
    SJU_SAP = dis_SAP(ST_SAp_JU)
    SJU_LVV = dis_LVV(ST_LVv_JU)
    ST_JU = [SJU_LVP SJU_SAP SJU_LVV]'

    ## Jansen Latin Hypercube
    SJLH_LVP = dis_LVP(ST_LVp_JLH)
    SJLH_SAP = dis_SAP(ST_SAp_JLH)
    SJLH_LVV = dis_LVV(ST_LVv_JLH)
    ST_JLH = [SJLH_LVP SJLH_SAP SJLH_LVV]'

    ## Janon Sobol 

    SJaS_LVP = dis_LVP(ST_LVp_JaS)
    SJaS_SAP = dis_SAP(ST_SAp_JaS)
    SJaS_LVV = dis_LVV(ST_LVv_JaS)

    ST_JaS = [SJaS_LVP SJaS_SAP SJaS_LVV]'

    ## Janon Lattice Rule 

    SJaLR_LVP = dis_LVP(ST_LVp_JaLR)
    SJaLR_SAP = dis_SAP(ST_SAp_JaLR)
    SJaLR_LVV = dis_LVV(ST_LVv_JaLR)

    ST_JaLR = [SJaLR_LVP SJaLR_SAP SJaLR_LVV]'

    ## Janon Golden

    SJaG_LVP = dis_LVP(ST_LVp_JaG)
    SJaG_SAP = dis_SAP(ST_SAp_JaG)
    SJaG_LVV = dis_LVV(ST_LVv_JaG)

    ST_JaG = [SJaG_LVP SJaG_SAP SJaG_LVV]'
    
    ## Janon Golden

    SJaU_LVP = dis_LVP(ST_LVp_JaU)
    SJaU_SAP = dis_SAP(ST_SAp_JaU)
    SJaU_LVV = dis_LVV(ST_LVv_JaU)
    
    ST_JaU = [SJaU_LVP SJaU_SAP SJaU_LVV]'
    
    ## Janon Uniform
    SJaU_LVP = dis_LVP(ST_LVp_JaU)
    SJaU_SAP = dis_SAP(ST_SAp_JaU)
    SJaU_LVV = dis_LVV(ST_LVv_JaU)
    ST_JaU = [SJaU_LVP SJaU_SAP SJaU_LVV]'

    ## Janon Latin Hypercube
    SJaLH_LVP = dis_LVP(ST_LVp_JaLH)
    SJaLH_SAP = dis_SAP(ST_SAp_JaLH)
    SJaLH_LVV = dis_LVV(ST_LVv_JaLH)
    ST_JaLH = [SJaLH_LVP SJaLH_SAP SJaLH_LVV]'
end



function hist_(S)

    Orth_heat1 = Matrix{Float64}(undef,9,9)
    for j in 1:9
        for i in 1:9
            if i==j
                Orth_heat1[i,j] = 0
            else 

            Orth_heat1[i,j] = sin(acos(((transpose(S[:,i])*S[:,j]))/(norm(S[:,i])*norm(S[:,j]))-1e-15))  #Slight numerical rounding error without the additional add on 
            end 
        end 
    end 
    data = Orth_heat1

    a1 = Orth_heat1[2:end,1]
    a2 = Orth_heat1[3:end,2]
    a3 = Orth_heat1[4:end,3]
    a4 = Orth_heat1[5:end,4]
    a5 = Orth_heat1[6:end,5]
    a6 = Orth_heat1[7:end,6]
    a7 = Orth_heat1[8:end,7]
    a8 = Orth_heat1[9:end,8]
    a = reduce(vcat, (a1,a2,a3,a4,a5,a6,a7,a8))

    return a
end 

function orth(S)

    Orth_heat1 = Matrix{Float64}(undef,9,9)
    for j in 1:9
        for i in 1:9
            if i==j
                Orth_heat1[i,j] = 0
            else 

            Orth_heat1[i,j] = sin(acos(((transpose(S[:,i])*S[:,j]))/(norm(S[:,i])*norm(S[:,j]))-1e-15))  #Slight numerical rounding error without the additional add on 
            end 
        end 
    end 
    data = Orth_heat1
    return data


end 


## Hist plot
begin
    f = Figure(resolution = (1400, 900),backgroundcolor = RGBf(0.98, 0.98, 0.98), fontsize = 20);
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

    ax = Axis(ga[1,1], xticks = 0.0:0.2:1.0)
    Label(ga[1, 1, Top()], "Homma", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    Label(ga[1, 1, Left()], "Sobol", valign = :center, font = :bold, rotation = pi/2, padding = (0, 40, 0, 0))
    hist!(ax, hist_(ST_HS), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gb[1,1], xticks = 0.0:0.2:1.0)  #ylabel = "Lattice Rule")
    Label(gb[1, 1, Left()], "Lattice Rule", valign = :center, font = :bold, rotation = pi/2, padding = (0, 40, 0, 0))
    hist!(ax, hist_(ST_HLR), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, font = :bold, normalization = :pdf)

    ax = Axis(gc[1,1], xticks = 0.0:0.2:1.0)
    Label(gc[1, 1, Left()], "Golden", valign = :center, font = :bold, rotation = pi/2, padding = (0, 40, 0, 0))
    hist!(ax, hist_(ST_HG), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gd[1,1], xticks = 0.0:0.2:1.0)
    Label(gd[1, 1, Left()], "Uniform", valign = :center, font = :bold, rotation = pi/2, padding = (0, 40, 0, 0))
    hist!(ax, hist_(ST_HU), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(ge[1,1], xticks = 0.0:0.2:1.0)
    Label(ge[1, 1, Left()], "Latin Hypercube", valign = :center, font = :bold, rotation = pi/2, padding = (0, 40, 0, 0))
    hist!(ax, hist_(ST_HLH), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gf[1,1], xticks = 0.0:0.2:1.0)
    Label(gf[1, 1, Top()], "Sobol", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    hist!(ax, hist_(ST_SS), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gg[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_SLR), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gh[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_SG), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gi[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_SU), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gj[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_SLH), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gk[1,1], xticks = 0.0:0.2:1.0)
    Label(gk[1, 1, Top()], "Jansen", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    hist!(ax, hist_(ST_JS), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gl[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JLR), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gm[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JG), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gn[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JU), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(go[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JLH), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gp[1,1], xticks = 0.0:0.2:1.0)
    Label(gp[1, 1, Top()], "Janon", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    hist!(ax, hist_(ST_JaS), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gq[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JaLR), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gr[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JaG), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gs[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JaU), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gt[1,1], xticks = 0.0:0.2:1.0)
    hist!(ax, hist_(ST_JaLH), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    for (label, layout) in zip(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"], [ga, gb, gc, gd, ge, gf, gg, gh, gi, gj, gk, gl, gm, gn, go, gp, gq, gr, gs, gt])
        Label(layout[1, 1, TopRight()], label,fontsize = 18,font = :bold,halign = :right)
    end
    
    resize_to_layout!(f)
    set_theme!(Theme(font = "Times New Roman bold"))

    f
end  


## Difference hist Plots
begin
    f = Figure(resolution = (1400, 900),backgroundcolor = RGBf(0.98, 0.98, 0.98), fontsize = 20);
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

    ax = Axis(ga[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(ga[1, 1, Top()], "Homma", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    Label(ga[1, 1, Left()], "Sobol", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, (hist_(ST_HS) - hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gb[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)  #ylabel = "Lattice Rule")
    Label(gb[1, 1, Left()], "Lattice Rule", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, (hist_(ST_HLR)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, font = :bold)

    ax = Axis(gc[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(gc[1, 1, Left()], "Golden", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, (hist_(ST_HG)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gd[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(gd[1, 1, Left()], "Uniform", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, (hist_(ST_HU)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(ge[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(ge[1, 1, Left()], "Latin Hypercube", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, (hist_(ST_HLH)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gf[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(gf[1, 1, Top()], "Sobol", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    hist!(ax, (hist_(ST_SS)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gg[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_SLR)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gh[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_SG)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gi[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_SU)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gj[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_SLH)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gk[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(gk[1, 1, Top()], "Jansen", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    hist!(ax, (hist_(ST_JS)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gl[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JLR)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gm[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JG)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gn[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JU)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(go[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JLH)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gp[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    Label(gp[1, 1, Top()], "Janon", valign = :center, font = :bold, padding = (0, 0, 10, 0))
    hist!(ax, (hist_(ST_JaS)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gq[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JaLR)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gr[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JaG)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gs[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JaU)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    ax = Axis(gt[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:2.0:10)
    hist!(ax, (hist_(ST_JaLH)- hist_(diff)), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black)

    for (label, layout) in zip(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"], [ga, gb, gc, gd, ge, gf, gg, gh, gi, gj, gk, gl, gm, gn, go, gp, gq, gr, gs, gt])
        Label(layout[1, 1, TopRight()], label,fontsize = 18,font = :bold,halign = :right)
    end
    
    resize_to_layout!(f)
    set_theme!(Theme(font = "Times New Roman bold"))

    f
end  
