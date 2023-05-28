using CairoMakie, HDF5, Statistics, LinearAlgebra, JLD
CairoMakie.activate!(type = "svg")
### Discrete
begin
    ## Homma ##
    ST_HS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma/Sobol_HS.jld")["data"]
    ST_HLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma/Sobol_HLR.jld")["data"]
    ST_HG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma/Sobol_HG.jld")["data"]
    ST_HU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma/Sobol_HU.jld")["data"]
    ST_HLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma/Sobol_HLH.jld")["data"]

    ## Sobol ##
    ST_SS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol/Sobol_SS.jld")["data"]
    ST_SLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol/Sobol_SLR.jld")["data"]
    ST_SG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol/Sobol_SG.jld")["data"]
    ST_SU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol/Sobol_SU.jld")["data"]
    ST_SLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol/Sobol_SLH.jld")["data"]

    ## Jansen ##
    ST_JS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen/Sobol_JS.jld")["data"]
    ST_JLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen/Sobol_JLR.jld")["data"]
    ST_JG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen/Sobol_JG.jld")["data"]
    ST_JU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen/Sobol_JU.jld")["data"]
    ST_JLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen/Sobol_JLH.jld")["data"]

    ## Janon ##
    ST_JaS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon/Sobol_JaS.jld")["data"]
    ST_JaLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon/Sobol_JaLR.jld")["data"]
    ST_JaG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon/Sobol_JaG.jld")["data"]
    ST_JaU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon/Sobol_JaU.jld")["data"]
    ST_JaLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon/Sobol_JaLH.jld")["data"]
end


## Continous 

begin
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
    
    ## Homma ##
    ST_HS =  load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma_C/Sobol_HS.jld")["data"]
    ST_HLR =  load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma_C/Sobol_HLR.jld")["data"]
    ST_HG =  load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma_C/Sobol_HG.jld")["data"]
    ST_HU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma_C/Sobol_HU.jld")["data"]
    ST_HLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Homma_C/Sobol_HLH.jld")["data"]

    ## Sobol ##
    ST_SS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol_C/Sobol_SS.jld")["data"]
    ST_SLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol_C/Sobol_SLR.jld")["data"]
    ST_SG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol_C/Sobol_SG.jld")["data"]
    ST_SU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol_C/Sobol_SU.jld")["data"]
    ST_SLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Sobol_C/Sobol_SLH.jld")["data"]

    ## Jansen ##
    ST_JS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen_C/Sobol_JS.jld")["data"]
    ST_JLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen_C/Sobol_JLR.jld")["data"]
    ST_JG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen_C/Sobol_JG.jld")["data"]
    ST_JU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen_C/Sobol_JU.jld")["data"]
    ST_JLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Jansen_C/Sobol_JLH.jld")["data"]

    ## Janon ##
    ST_JaS = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon_C/Sobol_JaS.jld")["data"]
    ST_JaLR = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon_C/Sobol_JaLR.jld")["data"]
    ST_JaG = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon_C/Sobol_JaG.jld")["data"]
    ST_JaU = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon_C/Sobol_JaU.jld")["data"]
    ST_JaLH = load("/home/harry/Desktop/PhD/Year 2/GSA project/Shi_Systemic/Data_for_Plots/Janon_C/Sobol_JaLH.jld")["data"]
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
        S_LVP = Vector{Float64}(undef,20)
        for i in 1:20
            for j in 1:length(t)
            S_LVP[i] = (sum(S[1:j,i]*var(sol[1,1:j])))/(sum(var(sol[1,1:j])))/length(t)
            end 
        end 
        return S_LVP 
    end 

    function dis_SAP(S)
        S_SAP = Vector{Float64}(undef,20)
        for i in 1:20
            for j in 1:length(t)
            S_SAP[i] = (sum(S[1:j,i]*var(sol[2,1:j])))/(sum(var(sol[2,1:j])))/length(t)
            end 
        end 
        return S_SAP
    end 

    function dis_LVV(S)
        S_LVV = Vector{Float64}(undef,20)
        for i in 1:20
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

    Orth_heat1 = Matrix{Float64}(undef,20,20)
    for j in 1:20
        for i in 1:20
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
    Label(ga[1, 1, Left()], "Sobol", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, hist_(ST_HS), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gb[1,1], xticks = 0.0:0.2:1.0)  #ylabel = "Lattice Rule")
    Label(gb[1, 1, Left()], "Lattice Rule", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, hist_(ST_HLR), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, font = :bold, normalization = :pdf)

    ax = Axis(gc[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:1.0:2.0)
    Label(gc[1, 1, Left()], "Golden", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, hist_(ST_HG), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(gd[1,1], xticks = 0.0:0.2:1.0)
    Label(gd[1, 1, Left()], "Uniform", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
    hist!(ax, hist_(ST_HU), color = :values, bins = 0.0:0.1:1.0, colormap= :plasma, strokewidth = 1, strokecolor = :black, normalization = :pdf)

    ax = Axis(ge[1,1], xticks = 0.0:0.2:1.0)
    Label(ge[1, 1, Left()], "Latin Hypercube", valign = :center, font = :bold, rotation = pi/2, padding = (0, 35, 0, 0))
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

    ax = Axis(gs[1,1], xticks = 0.0:0.2:1.0, yticks = 0.0:0.5:2.0)
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
