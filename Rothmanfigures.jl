## --- Load required packages
    using Plots; gr(); default(fmt=:svg)
    using Statistics, StatsBase, StatGeochem, SpecialFunctions, LoopVectorization
    using LaTeXStrings, Formatting
    using GLM, DataFrames

## --- Read datasets and define constants / subsets
    # Constants
    TPhanerozoic = 541.0 #Duration of the Phanerozoic
    Φc = 0.23 #Dimensionless constant specifying flux (from Rothman)
    Φc_error = 0.07 #Uncertainty on dimensionless constant specifying flux (from Rothman)
    τ0 = 140 #Oceanic turnover time in ky (from Rothman)
    m =  1.392*(10.0^20) #AKA m* from Rothman. Steady state value of mass of inorganic CO2 (g) in preindustrial reservoir (converted from 38,000Pg of C listed in Rothman, from Emerson and Hedges 2008)

    #Datasets
    cd(@__DIR__)
    LIPs = importdataset("data/LIPs Including Oceanic.tsv",'\t', importas=:Tuple)
    Boundaries = importdataset("data/boundaries.tsv",'\t', importas=:Tuple)

## --- LIP Rates: Rothman Figure Requirement
    #Set up dataframe to hold results of Rothman calculations
    LIP_calculations = DataFrame(LIP_name = [], Continental_or_Oceanic = [], Mass_change = [], Mass_change_low =[], Dimensionless_duration = [], Duration_error_min = [], Duration_error_max = [])

    #For each LIP, make the calculations from Rothman 2017
    for i = 1:length(LIPs.Name)
        τenv = LIPs.Duration_Ma[i] * 1000 #LIP duration in ky
        τenv_error = sqrt.(LIPs.Start_Age_sigma_Ma[i].^2 .+ LIPs.End_Age_sigma_Ma[i].^2) * 1000 #LIP duration uncertainty in ky
        LIP_volume = LIPs.Volume_km3[i] #LIP volume in km3
        Δm = (0.5/100) .* LIP_volume .* (1 * 10^9) .* (2.8 * 10^6) #Mass of inorganic carbon degassing from volcanics(g) with 0.5wt% CO2 in LIPs (taken as a reasonable value after Self et al. 2005). wt%CO2*LIP volume (convert to m3)*density of basalt(g/m3)
        Δm_low = (0.2/100) .* LIP_volume .* (1 * 10^9) .* (2.8 * 10^6) #Mass of inorganic carbon degassing from volcanics(g) with 0.2wt% CO2 in LIPs (lower part of range of 0.2-0.5 wt% CO2 in MORBS from Wallace et al. 2015) wt%CO2*LIP volume (convert to m3)*density of basalt(g/m3)
        Mass_change = abs(Δm) ./ m #Uses 0.5wt% CO2
        Mass_change_low = abs(Δm_low) ./ m #Use 0.2wt% CO2
        Dimensionless_duration = (Φc .* τenv) ./ (2 .* τ0)

        #Dimensionless duration error propagation
        n = 10^5 #Number of Monte Carlo simulations for error propagation
        pop_Φc = rand(Normal(Φc, Φc_error), n)
        pop_τenv = rand(Normal(τenv, τenv_error), n)
        pop_dimensionless_duration = Vector{Float64}(undef, n)
        for j in 1:n
            pop_dimensionless_duration[j] = (pop_Φc[rand(1:n)] .* pop_τenv[rand(1:n)]) ./ (2 .* τ0)
        end
        duration_error = std(pop_dimensionless_duration)

        #Some LIPs (like Emeishan) have durations such that, with uncertainties, it can be negative. Avoid this by setting a reasonable minimum duration of 10,000 years. All other LIPs have 1σ uncertainties
        if duration_error .> Dimensionless_duration
            duration_error_min = (Φc .* 10) ./ (2 .* τ0)
            duration_error_max = Dimensionless_duration .+ duration_error
        else
            duration_error_min = Dimensionless_duration .- duration_error
            duration_error_max = Dimensionless_duration .+ duration_error
        end

        if !isnan(Mass_change) #Only consider the LIPs that actually have volume information and therefore calculable Mass changes
            push!(LIP_calculations, [LIPs.Name[i], LIPs.Continental_or_oceanic_[i], Mass_change, Mass_change_low, Dimensionless_duration, duration_error_min, duration_error_max])
        end
    end
## --- LIP Rates: Match LIPs with corresponding extinction percents
    nlips = length(LIPs.Start_Age_Ma)
    nboundaries = length(Boundaries.Age_Ma)

    # Resample Lip boundaries using widest date range
    lip_starts = LIPs.Start_Age_Ma + (2 .* LIPs.Start_Age_sigma_Ma)
    lip_ends = LIPs.End_Age_Ma - (2 .* LIPs.End_Age_sigma_Ma)

    #Set up dataframe to hold LIPs that match with extinctions
    LIPs_with_extinction = DataFrame(LIP_name = [], Corresponding_extinction_rate = [])

    for i = 1:nboundaries
        boundary_age = (Boundaries.Age_Ma[i])
        boundary_age_sigma = (Boundaries.Age_sigma_Ma[i])
        t = (lip_ends .< (boundary_age + boundary_age_sigma)) .& ((boundary_age - boundary_age_sigma) .< lip_starts) #Do any LIPs overlap this boundary?
        if any(t) .&& !isnan(Boundaries.Mean_Extinction_Rates[i]) #If a LIP overlaps a boundary with a calculated extinction percent
            j = argmax(LIPs.Eruption_rate_km3_yr[t])
            lip_index = findall(t)[j] # Which is the highest-rate LIP that overlaps the boundary
            push!(LIPs_with_extinction, [LIPs.Name[lip_index], 100 .* Boundaries.Mean_Extinction_Rates[i]])
        end
    end

## --- LIP Rates: When multiple boundaries for same LIP, choose just the highest extinction rate one
    Highest_extinction_rate = DataFrame(LIP_name = [], Highest_extinction_rate_boundary = [])

    for i = 1:length(LIPs_with_extinction.LIP_name)
        extinction_rates_for_given_LIP = Vector{Float64}()
        push!(extinction_rates_for_given_LIP, LIPs_with_extinction.Corresponding_extinction_rate[i]) #Hold all the extinction rates for a given LIP
        for j in 1:length(LIPs_with_extinction.LIP_name)
            if LIPs_with_extinction.LIP_name[i] == LIPs_with_extinction.LIP_name[j] #To catch LIPs with more than one boundary, add all the extinction rates corresponding to that LIP to above list
                push!(extinction_rates_for_given_LIP, LIPs_with_extinction.Corresponding_extinction_rate[j])
            end
        end
        push!(Highest_extinction_rate, [LIPs_with_extinction.LIP_name[i], nanmaximum(extinction_rates_for_given_LIP)]) #Find highest extinction rate boundary that overlaps a given LIP
    end

    Highest_extinction_rate = unique(Highest_extinction_rate) #Remove any duplicated LIP-extinction pairs
## --- LIP Rates: Plots
    rothmanfigure = plot(xlabel="Dimensionless duration Φcτenv/2τ₀", ylabel="Direct Igneous CO₂ Release |Δm|/m*", ylims=(10^(-2),10^0), xlims=(10^(-1),10^(2)), colorbar_title="Extinction magnitude (%)", framestyle=:box, legend=false, colorbar=true)
    rothmanfigure_allLIPs = plot(xlabel="Dimensionless duration Φcτenv/2τ₀", ylabel="Direct Igneous CO₂ Release |Δm|/m*", ylims=(10^(-2),10^1), xlims=(10^(-1),10^(2)), colorbar_title="Extinction magnitude (%)", framestyle=:box, legend=false, colorbar=true)

    for i = 1:length(LIP_calculations.LIP_name)
        for j in 1:length(Highest_extinction_rate.LIP_name)
            if LIP_calculations.LIP_name[i] == Highest_extinction_rate.LIP_name[j]
                xerr = (([LIP_calculations.Dimensionless_duration[i] .- LIP_calculations.Duration_error_min[i]]), ([LIP_calculations.Duration_error_max[i] .- LIP_calculations.Dimensionless_duration[i]]))
                n = LIP_calculations.LIP_name[i]
                if LIP_calculations.Continental_or_Oceanic[i] == "Oceanic"  #Add to supplementary figure
                    plot!(rothmanfigure_allLIPs, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xerr = xerr, yerr=([LIP_calculations.Mass_change[i] .- LIP_calculations.Mass_change_low[i]], [0]), xscale=:log10, yscale=:log10, seriestype=:scatter, markershape=:square)
                    plot!(rothmanfigure_allLIPs, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xscale=:log10, yscale=:log10, seriestype=:scatter, marker_z=Highest_extinction_rate.Highest_extinction_rate_boundary[j], color=:lighttest, series_annotations=[text("$n",3)], markershape=:square)
                else #Add all other points to main text plot and supplmental
                    plot!(rothmanfigure, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xerr = xerr, yerr=([LIP_calculations.Mass_change[i] .- LIP_calculations.Mass_change_low[i]], [0]), xscale=:log10, yscale=:log10, seriestype=:scatter)
                    plot!(rothmanfigure, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xscale=:log10, yscale=:log10, seriestype=:scatter, marker_z=Highest_extinction_rate.Highest_extinction_rate_boundary[j], color=:lighttest, series_annotations=[text("$n",3)])
                    plot!(rothmanfigure_allLIPs, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xerr = xerr, yerr=([LIP_calculations.Mass_change[i] .- LIP_calculations.Mass_change_low[i]], [0]), xscale=:log10, yscale=:log10, seriestype=:scatter)
                    plot!(rothmanfigure_allLIPs, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xscale=:log10, yscale=:log10, seriestype=:scatter, marker_z=Highest_extinction_rate.Highest_extinction_rate_boundary[j], color=:lighttest, series_annotations=[text("$n",3)])
                end
            end
        end
    end

#Critical thresholds for eruptions
    #From Rothman Figure 2, "the straight (identity) line denotes the equality predicted by Eq. 5 when Φ = Φc = 0.23 ± 0.07." So plot identity line as critical threshold 100% then multiply by various factors to get the other lines
    x = 0.01:0.01:500
    #Main text
    plot!(rothmanfigure, x, x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[x; reverse(x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfigure, x, 0.1 .* x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[0.1 .* x; reverse(0.1 .* x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfigure, x, 0.01 .* x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[0.01 .* x; reverse(0.01 .* x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfigure, x, 0.001 .* x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[0.001 .* x; reverse(0.001 .* x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)

    #Supplemental
    plot!(rothmanfigure_allLIPs, x, x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure_allLIPs, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[x; reverse(x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfigure_allLIPs, x, 0.1 .* x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure_allLIPs, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[0.1 .* x; reverse(0.1 .* x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfigure_allLIPs, x, 0.01 .* x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure_allLIPs, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[0.01 .* x; reverse(0.01 .* x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfigure_allLIPs, x, 0.001 .* x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfigure_allLIPs, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[0.001 .* x; reverse(0.001 .* x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)

    #Show and save plots
    display(rothmanfigure)
    savefig(rothmanfigure, "CFBRothmanFigure.pdf")

    display(rothmanfigure_allLIPs)
    savefig(rothmanfigure_allLIPs, "RothmanFigureAll.pdf")

## --- Rothman Figure with Anthropogenic CO2 inputs (based on IPCC report 2021)
    #Anthropogenic Forcings
    emission_scenario = ["present","1.9", "2.6", "4.5", "7.0", "8.5"]
    oceanic_carbon_uptake = [637, 1020, 1156, 1472, 1743, 1973] #Pg of carbon dioxide uptake in each emission scenario from IPCC report
    Mass_change_modern = Vector{Float64}(undef,length(emission_scenario))
    Dimensionless_duration_modern = Vector{Float64}(undef,length(emission_scenario))
    Dimensionless_duration_modern_error_min = Vector{Float64}(undef,length(emission_scenario))
    Dimensionless_duration_modern_error_max = Vector{Float64}(undef,length(emission_scenario))

    for i = 1:length(emission_scenario)
        if i == 1
            τenv_modern = 0.171 #Present(in ky years since 1850. Taken as 2021 because that is year of IPCC report values). Minimum effective duration is 10ky
        else
            τenv_modern = 0.250 #2100 (in ky since 1850). Minimum effective duration is 10ky
        end
        Δm = oceanic_carbon_uptake[i] .* 10^15 #Mass of inorganic carbon dioxide (g) from anthropogenic emissions into the ocean
        Mass_change_modern[i] = abs(Δm) ./ m  #Mass change and convert to g
        Dimensionless_duration_modern[i] = (Φc .* τenv_modern) ./ (2 .* τ0)
        #Duration error propagation
        n = 10^5 #Number of Monte Carlo simulations for error propagation
        pop_Φc = rand(Normal(Φc, Φc_error), n)
        pop_dimensionless_duration_modern = Vector{Float64}(undef, n)
        for j in 1:n
            pop_dimensionless_duration_modern[j] = (pop_Φc[rand(1:n)] .* τenv_modern) ./ (2 .* τ0)
        end
        Dimensionless_duration_modern_error = std(pop_dimensionless_duration_modern)
        Dimensionless_duration_modern_error_min[i] = Dimensionless_duration_modern[i] .- Dimensionless_duration_modern_error
        Dimensionless_duration_modern_error_max[i] = Dimensionless_duration_modern[i] .+ Dimensionless_duration_modern_error
    end

    #Plot the LIP-extinction pairs frm the main text figure
    rothmanfiguremodern = plot(xlabel="Dimensionless duration Φcτenv/2τ₀", ylabel="Direct Igneous CO₂ Release |Δm|/m*", ylims=(10^(-2.5),10^0), xlims=(10^(-4),10^(2)), framestyle=:box, legend=false)
    for i = 1:length(LIP_calculations.LIP_name)
        for j in 1:length(Highest_extinction_rate.LIP_name)
            if LIP_calculations.LIP_name[i] == Highest_extinction_rate.LIP_name[j]
                xerr = (([LIP_calculations.Dimensionless_duration[i] .- LIP_calculations.Duration_error_min[i]]), ([LIP_calculations.Duration_error_max[i] .- LIP_calculations.Dimensionless_duration[i]]))
                n = LIP_calculations.LIP_name[i]
                if LIP_calculations.Continental_or_Oceanic[i] !== "Oceanic"  #Add to supplementary figure
                    plot!(rothmanfiguremodern, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xerr = xerr, yerr=([LIP_calculations.Mass_change[i] .- LIP_calculations.Mass_change_low[i]], [0]), xscale=:log10, yscale=:log10, seriestype=:scatter)
                    plot!(rothmanfiguremodern, [LIP_calculations.Dimensionless_duration[i]], [LIP_calculations.Mass_change[i]], xscale=:log10, yscale=:log10, seriestype=:scatter, color=:orange, series_annotations=[text("$n",3)])
                end
            end
        end
    end

    #Plot the modern and projected emission scenarios from IPCC report
    for i = 1:length(Mass_change_modern)
        e = emission_scenario[i]
        xerr = (([Dimensionless_duration_modern[i] .- Dimensionless_duration_modern_error_min[i]]), ([Dimensionless_duration_modern_error_max[i] .- Dimensionless_duration_modern[i]]))
        if i == 1
            plot!(rothmanfiguremodern, [Dimensionless_duration_modern[i]], [Mass_change_modern[i]], xerr = xerr, xscale=:log10, yscale=:log10, seriestype=:scatter, color=:green1, series_annotations=[text("Present",3)])
        else
            plot!(rothmanfiguremodern, [Dimensionless_duration_modern[i]], [Mass_change_modern[i]], xerr = xerr, xscale=:log10, yscale=:log10, seriestype=:scatter, color=:blue, series_annotations=[text("SSP-$e",3)])
        end
    end

    #Critical threshold for durations above 10ky (minimum effective uptake, which has dimensionless duration of 0.082)
    x = 0.001:0.0001:60 #Overlap the other critical threshold section so that the uncertainty lines can all intersect
    plot!(rothmanfiguremodern, x, x, label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfiguremodern, [x*1/(.3/.23); reverse(x*1/(.16/.23))],[x; reverse(x)], linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)

    #Critical threshold for durations below 10ky (minimum effective uptake, which has dimensionless duration of 0.082). Modern critical size for changes to the ocean carbon dioxide resevoir (310±155 PgC (1140+/570 Pg CO2) is modern critical threshold)
    x2 = 0.0001:0.0001:0.01 #Overlap the other critical threshold section so that the uncertainty lines can all intersect
    #Modern critical size of perturbation of the marine CO2 cycle for durations below 10ky is M_c = 0.0082±0.0041
    plot!(rothmanfiguremodern, x2, 0.0082*ones(length(x2)), label="", xscale=:log10, yscale=:log10, series_type=:line, color=:black)
    plot!(rothmanfiguremodern, x2, 0.0123*ones(length(x2)), linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)
    plot!(rothmanfiguremodern, x2, 0.0041*ones(length(x2)), linealpha = 0.3, series_type=:line, xscale=:log10, yscale=:log10, color=:black)

    display(rothmanfiguremodern)
    savefig(rothmanfiguremodern, "RothmanFigurewithModernInputs.pdf")

## --- End of File
