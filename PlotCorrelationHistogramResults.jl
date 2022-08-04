## -- Load Packages and setup
    using StatGeochem, DelimitedFiles, Plots, LoopVectorization
    # Move to directory where script is present
    cd(@__DIR__)
    include("CorrelationUtilities.jl")

## --- Data about the parallel coincidence histogram results we're loading. Uncomment one file at a time to plot
    #Continental and Oceanic LIPs
        #filenamebase = "histcounts_LIPs_BigFive_with_SiberianTraps" #Mass extinction boundaries and LIPs including Siberian Traps
        #filenamebase = "histcounts_LIPs_BigFive_without_SiberianTraps" #Mass extinction boundaries and LIPs excluding Siberian Traps
        #filenamebase = "histcounts_LIPs_Stage_with_SiberianTraps" #Stage boundaries and LIPs including Siberian Traps
        #filenamebase = "histcounts_LIPs_Stage_without_SiberianTraps" #Stage boundaries and LIPs excluding Siberian Traps

    #Continental LIPs only
        #filenamebase = "histcounts_LIPs_Continental_BigFive_with_SiberianTraps" #Mass extinction boundaries and LIPs including Siberian Traps
        #filenamebase = "histcounts_LIPs_Continental_BigFive_without_SiberianTraps" #Mass extinction boundaries and LIPs excluding Siberian Traps
        #filenamebase = "histcounts_LIPs_Continental_Stage_with_SiberianTraps" #Stage boundaries and LIPs including Siberian Traps
        #filenamebase = "histcounts_LIPs_Continental_Stage_without_SiberianTraps" #Stage boundaries and LIPs excluding Siberian Traps

    #LIPs minus the ones proposed to overlap mass extinctions (Deccan Traps, CAMP, Siberian Traps, Viluy, Kola-Dneiper)
        #filenamebase = "histcounts_LIPs_Stage_minus_MassExtinctionsProvinces"
        #filenamebase = "histcounts_LIPs_Continental_Stage_minus_MassExtinctionProvinces"

    #Impacts >= 20 km
        #filenamebase = "histcounts_Impacts_BigFive_with_Chicxulub_20km" #Mass extinction boundaries and Impacts >= 20km including Chicxulub
        #filenamebase = "histcounts_Impacts_BigFive_without_Chicxulub_20km" #Mass extinction boundaries and Impacts >= 20km excluding Chicxulub
        #filenamebase = "histcounts_Impacts_Stage_with_Chicxulub_20km" #Stage boundaries and Impacts >= 20km including Chicxulub
        #filenamebase = "histcounts_Impacts_Stage_without_Chicxulub_20km" #Stage boundaries and Impacts >= 20km excluding Chicxulub

    #Impacts >= 40 km
        #filenamebase = "histcounts_Impacts_BigFive_with_Chicxulub_40km" #Mass extinction boundaries and Impacts >= 40km including Chicxulub
        #filenamebase = "histcounts_Impacts_BigFive_without_Chicxulub_40km" #Mass extinction boundaries and Impacts >= 40km excluding Chicxulub
        #filenamebase = "histcounts_Impacts_Stage_with_Chicxulub_40km" #Stage boundaries and Impacts >= 40km including Chicxulub
        #filenamebase = "histcounts_Impacts_Stage_without_Chicxulub_40km" #Stage boundaries and Impacts >= 40km excluding Chicxulub

## -- LIPs
    if occursin("LIPs", filenamebase)
    #Calculate corresponding observed coincidence
        Δσ = 10 # How many standard deviations away from the mean do we sum up to?
        nbins = 100 # Number of bins used in numerical integration / summation of distribution product
        histmin = 0
        histmax = 20
        histbins = 400

        # Import LIP and bounary datasets from file
        if occursin("Continental", filenamebase)
            LIPs = importdataset("data/LIPs.tsv", '\t', importas=:Tuple)
        else
            LIPs = importdataset("data/LIPs Including Oceanic.tsv",'\t', importas=:Tuple)
        end

        boundaries = importdataset("data/boundaries.tsv",'\t', importas=:Tuple)
        TPhanerozoic = 541.0

        #Choose boundary list (either Big Five extinctions or all Stage boundaries)
        if occursin("BigFive", filenamebase) # Big Five boundaries
            boundary_t = [23,46,53,70,84,]
            if occursin("without", filenamebase)
                listname = "BigFive_without_SiberianTraps"
            else
                listname = "BigFive_with_SiberianTraps"
            end
        else
            boundary_t = 7:101 #Phanerozoic Stage boundaries older than the Quaternary
            if occursin("without", filenamebase)
                listname = "Stage_without_SiberianTraps"
            else
                listname = "Stage_without_SiberianTraps"
            end
        end

        # Grab the relevant data for our chosen boundary list
        boundary_age = boundaries.Age_Ma[boundary_t]
        boundary_age_σ = boundaries.Age_sigma_Ma[boundary_t]

        # And from our chosen LIP list
        if occursin("without", listname) #Exclude Siberian Traps
            lip_t = (2.588 .< LIPs.Start_Age_Ma .< 541) .& (LIPs.Start_Age_Ma .!= 252.24) #Exclude Siberian Traps
        elseif occursin("minus", listname) #Exclude the LIPs assoiated with mass extinctions
                lip_t = (2.588 .< LIPs.Start_Age_Ma .< 541) .& (LIPs.Start_Age_Ma .!= 66.36) .& (LIPs.Start_Age_Ma .!= 201.64) .& (LIPs.Start_Age_Ma .!= 252.24) .& (LIPs.Start_Age_Ma .!= 376.7) .& (LIPs.Start_Age_Ma .!= 378.64)
        else occursin("with_", listname) # All Phanerozoic LIPs older than Quaternary
            lip_t = (2.588 .< LIPs.Start_Age_Ma .< 541)
        end

        name = LIPs.Name[lip_t]
        startage = LIPs.Start_Age_Ma[lip_t]
        startage_σ = LIPs.Start_Age_sigma_Ma[lip_t]
        endage = LIPs.End_Age_Ma[lip_t]
        endage_σ = LIPs.End_Age_sigma_Ma[lip_t]

        # Calculate the actual observed total coincidence product
        coincidence_ratio_observed = interval_coincidence(startage, startage_σ, endage, endage_σ, boundary_age, boundary_age_σ, nbins, Δσ, TPhanerozoic)

## --- Impacts
    elseif occursin("Impacts", filenamebase)

    #Calculate corresponding observed coincidence

        #Import impact and boundary datasets from file
        impacts = importdataset("data/impacts.tsv",'\t', importas=:Tuple)
        boundaries = importdataset("data/boundaries.tsv",'\t', importas=:Tuple)
        TPhanerozoic = 541.0

        #Choose boundary list (either Big Five extinctions or all Stage boundaries) and whether ≥20km or ≥40km impacts
        if occursin("BigFive", filenamebase)
            boundary_t = [23,46,53,70,84,] # Big Five boundaries
            if occursin("20km", filenamebase)
                if occursin("without", filenamebase)
                    listname = "Impacts_BigFive_without_Chicxulub_20km"
                else
                    listname = "Impacts_BigFive_with_Chicxulub_20km"
                end
            elseif occursin("40km", filenamebase)
                if occursin("without", filenamebase)
                    listname = "Impacts_BigFive_without_Chicxulub_40km"
                else
                    listname = "Impacts_BigFive_without_Chicxulub_40km"
                end
            end
        else
            boundary_t = 7:101 #Phanerozoic Stage boundaries older than the Quaternary
            if occursin("20km", filenamebase)
                if occursin("without", filenamebase)
                    listname = "Impacts_Stage_without_Chicxulub_20km"
                else
                    listname = "Impacts_Stage_with_Chicxulub_20km"
                end
            else occursin("40km", filenamebase)
                if occursin("without", filenamebase)
                    listname = "Impacts_Stage_without_Chicxulub_40km"
                else
                    listname = "Impacts_Stage_with_Chicxulub_40km"
                end
            end
        end

        # Grab the relevant data for our chosen boundary list
        boundary_age = boundaries.Age_Ma[boundary_t]
        boundary_age_σ = boundaries.Age_sigma_Ma[boundary_t]

        # And from our chosen impact list
        if occursin("20km", filenamebase)
            if occursin("without", filenamebase)
                impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (20 .<= impacts.Diameter .< 150) #Exclude Chicxulub
            else
                impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (20 .<= impacts.Diameter) #All impacts older than Quaternary and ≥20km
            end
        elseif occursin("40km", filenamebase)
            if occursin("without", filenamebase)
                impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (40 .<= impacts.Diameter .< 150) #Exclude Chicxulub
            else
                impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (40 .<= impacts.Diameter) #All impacts older than quaternary and ≥40km
            end
        end
        name = impacts.Name[impact_t]
        event_age = impacts.Age_Ma[impact_t]
        event_age_σ = impacts.Age_sigma_Ma[impact_t]

        # Some relevant parameters
        ntests = 2*10^7 # Number of uniform random tests
        histmin = 0
        if occursin("BigFive_with_", filenamebase)
            histmax = 150
        elseif occursin("BigFive_without_", filenamebase)
            histmax = 100
        else
            histmax = 20
        end
        histbins = 400

        # Calculate the actual observed coincidence ratio
        coincidence_ratio_observed = gaussian_coincidence(event_age, event_age_σ, boundary_age, boundary_age_σ, TPhanerozoic)
    end
## --- Plot results compared to observed coincidence

    histcounts = readdlm("results/"*filenamebase*".csv", ',', Int64)
    x = range(histmin, histmax, length=histbins)

    #Log Scale
    y = log10.(histcounts ./ (sum(histcounts) * (histmax-histmin)/histbins))
    histogram = plot(xlabel="Total coincidence product Ψ", ylabel="Log relative frequency", framestyle=:box, fg_color_legend=:white)
    plot!(histogram, x, y, label="Stochastic distribution", xlims=(histmin, histmax),color=lines[1])
    plot!(histogram, x, y, label="", fillrange=ylims()[1]*ones(size(histcounts)), fillalpha=0.05, fillcolor=lines[1], linealpha=0, ylims=ylims())
    plot!(histogram, x[x.>coincidence_ratio_observed], y[x.>coincidence_ratio_observed], label="", fillrange=ylims()[1]*ones(count(x.>coincidence_ratio_observed)), fillalpha=0.3, fillcolor=lines[1], linealpha=0, ylims=ylims())
    vline!(histogram, [coincidence_ratio_observed], label="Observed coincidence", color=:firebrick, linestyle=:dash)
    xlims!(histmin,histmax)

    #Linear scale
    #y = histcounts ./ (sum(histcounts) * (histmax-histmin)/histbins)
    #histogram = plot(xlabel="Total coincidence product Ψ", ylabel="Likelihood", framestyle=:box, fg_color_legend=:white)

    #plot!(histogram, x, y, label="Stochastic distribution", xlims=(histmin, round(coincidence_ratio_observed) + 1), ylims=(0,Inf),color=lines[1])
    #plot!(histogram, x, y, label="", fillrange=ylims()[1]*ones(size(histcounts)), fillalpha=0.8, fillcolor=lines[1], linealpha=0, ylims=ylims())
    #plot!(histogram, x[x.>coincidence_ratio_observed], y[x.>coincidence_ratio_observed], label="", fillrange=ylims()[1]*ones(count(x.>coincidence_ratio_observed)), fillalpha=0.3, fillcolor=lines[1], linealpha=0, ylims=ylims())
    #vline!(histogram, [coincidence_ratio_observed], label="Observed coincidence", color=:firebrick, linestyle=:dash)
    #xlims!(histmin,round(coincidence_ratio_observed) + 1)

    savefig(histogram, "results/"*filenamebase*".pdf")
    display(histogram)

## --- How many are greater than the observed total coincidence product?

    t = x .> coincidence_ratio_observed
    p = sum(histcounts[t])/sum(histcounts)

    one_in_every = round(Int,1/p)
    println("$filenamebase: One in every $one_in_every, p = $p")

## --- End of file
