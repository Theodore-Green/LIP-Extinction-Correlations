
## --- Load required packages
    using Plots; gr()
    using Statistics, StatsBase, StatGeochem, SpecialFunctions, LoopVectorization
    using LaTeXStrings, Formatting

## --- Stats functions
    # Move to directory of script
    cd(@__DIR__)
    include("CorrelationUtilities.jl") #Functions for correlation calculations

## --- Read datasets and define constants / subsets

    LIPs = importdataset("data/LIPs Including Oceanic.tsv", '\t', importas=:Tuple)
    Boundaries = importdataset("data/boundaries.tsv", '\t', importas=:Tuple)

    # Duration of the Phanerozoic
    TPhanerozoic = 541.0

## --- Find observed coincidence products between LIPs and extinction percentage at stage boundaries

    eruption_rate = 0:0.5:3 #Eruption rate in km^3/yr for LIPs
    extinction_percentage = 0:10:70 #Percent of species going extinct at a given stage boundary
    coincidence_products_observed = Array{Float64,2}(undef, length(extinction_percentage), length(eruption_rate))

    # Some relevant parameters
    Δσ = 10 # How many standard deviations away from the mean do we sum up to?
    nbins = 100 # Number of bins used in numerical integration / summation of distribution product

    #Cycle through the above LIP eruption rate lists
    for i = 1:length(eruption_rate)

        #Grab relevant info for our chosen set of LIPs
        lip_list_observed = LIPs.Eruption_rate_km3_yr .> eruption_rate[i]
        name = LIPs.Name[lip_list_observed]
        startage = LIPs.Start_Age_Ma[lip_list_observed]
        startage_σ = LIPs.Start_Age_sigma_Ma[lip_list_observed]
        endage = LIPs.End_Age_Ma[lip_list_observed]
        endage_σ = LIPs.End_Age_sigma_Ma[lip_list_observed]


        # Cycle through the above boundary lists
        for j = 1:length(extinction_percentage)

            # Grab the relevant data for our chosen set of boundaries
            boundary_list_observed = Boundaries.Mean_Extinction_Rates .> extinction_percentage[j] / 100 # List of boundaries with extinction percent higher than any given percent in above list
            boundary_list_observed .&= Boundaries.Age_Ma .> 2.6 # Exclude very recent boundaries where hominids might mess stuff up
            boundary_age = Boundaries.Age_Ma[boundary_list_observed]
            boundary_age_σ = Boundaries.Age_sigma_Ma[boundary_list_observed]

            # Calculate the observed coincidence product
            coincidence_products_observed[j,i] = interval_coincidence(startage, startage_σ, endage, endage_σ, boundary_age, boundary_age_σ, nbins, Δσ, TPhanerozoic)
        end
    end

    ## --- Plot results (lip eruption rate vs extinction %) as a heatmap

        # X and Y axes labels
        xloc = [1,2,3,4,5,6,7]
        xlab = ["0.0","0.5","1.0","1.5","2.0","2.5","3.0"]
        yloc = [1,2,3,4,5,6,7,8]
        ylab = ["0","10","20","30","40","50","60","70"]

        #Plot results as heatmap
        heatmap_observed = heatmap(coincidence_products_observed, xlabel="Bulk eruption rate (km³/yr)", xticks=(xloc,xlab), ylabel="Mean extinction rate (%/interval)", yticks=(yloc,ylab), color = :gnuplot)
        savefig(heatmap_observed, "Observed Extinction Percent LIP Size Coincidence Products.pdf")
        display(heatmap_observed)

## --- Find max possible coincdience products between LIPs and extinction percentage at stage boundaries

    eruption_rate = 0:0.5:3 #Eruption rate in km^3/yr for LIPs
    extinction_percentage = 0:10:70 #Percent of species going extinct at a given stage boundary
    coincidence_products_max = Array{Float64,2}(undef, length(extinction_percentage), length(eruption_rate))

    Δσ = 10 # How many standard deviations away from the mean do we sum up to?
    nbins = 100 # Number of bins used in numerical integration / summation of distribution product

    #For each lip list, cycle through the extinction percentages. For each extinction percentage, sort lips and boundaries by highest precision. Assume most precise lip is coincident with most precise boundary and onward through list. Find coincidence product for each pair
    for i = 1:length(eruption_rate)

        #Grab relevant info for our chosen lip list
        lip_list_max = LIPs.Eruption_rate_km3_yr  .> eruption_rate[i]
        name = LIPs.Name[lip_list_max]
        startage = LIPs.Start_Age_Ma[lip_list_max]
        startage_σ = LIPs.Start_Age_sigma_Ma[lip_list_max]
        endage = LIPs.End_Age_Ma[lip_list_max]
        endage_σ = LIPs.End_Age_sigma_Ma[lip_list_max]

        #Cycle through each boundary
        for j = 1:length(extinction_percentage)

            # Grab the relevant data for our chosen boundary list
            boundary_list_max = Boundaries.Mean_Extinction_Rates .> extinction_percentage[j] / 100 # List of boundaries with extinction percent higher than any given percent in above list
            boundary_list_max .&= Boundaries.Age_Ma .> 5.5 # Exclude very recent boundaries where hominids might mess stuff up
            boundary_age = Boundaries.Age_Ma[boundary_list_max]
            boundary_age_σ = Boundaries.Age_sigma_Ma[boundary_list_max]

            #Sort boundaries and lips by precision of their ages
            boundary_perm = sortperm(boundary_age_σ) #sort boundary list by date precision
            LIP_precision = sqrt.(startage_σ .^ 2 + endage_σ .^ 2) .+ startage .- endage #Calculate precision of each LIP in list
            LIP_perm = sortperm(LIP_precision) # Sort lips by date precision

            coincidence_product = 0

            #For whichever _perm list is shorter, find the coincidence product for those LIP-boundary pairs
            for k = 1:min(length(boundary_perm), length(LIP_perm))

                #Grab relevant info for the LIP-boundary pairs
                name_perm = LIPs.Name[lip_list_max][LIP_perm[k]]
                startage_perm = LIPs.Start_Age_Ma[lip_list_max][LIP_perm[k]]
                startage_σ_perm = LIPs.Start_Age_sigma_Ma[lip_list_max][LIP_perm[k]]
                endage_perm = LIPs.End_Age_Ma[lip_list_max][LIP_perm[k]]
                endage_σ_perm = LIPs.End_Age_sigma_Ma[lip_list_max][LIP_perm[k]]
                boundary_age_perm = (mean(startage_perm) + mean(endage_perm))/2 #boundary age is the mean of the start and end ages of the LIP instead of the real boundary age (boundary is assumed to fall within the LIP eruption interval)
                boundary_age_σ_perm = Boundaries.Age_sigma_Ma[boundary_list_max][boundary_perm[k]]

                # Calculate the maximum coincidence product of each Lip-boundary pair
                coincidence_product += interval_coincidence(startage_perm, startage_σ_perm, endage_perm, endage_σ_perm, boundary_age_perm, boundary_age_σ_perm, nbins, Δσ, TPhanerozoic)
            end

            coincidence_products_max[j,i] = coincidence_product
        end
    end


    ## --- Plot maximum possible coincidence product results as heatmap

        # X and Y axes labels
        xloc = [1,2,3,4,5,6,7]
        xlab = ["0.0","0.5","1.0","1.5","2.0","2.5","3.0"]
        yloc = [1,2,3,4,5,6,7,8]
        ylab = ["0","10","20","30","40","50","60","70"]

        #Plot results as heatmap
        heatmap_max = heatmap(coincidence_products_max, xlabel="Bulk eruption rate (km³/yr)", xticks=(xloc,xlab), ylabel="Mean extinction rate (%/interval)", yticks=(yloc,ylab), color = :gnuplot)
        savefig(heatmap_max, "Max Possible Extinction Percent LIP Size Coincidence Product.pdf")
        display(heatmap_max)

## --- Heatmap of actual coincidence product divided by max possible

        #Ratio of calculated coincidence products to maximum possible coincidence products
        c_ratio = coincidence_products_observed ./ coincidence_products_max

        # X and Y axes labels
        xloc = [1,2,3,4,5,6,7]
        xlab = ["0.0","0.5","1.0","1.5","2.0","2.5","3.0"]
        yloc = [1,2,3,4,5,6,7,8]
        ylab = ["0","10","20","30","40","50","60","70"]

        #Plot results as heatmap
        heatmap_normalized = heatmap(c_ratio, xlabel="Bulk eruption rate (km³/yr)", xticks=(xloc,xlab), ylabel="Mean extinction rate (%/interval)", yticks=(yloc,ylab), color = :gnuplot)
        savefig(heatmap_normalized, "ObservedExtinctionPercentLIPSizeRelativeToMax.pdf")
        display(heatmap_normalized)


## --- End of File
