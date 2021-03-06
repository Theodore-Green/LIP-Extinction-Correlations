## --- Load required packages
    using Plots; gr()
    using Statistics, StatsBase, StatGeochem, SpecialFunctions, LoopVectorization
    using LaTeXStrings, Formatting
    using GLM, DataFrames
    using LsqFit: curve_fit
    using Measurements # For error propagation

## --- Read datasets and define constants / subsets

    # Move to directory of script
    cd(@__DIR__)
    # Read data files
    LIPs = importdataset("data/LIPs.tsv", '\t', importas=:Tuple) #Use list of continental LIPs only
    Boundaries = importdataset("data/boundaries.tsv", '\t', importas=:Tuple)

## --- Observed Correlation between Extinction and LIP Eruption Rate. Define categories based on how two_regressionly dated LIPs are by averaging their start and end precision
    nlips = length(LIPs.Start_Age_Ma)
    nboundaries = length(Boundaries.Age_Ma)

    extinction_and_rate = plot(xlabel="Bulk eruptive rate (km³/yr)", ylabel="Extinction magnitude (%)", xlims=(0,3.7), framestyle=:box)

    # Resample Lip boundaries using widest date range
    lip_starts = LIPs.Start_Age_Ma + (2 .* LIPs.Start_Age_sigma_Ma)
    lip_ends = LIPs.End_Age_Ma - (2 .* LIPs.End_Age_sigma_Ma)

    #Boundary parameters
    boundary_age = Boundaries.Age_Ma
    boundary_age_sigma = Boundaries.Age_sigma_Ma

    LIPs_with_extinction = DataFrame(LIP_name=[], Eruption_rate=[], Eruption_rate_sigma=[], Mean_extinction_rate=[], Extinction_rate_error_lower=[], Extinction_rate_error_upper=[], boundary_number=[])

    for i = 1:nboundaries
        t = (lip_ends .< (boundary_age[i] .+ boundary_age_sigma[i])) .& ((boundary_age[i] .- boundary_age_sigma[i]) .< lip_starts)
        if any(t)
            j = argmax(LIPs.Eruption_rate_km3_yr[t])
            lip_index = findall(t)[j] # Which is the highest-eruption-rate LIP that overlaps the boundary

            #Pull relevant information about that LIP, including its duration and eruptive rate error bars
            name = LIPs.Name[lip_index] # What is the name of the LIP we're plotting
            LIP_volume = LIPs.Volume_km3[lip_index] #LIP volume in km3
            LIP_duration = LIPs.Duration_Ma[lip_index] #LIP duration in Ma
            duration_sigma = sqrt.(LIPs.Start_Age_sigma_Ma[lip_index].^2 .+ LIPs.End_Age_sigma_Ma[lip_index].^2)
            LIP_rate = LIPs.Eruption_rate_km3_yr[lip_index] #LIP eruption rate
            LIP_rate_sigma = (duration_sigma ./ LIP_duration) .* LIP_rate #1 sigma uncertainty on LIP eruption rate

            #If the error on the duration is longer than the duration (producing ∞ eruptive rate), the error on the bulk eruptive rate is set equal to the bulk eruptive rate
            #This makes the minimum eruptive rate zero and the max twice the calculated rate
            if LIP_rate_sigma .> LIP_rate
                LIP_rate_sigma = LIP_rate
            end

            #Pull relevant information about the boundary that the LIP overlaps
            Extinction_rate = Boundaries.Mean_Extinction_Rates[i] .* 100
            Extinction_rate_error_lower = Extinction_rate .- (Boundaries.Min_extinction_rates[i] .* 100)
            Extinction_rate_error_upper = (Boundaries.Max_extinction_rates[i] .* 100) .- Extinction_rate

            #LIP eruption rates or boundary extinction rates that are NaN cannot be plotted, so add the non-NaNs to data table
            if !isnan(LIP_rate) .&& !isnan(Extinction_rate) .&& !isnan(Extinction_rate_error_lower) .&& !isnan(Extinction_rate_error_upper)
                push!(LIPs_with_extinction, [name, LIP_rate, LIP_rate_sigma, Extinction_rate, Extinction_rate_error_lower, Extinction_rate_error_upper, i])
            end
        end
    end

    #Some LIPs overlap multiple boundaries, so just take the ones with the highest extinction percent
    Highest_extinction_rate = DataFrame(LIP_name=[], Eruption_rate=[], Eruption_rate_sigma=[], Mean_extinction_rate=[], Extinction_rate_error_lower=[], Extinction_rate_error_upper=[], boundary_number=[])

    for i = 1:length(LIPs_with_extinction.LIP_name)
        extinction_rates_for_given_LIP = DataFrame() #Hold the results for all the LIPs that overlap more than one boundary
        push!(extinction_rates_for_given_LIP, LIPs_with_extinction[i,1:end]) #Add all the data for the first time a LIP's name appears
        for j in 1:length(LIPs_with_extinction.LIP_name)
            if LIPs_with_extinction.LIP_name[i] == LIPs_with_extinction.LIP_name[j] #To catch LIPs with more than one boundary, separate out the list based on names that repeat
                push!(extinction_rates_for_given_LIP, LIPs_with_extinction[j,1:end])
            end
        end
        highest_rate = findmax(extinction_rates_for_given_LIP.Mean_extinction_rate) #Find highest extinction rate boundary that overlaps a given LIP
        position = highest_rate[2]
        push!(Highest_extinction_rate, extinction_rates_for_given_LIP[position,1:end]) #Knowing what row the highest extinction rate is in for a given LIP, just keep that row
    end

    Highest_extinction_rate = unique(Highest_extinction_rate) #Remove any duplicates from above proess
    #This dataframe now holds the list of the highest rate LIPs that overlap a given boundary further refined to include each LIP only once, at the highest extinction rate boundary it crosses

    #Use for curve fitting
    boundary_match_index = Array{Int64}(undef,0)
    lip_match_index = Array{Int64}(undef,0)

    #Plotting
    for i in 1:length(LIPs.Name)
        n = LIPs.Name[i]
        for j in 1:length(Highest_extinction_rate.LIP_name)
            if LIPs.Name[i] == Highest_extinction_rate.LIP_name[j] #Plot when the dataframes align, so that info can be pulled for the same LIP from each
                #Pull the relevant data out of our dataframe
                x = Highest_extinction_rate.Eruption_rate[j]
                y = Highest_extinction_rate.Mean_extinction_rate[j]
                x_error = Highest_extinction_rate.Eruption_rate_sigma[j]
                y_min = Highest_extinction_rate.Extinction_rate_error_lower[j]
                y_max = Highest_extinction_rate.Extinction_rate_error_upper[j]
                #Plot according to the precision to which the LIP is dated
                if (LIPs.Start_Age_sigma_Ma[i] .+ LIPs.End_Age_sigma_Ma[i]) / (LIPs.Start_Age_Ma[i] .+ LIPs.End_Age_Ma[i]) .< 0.002
                    # Precisely Dated LIPs. 1σ error < 0.2%
                    push!(boundary_match_index, Highest_extinction_rate.boundary_number[j])
                    push!(lip_match_index, i)
                    plot!(extinction_and_rate, [x], [y], xerror = [x_error], yerror =([y_min],[y_max]), legend=false, markerstrokecolor=:firebrick)
                    plot!(extinction_and_rate, [x], [y], seriestype=:scatter, color=:firebrick, markerstrokecolor=:firebrick, label="", markersize=4, series_annotations=[text("$n",9)])
                elseif (LIPs.Start_Age_sigma_Ma[i] .+ LIPs.End_Age_sigma_Ma[i]) / (LIPs.Start_Age_Ma[i] .+ LIPs.End_Age_Ma[i]) .< 0.0085
                    # Moderately well-dated LIPs. 1σ error < 0.85%
                    plot!(extinction_and_rate, [x], [y], xerror = [x_error], yerror =([y_min],[y_max]), legend=false, markerstrokecolor=:purple4)
                    plot!(extinction_and_rate, [x], [y], seriestype=:scatter, color=:purple4, markerstrokecolor=:purple4, label="", markersize=4, series_annotations=[text("$n",4)])
                else
                    # All other LIPs. Not well dated
                    plot!(extinction_and_rate, [x], [y], seriestype=:scatter, color=:blue, markerstrokecolor=:blue, label="", markersize=2.5)
                end
            end
        end
    end

    #Fits and Regression
        matches = Dict()
        matches["lip_rates"] = LIPs.Eruption_rate_km3_yr[lip_match_index]
        matches["extinction_percentages"] = Boundaries.Mean_Extinction_Rates[boundary_match_index] .* 100

    # Ordinary Least Squares Regression
        ols_mdl = lm(@formula(extinction_percentages ~ lip_rates), DataFrame(matches))
        line = GLM.coef(ols_mdl)[2] .* matches["lip_rates"] .+ GLM.coef(ols_mdl)[1] # Linear regression line

        # Show std error, t value, p value
        show(ols_mdl)

        # Calculate r-squared statistic (AKA coefficient of determination)
        rsquared(x,y,line_y) = 1 - sum((y .- line_y).^2)  / sum((y .- mean(y)).^2) # Function for rsquared value
        r = rsquared(matches["lip_rates"],matches["extinction_percentages"],line) # R-squared calculation
        println("\nr-squared: $r")

        # Plot regression line across a wider range of points (to see intercept)
        x_line = (0:0.1:4)
        longer_line = GLM.coef(ols_mdl)[2] .* x_line .+ GLM.coef(ols_mdl)[1] #Linear Regression line
        plot!(extinction_and_rate, x_line, longer_line, label="", line=(:dash, 1, :firebrick))

    # Logistic Regression
        # Logistic Regression Function
        function logistic(x, param)
            k,x₀,y₀ = param
            xₗ = log.(x)
            y = @. ((100 - y₀) / (1 + exp(-k*(xₗ - x₀)))) + y₀
        end

        # Fitting
        param₀ = [1.0,1.0,10.0] # Initial guess
        fitted_logistic = curve_fit(logistic, Highest_extinction_rate.Eruption_rate, Highest_extinction_rate.Mean_extinction_rate, param₀)
        fitted_logistic.param

        # Results and standard error
        println("\nLogistic fit parameters:")
        show(fitted_logistic.param)
        println("\nParameters standard error:")
        show(stderror(fitted_logistic))

        # Calculate r-squared statistic (AKA coefficient of determination)
        ŷ = logistic(Highest_extinction_rate.Eruption_rate, fitted_logistic.param)
        logistic_r = rsquared(Highest_extinction_rate.Eruption_rate, Highest_extinction_rate.Mean_extinction_rate, ŷ) #R-squared calculation
        println("\n Logistic regression r-squared: $logistic_r")

        # Plot regression line
        plot!(extinction_and_rate, 0:0.1:4, logistic(0:0.1:4, fitted_logistic.param), label="", line=(:dash, 1, :black))

    savefig(extinction_and_rate, "LIPRatevsExtinctionPercent.pdf")
    display(extinction_and_rate)

## --- Calculate predicted extinction rate from volcanism alone for precisely and moderately well-dated LIPs
    for i in 1:length(Highest_extinction_rate.LIP_name)
        name = Highest_extinction_rate.LIP_name[i]
        n = 10^5 # For Monte Carlo simulations for error propagation
        pop_volcanism_eruption_rate = rand(Normal(Highest_extinction_rate.Eruption_rate[i], Highest_extinction_rate.Eruption_rate_sigma[i]), n)
        pop_slope = rand(Normal(GLM.coef(ols_mdl)[2], GLM.stderror(ols_mdl)[2]), n)
        pop_intercept = rand(Normal(GLM.coef(ols_mdl)[1], GLM.stderror(ols_mdl)[1]), n)

        pop_predicted_extinction_percent = Vector{Float64}(undef, n)
        for j in 1:n
            pop_predicted_extinction_percent[j] = pop_intercept[rand(1:n)] + (pop_slope[rand(1:n)] * pop_volcanism_eruption_rate[rand(1:n)])
        end
        predicted_extinction_percent = mean(pop_predicted_extinction_percent) ± std(pop_predicted_extinction_percent)

        println("\nPredicted extinction magnitude from $name eruption rate (± one-sigma):\n $predicted_extinction_percent")
    end

## --- Observed Correlation between Extinction and LIP Eruption Rate, showing only the most precisely dated LIPs (1σ error < 0.2%). Same method as above, just fewer points plotted
    precise_extinction_and_rate = plot(xlabel="Bulk eruptive rate (km³/yr)", ylabel="Extinction magnitude (%)", xlims=(0.0,3.7), ylims=(0,90), framestyle=:box)

    #Use for linear regression
    boundary_match_index = Array{Int64}(undef,0)
    lip_match_index = Array{Int64}(undef,0)

    #Plotting
    for i in 1:length(LIPs.Name)
        n = LIPs.Name[i]
        for j in 1:length(Highest_extinction_rate.LIP_name)
            if LIPs.Name[i] == Highest_extinction_rate.LIP_name[j] #Plot when the dataframes align, so that info can be pulled for the same LIP from each
                #Pull the relevant data out of our dataframe
                x = Highest_extinction_rate.Eruption_rate[j]
                y = Highest_extinction_rate.Mean_extinction_rate[j]
                x_error = Highest_extinction_rate.Eruption_rate_sigma[j]
                y_min = Highest_extinction_rate.Extinction_rate_error_lower[j]
                y_max = Highest_extinction_rate.Extinction_rate_error_upper[j]
                #Plot according to the precision to which the LIP is dated
                if (LIPs.Start_Age_sigma_Ma[i] .+ LIPs.End_Age_sigma_Ma[i]) / (LIPs.Start_Age_Ma[i] .+ LIPs.End_Age_Ma[i]) .< 0.002
                    # Precisely Dated LIPs. 1σ error < 0.2%
                    push!(boundary_match_index, Highest_extinction_rate.boundary_number[j])
                    push!(lip_match_index, i)
                    plot!(precise_extinction_and_rate, [x], [y], xerror = [x_error], yerror =([y_min],[y_max]), legend=false, markerstrokecolor=:firebrick)
                    plot!(precise_extinction_and_rate, [x], [y], seriestype=:scatter, color=:firebrick, markerstrokecolor=:firebrick, label="", markersize=4, series_annotations=[text("$n",9)])
                end
            end
        end
    end

    #Fits and Regression
        matches = Dict()
        matches["lip_rates"] = LIPs.Eruption_rate_km3_yr[lip_match_index]
        matches["extinction_percentages"] = Boundaries.Mean_Extinction_Rates[boundary_match_index] .* 100

    # Ordinary Least Squares Regression
        ols_mdl = lm(@formula(extinction_percentages ~ lip_rates), DataFrame(matches))
        line = GLM.coef(ols_mdl)[2] .* matches["lip_rates"] .+ GLM.coef(ols_mdl)[1] # Linear regression line

        #Show std error, t value, p value
        show(ols_mdl)

        # Calculate r-squared statistic (AKA coefficient of determination)
        rsquared(x,y,line_y) = 1 - sum((y .- line_y).^2)  / sum((y .- mean(y)).^2) #Function for rsquared value
        r = rsquared(matches["lip_rates"],matches["extinction_percentages"],line) #Rsquared calculation
        println("\nr-squared: $r")

        # Plot regression line across a wider range of points (to see intercept)
        x_line = (0:0.1:4)
        longer_line = GLM.coef(ols_mdl)[2] .* x_line .+ GLM.coef(ols_mdl)[1]
        plot!(precise_extinction_and_rate, x_line, longer_line, label="", line=(:dash, 1, :gray))

    savefig(precise_extinction_and_rate, "PreciseLIPRatevsExtinctionPercent.pdf")
    display(precise_extinction_and_rate)

## --- Run regression line through <1km3/yr points separately from >1km3/yr
    nlips = length(LIPs.Start_Age_Ma)
    nboundaries = length(Boundaries.Age_Ma)

    two_regression_extinction_and_rate = plot(xlabel="Bulk eruptive rate (km³/yr)", ylabel="Extinction magnitude (%)", xlims=(0.0,3.7), ylims=(-10,90), framestyle=:box)

    #Use for linear regression
    boundary_match_index_lowrate = Array{Int64}(undef,0)
    lip_match_index_lowrate = Array{Int64}(undef,0)
    boundary_match_index_highrate = Array{Int64}(undef,0)
    lip_match_index_highrate = Array{Int64}(undef,0)

    for i in 1:length(LIPs.Name)
        n = LIPs.Name[i]
        for j in 1:length(Highest_extinction_rate.LIP_name)
            if LIPs.Name[i] == Highest_extinction_rate.LIP_name[j] #Plot when the dataframes align, so that info can be pulled for the same LIP from each
                #Pull the relevant data out of our dataframe
                x = Highest_extinction_rate.Eruption_rate[j]
                y = Highest_extinction_rate.Mean_extinction_rate[j]
                x_error = Highest_extinction_rate.Eruption_rate_sigma[j]
                y_min = Highest_extinction_rate.Extinction_rate_error_lower[j]
                y_max = Highest_extinction_rate.Extinction_rate_error_upper[j]
                #Plot according to the rate of the precise LIPs
                if (LIPs.Start_Age_sigma_Ma[i] .+ LIPs.End_Age_sigma_Ma[i]) / (LIPs.Start_Age_Ma[i] .+ LIPs.End_Age_Ma[i]) .< 0.002
                    # Precisely Dated LIPs. 1σ error < 0.2%
                    if x .< 1.0
                        push!(boundary_match_index_lowrate, Highest_extinction_rate.boundary_number[j])
                        push!(lip_match_index_lowrate, i)
                        plot!(two_regression_extinction_and_rate, [x], [y], xerror = [x_error], yerror =([y_min],[y_max]), legend=false, markerstrokecolor=:orange)
                        plot!(two_regression_extinction_and_rate, [x], [y], seriestype=:scatter, color=:orange, markerstrokecolor=:orange, label="", markersize=4, series_annotations=[text("$n",9)])
                    elseif x .>= 1.0
                        push!(boundary_match_index_highrate, Highest_extinction_rate.boundary_number[j])
                        push!(lip_match_index_highrate, i)
                        plot!(two_regression_extinction_and_rate, [x], [y], xerror = [x_error], yerror =([y_min],[y_max]), legend=false, markerstrokecolor=:blue2)
                        plot!(two_regression_extinction_and_rate, [x], [y], seriestype=:scatter, color=:blue2, markerstrokecolor=:blue2, label="", markersize=4, series_annotations=[text("$n",9)])
                    end
                end
            end
        end
    end

    #Fits and Regression
        #Low rate regression
            lowrate_matches = Dict()
            lowrate_matches["lip_rates"] = LIPs.Eruption_rate_km3_yr[lip_match_index_lowrate]
            lowrate_matches["extinction_percentages"] = Boundaries.Mean_Extinction_Rates[boundary_match_index_lowrate] .* 100

            # Ordinary Least Squares Regression
            ols_mdl = lm(@formula(extinction_percentages ~ lip_rates), DataFrame(lowrate_matches))
            line = GLM.coef(ols_mdl)[2] .* lowrate_matches["lip_rates"] .+ GLM.coef(ols_mdl)[1] # Linear regression line

            #Show std error, t value, p value
            show(ols_mdl)
            slope_lowrate = GLM.coef(ols_mdl)[2] ± GLM.stderror(ols_mdl)[2]
            intercept_lowrate = GLM.coef(ols_mdl)[1] ± GLM.stderror(ols_mdl)[1]
            println("\nLow Rate slope = $slope_lowrate, intercept = $intercept_lowrate")

            # Calculate r-squared statistic (AKA coefficient of determination)
            rsquared(x,y,line_y) = 1 - sum((y .- line_y).^2)  / sum((y .- mean(y)).^2) #Function for rsquared value
            r = rsquared(lowrate_matches["lip_rates"],lowrate_matches["extinction_percentages"],line) #Rsquared calculation
            println("\nLow Rate LIPs r-squared: $r")

            # Plot regression line across a wider range of points (to see intercept)
            x_line = (0:0.1:4)
            longer_line = GLM.coef(ols_mdl)[2] .* x_line .+ GLM.coef(ols_mdl)[1]
            plot!(two_regression_extinction_and_rate, x_line, longer_line, label="", line=(:dash, 1, :orange))
        #High rate regression
            highrate_matches = Dict()
            highrate_matches["lip_rates"] = LIPs.Eruption_rate_km3_yr[lip_match_index_highrate]
            highrate_matches["extinction_percentages"] = Boundaries.Mean_Extinction_Rates[boundary_match_index_highrate] .* 100

            # Ordinary Least Squares Regression
            ols_mdl = lm(@formula(extinction_percentages ~ lip_rates), DataFrame(highrate_matches))
            line = GLM.coef(ols_mdl)[2] .* highrate_matches["lip_rates"] .+ GLM.coef(ols_mdl)[1] # Linear regression line

            #Show std error, t value, p value
            show(ols_mdl)
            slope_highrate = GLM.coef(ols_mdl)[2] ± GLM.stderror(ols_mdl)[2]
            intercept_highrate = GLM.coef(ols_mdl)[1] ± GLM.stderror(ols_mdl)[1]
            println("\nHigh Rate slope = $slope_highrate, intercept = $intercept_highrate")

            # Calculate r-squared statistic (AKA coefficient of determination)
            rsquared(x,y,line_y) = 1 - sum((y .- line_y).^2)  / sum((y .- mean(y)).^2) #Function for rsquared value
            r = rsquared(highrate_matches["lip_rates"],highrate_matches["extinction_percentages"],line) #Rsquared calculation
            println("\nHigh Rate LIPs r-squared: $r")

            # Plot regression line across a wider range of points (to see intercept)
            x_line = (0:0.1:4)
            longer_line = GLM.coef(ols_mdl)[2] .* x_line .+ GLM.coef(ols_mdl)[1]
            plot!(two_regression_extinction_and_rate, x_line, longer_line, label="", line=(:dash, 1, :blue2))

    savefig(two_regression_extinction_and_rate, "TwoRegressionLIPRatevsExtinctionPercent.pdf")
    display(two_regression_extinction_and_rate)

## --- End of File
