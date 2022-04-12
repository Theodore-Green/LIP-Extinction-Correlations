## --  General math
    function within_interval(μ1::Number, σ1::Number, μ2::Number, σ2::Number, x::Number)
        return normcdf(μ1, σ1, x) * normcdf(-μ2, σ2, -x)
    end

    function normcdf_per_bin(μ::Number, σ::Number, edges::AbstractArray)
        # Allocating version
        cdf_array = normcdf(μ,σ,edges)
        return @views cdf_array[2:end] - cdf_array[1:end-1]
    end

    function normcdf_per_bin!(result::Array, μ::Number, σ::Number, edges::AbstractArray)
        # In-place version
        copyto!(result, 1, edges, 2, length(result))
        normcdf!(result, μ, σ, result)
        @inbounds for i=length(result):-1:2
            result[i] -= result[i-1]
        end
        result[1] -= normcdf(μ,σ,edges[1])
        return result
    end

    function interval_coincidence(startage, startage_σ, endage, endage_σ, boundary_age, boundary_age_σ, nbins, Δσ, ΔT)
        # Allocate required arrays
        boundary_probability = Array{Float64}(undef, nbins)
        cdfperbin = Array{Float64}(undef, nbins)
        xedges = Array{Float64}(undef, nbins+1)
        # Run non-allocating version
        return interval_coincidence!(boundary_probability, cdfperbin, xedges, startage, startage_σ, endage, endage_σ, boundary_age, boundary_age_σ, ΔT; Δσ=10, nbins=100)
    end

    function interval_coincidence!(boundary_probability, cdfperbin, xedges, startage, startage_σ, endage, endage_σ, boundary_age, boundary_age_σ, ΔT; Δσ=10, nbins=100)
        coincidence_sum = 0.0
        @inbounds for i = 1:length(startage)
            # Time range to integrate over
            t1 = endage[i] - Δσ*endage_σ[i]
            t0 = startage[i] + Δσ*startage_σ[i]
            dt = (t0-t1)/length(boundary_probability) # Timestep
            # Fill xedges vector with time bin edges
            xedges[1] = t1
            for i=2:length(xedges)
                xedges[i] = xedges[i-1] + dt
            end

            # Add up probability distributions for all the boundaries
            fill!(boundary_probability, 0.0) # Zero out the array to start
            for j = 1:length(boundary_age)
                # Coincidence contribution is virtually zero for boundaries farther away than 50 sigma
                if (t1 < boundary_age[j] < t0) || abs(boundary_age[j] - t0) < Δσ*boundary_age_σ[j] || abs(boundary_age[j] - t1) < Δσ*boundary_age_σ[j]
                    normcdf_per_bin!(cdfperbin, boundary_age[j], boundary_age_σ[j], xedges)
                    boundary_probability .+= cdfperbin
                end
            end

            # Calculate coincidence as product of boundary distribution and LIP distribution
            coincidence = 0.0
            interval_normconst = 0.0
            @simd for j = 1:length(boundary_probability)
                # Probability that a LIP is ongoing at a given time
                ongoing = within_interval(endage[i], endage_σ[i], startage[i], startage_σ[i], xedges[j]+dt/2)
                interval_normconst += ongoing
                coincidence += ongoing * boundary_probability[j]
            end

            # Add standardized coincidence factor to sum
            coincidence_sum += coincidence / interval_normconst * ΔT / dt
        end
        return coincidence_sum / length(startage) / length(boundary_age)
    end


    function gaussian_coincidence(event_age, event_age_σ, boundary_age, boundary_age_σ, ΔT)
           coincidence_sum = 0.0
           for i = 1:length(event_age)
               coincidence = 0.0
               @avx for j = 1:length(boundary_age)
                   coincidence += normproduct(event_age[i], event_age_σ[i], boundary_age[j], boundary_age_σ[j])
               end
               coincidence_sum += coincidence * ΔT
           end
           return coincidence_sum / length(event_age) / length(boundary_age)
    end

## -- Histogram functions

    # Δσ = 10 # How many standard deviations away from the mean do we sum up to?
    # nbins = 100 # Number of bins used in numerical integration / summation of distribution product

    # Define histogram function for LIPs
    function uniform_LIP_histogram(ntests, startage, startage_σ, endage, endage_σ, boundary_age_σ, ΔT, histmin, histmax, histbins; Δσ=10, nbins=100)
        # Initialize array of histogram bins, with zero counts per bin to start
        histcounts = fill(0, histbins)

        # Allocate temporary/intermediary arrays
        boundary_probability = Array{Float64}(undef, nbins)
        cdfperbin = Array{Float64}(undef, nbins)
        xedges = Array{Float64}(undef, nbins+1)
        rand_boundary_age = Array{Float64}(undef, length(boundary_age_σ))
        for n = 1:ntests
            # Draw N boundaries from a random uniform distribuion unif(0, ΔT)
            rand!(rand_boundary_age) # Fill array in-place with VectorizedRNG
            rand_boundary_age .*= ΔT # Scale appropriately
            # Calculate the coicidence ratio for these uniform random boundaries
            coincidenceₙ = interval_coincidence!(boundary_probability, cdfperbin, xedges, startage, startage_σ, endage, endage_σ, rand_boundary_age, boundary_age_σ, ΔT; Δσ=Δσ)
            histindex_float = (coincidenceₙ - histmin) * histbins / (histmax - histmin)
            if 0 < histindex_float < histbins
                histindex = ceil(Int, histindex_float)
                histcounts[histindex] += 1
            elseif histindex_float <= 0
                histcounts[1] += 1
            elseif histindex_float >= histbins
                histcounts[end] += 1
            end
        end
        return histcounts
    end

    # Define histogram function for impacts
    function uniform_impact_histogram(ntests, event_age, event_age_σ, boundary_age_σ, ΔT, histmin, histmax, histbins)
        # Initialize array of histogram bins, with zero counts per bin to start
        histcounts = fill(0, histbins)

        # Allocate temporary/intermediary arrays
        rand_boundary_age = Array{Float64}(undef, length(boundary_age_σ))
        for n = 1:ntests
            # Draw N boundaries from a random uniform distribuion unif(0, ΔT)
            rand!(rand_boundary_age) # Fill array in-place with VectorizedRNG
            rand_boundary_age .*= ΔT # Scale appropriately
            # Calculate the coicidence ratio for these uniform random boundaries
            coincidenceₙ = gaussian_coincidence(event_age, event_age_σ, rand_boundary_age, boundary_age_σ, ΔT)
            histindex_float = (coincidenceₙ - histmin) * histbins / (histmax - histmin)
            if 0 < histindex_float < histbins
                histindex = ceil(Int, histindex_float)
                histcounts[histindex] += 1
            elseif histindex_float <= 0
                histcounts[1] += 1
            elseif histindex_float >= histbins
                histcounts[end] += 1
            end
        end
        return histcounts
    end


## -- End of File
