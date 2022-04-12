using StatGeochem, StatsBase, Plots, HypothesisTests, Distributions, LoopVectorization, VectorizedRNG

## --- Plot distribution of stage boundaries in the Phanerozoic
    cd(@__DIR__)
    Boundaries = importdataset("data/boundaries.tsv",'\t', importas=:Tuple)
    t = 7:101 #Only consider boundaries older than the Quaternary

    binedges = 2.588:10:541
    h = fit(Histogram, Boundaries.Age_Ma[t], binedges)
    counts = h.weights
    bar(cntr(binedges), counts, barwidths=binedges[2]-binedges[1], linealpha=0, label="")
    stages_histogram = plot!(xlabel="Age [Ma]", ylabel="Number of boundaries", framestyle=:box, xlims=(0,540), ylims=(0,8))
    savefig(stages_histogram, "StageBoundaryDistribution.pdf")
    display(stages_histogram)

## --- Kolmogorov–Smirnov test for the Phanerozoic stage boundary distribution
    # Draw N boundaries from a random uniform distribuion unif(0, ΔT)
    TPhanerozoic = 541.0
    rand_boundary_age = Array{Float64}(undef, length(Boundaries.Age_Ma[t]))
    rand!(rand_boundary_age) # Fill array in-place with VectorizedRNG
    rand_boundary_age .*= TPhanerozoic # Scale appropriately

    #Perform one-sample exact Kolmogorov–Smirnov test.
    #Null hypothesis: the Phanerozoic stage boundaries (excluding the Quaternary) are drawn from a uniform distribtion.
    #Alternative hypothesis the boundaries are not drawn from a uniform distribution.
    println(ExactOneSampleKSTest(Boundaries.Age_Ma[t], Uniform(Boundaries.Age_Ma[t[1]],TPhanerozoic)))

    #Perform an approximate two-sample Kolmogorov-Smirnov test.
    #Null hypothesis: the Phanerozoic stage boundaries (excluding the Quaternary) are drawn from the same distribution as a random uniform distribution of stage boundaries.
    #Alternative hypothesis: they are drawn from different distributions.
    println(ApproximateTwoSampleKSTest(Boundaries.Age_Ma[t], rand_boundary_age))

## --- End of file
