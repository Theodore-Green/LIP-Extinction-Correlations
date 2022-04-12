using MPI
MPI.Init()

size = MPI.Comm_size(MPI.COMM_WORLD)
rank = MPI.Comm_rank(MPI.COMM_WORLD)
print("Hello from $rank of $size processors!\n")

using StatGeochem, DelimitedFiles, LoopVectorization, VectorizedRNG

include("CorrelationUtilities.jl") #Functions for correlation calculations

#Listname. Choose which experiment to run by uncommenting one of the listname options below
    #All LIPs (oceanic and continetnal)
        #listname = "BigFive_with_SiberianTraps"
        #listname = "BigFive_without_SiberianTraps"
        #listname = "Stage_with_SiberianTraps"
        #listname = "Stage_without_SiberianTraps"
    #Continental LIPs only
        #listname = "Continental_BigFive_with_SiberianTraps"
        #listname = "Continental_BigFive_without_SiberianTraps"
        #listname = "Continental_Stage_with_SiberianTraps"
        #listname = "Continental_Stage_without_SiberianTraps"
    #LIPs minus the ones proposed to overlap mass extinctions (Deccan Traps, CAMP, Siberian Traps, Viluy, Kola-Dneiper)
        #listname = "Stage_minus_MassExtinctionsProvinces"
        #listname = "Continental_Stage_minus_MassExtinctionProvinces"

# Import data
if occursin("Continental", listname)
    LIPs = importdataset("data/LIPs.tsv", '\t', importas=:Tuple) #Just continental LIPs
end
if occursin("Continental", listname) == false
    LIPs = importdataset("data/LIPs Including Oceanic.tsv", '\t', importas=:Tuple) #Continental and oceanic LIPs
end

boundaries = importdataset("data/boundaries.tsv", '\t', importas=:Tuple)
TPhanerozoic = 541.0

# Based on chosen listname, selects which set of boundaries to consider
if occursin("BigFive", listname) # Big Five boundaries
    boundary_t = [23,46,53,70,84,]
end
if occursin("Stage", listname) #Phanerozoic Stage boundaries older than the Quaternary
    boundary_t = 7:101
end

# Grab the relevant data for our chosen boundary list
boundary_age = boundaries.Age_Ma[boundary_t]
boundary_age_σ = boundaries.Age_sigma_Ma[boundary_t]

# And from our chosen LIP list based on the listname selected
if occursin("without", listname)
    lip_t = (2.588 .< LIPs.Start_Age_Ma .< 541) .& (LIPs.Start_Age_Ma .!= 252.24) #Exclude Siberian Traps
end
if occursin("with_", listname)
    lip_t = (2.588 .< LIPs.Start_Age_Ma .< 541) # All Phanerozoic LIPs older than Quaternary
end
if occursin("minus", listname)
    lip_t = (2.588 .< LIPs.Start_Age_Ma .< 541) .& (LIPs.Start_Age_Ma .!= 66.36) .& (LIPs.Start_Age_Ma .!= 201.64) .& (LIPs.Start_Age_Ma .!= 252.24) .& (LIPs.Start_Age_Ma .!= 376.7) .& (LIPs.Start_Age_Ma .!= 378.64)
end
name = LIPs.Name[lip_t]
startage = LIPs.Start_Age_Ma[lip_t]
startage_σ = LIPs.Start_Age_sigma_Ma[lip_t]
endage = LIPs.End_Age_Ma[lip_t]
endage_σ = LIPs.End_Age_sigma_Ma[lip_t]

# Some relevant parameters
Δσ = 10 # How many standard deviations away from the mean do we sum up to?
nbins = 100 # Number of bins used in numerical integration / summation of distribution product
ntests = 2*10^7 # Number of uniform random tests
histmin = 0
histmax = 20
histbins = 400

# Calculate the actual observed coincidence ratio
coincidence_ratio_observed = interval_coincidence(startage, startage_σ, endage, endage_σ, boundary_age, boundary_age_σ, nbins, Δσ, TPhanerozoic)

# Calculate histogram counts for our random slice of the problem on our node
histcounts = uniform_LIP_histogram(ntests, startage, startage_σ, endage, endage_σ, boundary_age_σ, TPhanerozoic, histmin, histmax, histbins)

# MPI_Reduce to combine our results!
combined_histcounts = MPI.Reduce(histcounts, +, 0, MPI.COMM_WORLD)

# If we're on the root node, let's print our combined histogram out to file
if rank==0
    writedlm("results/histcounts_LIPs_$listname.csv", combined_histcounts, ',')
end

MPI.Finalize()


## ---
