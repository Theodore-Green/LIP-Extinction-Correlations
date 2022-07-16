using MPI
MPI.Init()

size = MPI.Comm_size(MPI.COMM_WORLD)
rank = MPI.Comm_rank(MPI.COMM_WORLD)
print("Hello from $rank of $size processors!\n")

using StatGeochem, DelimitedFiles, LoopVectorization, VectorizedRNG

include("CorrelationUtilities.jl") #Functions for correlation calculations

#Listname. Choose which experiment to run by uncommenting one of the listname options below
    #Impacts ≥20km
        #listname = "Impacts_BigFive_with_Chicxulub_20km"
        #listname = "Impacts_BigFive_without_Chicxulub_20km"
        #listname = "Impacts_Stage_with_Chicxulub_20km"
        #listname = "Impacts_Stage_without_Chicxulub_20km"
    #Impacts ≥40km
        #listname = "Impacts_BigFive_with_Chicxulub_40km"
        #listname = "Impacts_BigFive_without_Chicxulub_40km"
        #listname = "Impacts_Stage_with_Chicxulub_40km"
        #listname = "Impacts_Stage_without_Chicxulub_40km"

# Import data
impacts = importdataset("data/impacts.tsv",'\t', importas=:Tuple)
boundaries = importdataset("data/boundaries.tsv",'\t', importas=:Tuple)
TPhanerozoic = 541.0

# Select the impacts and boundaries we want

# Based on chosen listname, selects which set of boundaries to consider
if occursin("BigFive", listname) # Big Five boundaries
    boundary_t = [23,46,53,70,84,]
end
if occursin("Stage", listname) # All Phanerozoic Stage boundaries
    boundary_t = 7:101
end

# Grab the relevant data for our chosen boundary list
boundary_age = boundaries.Age_Ma[boundary_t]
boundary_age_σ = boundaries.Age_sigma_Ma[boundary_t]

# And from our chosen impact list. Only include impacts that have precise radiometric ages, not ones listed without uncertainties
if occursin("20km", listname)
    if occursin("without", listname)
        impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (20 .<= impacts.Diameter .< 150) #Exclude Chicxulub
    end
    if occursin("with_", listname)
        impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (20 .<= impacts.Diameter) #All impacts older than Quaternary and ≥20km
    end
end
if occursin("40km", listname)
    if occursin("without", listname)
        impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (40 .<= impacts.Diameter .< 150) #Exclude Chicxulub
    end
    if occursin("with_", listname)
        impact_t = (2.588 .< impacts.Age_Ma .< 541) .& (40 .<= impacts.Diameter) #All impacts older than quaternary and ≥ 40km
    end
end
name = impacts.Name[impact_t]
event_age = impacts.Age_Ma[impact_t]
event_age_σ = impacts.Age_sigma_Ma[impact_t]

# Some relevant parameters for number of tests and min/max of histogram
ntests = 2*10^7 # Number of uniform random tests
histmin = 0
if occursin("BigFive_with_", listname)
    histmax = 150 #Observed coincidence products for Big Five with Chicxulub are higher than other lists
end
if occursin("Stage", listname)
    histmax = 20
end
if occursin("BigFive_without", listname)
    histmax = 100
end
histbins = 400

# Calculate the actual observed coincidence ratio
coincidence_ratio_observed = gaussian_coincidence(event_age, event_age_σ, boundary_age, boundary_age_σ, TPhanerozoic)

# Calculate histogram counts for our random slice of the problem on our node
histcounts = uniform_impact_histogram(ntests, event_age, event_age_σ, boundary_age_σ, TPhanerozoic, histmin, histmax, histbins)

# MPI_Reduce to combine our results!
combined_histcounts = MPI.Reduce(histcounts, +, 0, MPI.COMM_WORLD)

# If we're on the root node, let's print our combined histogram out to file
if rank==0
    writedlm("results/histcounts_$listname.csv", combined_histcounts, ',')
end

MPI.Finalize()


## ---
