## --- Load required pacages
    using Plots; gr(); default(fmt=:svg)
    using Statistics, StatsBase, StatGeochem, SpecialFunctions
    using ProgressMeter: @showprogress
    using LaTeXStrings, Formatting

## --- Stats functions
    # Move to directory of script
    cd(@__DIR__)
    include("CorrelationUtilities.jl")

## --- K-Pg bounday and Deccan Traps
    yloc = [65.0,65.5,66.0,66.5,67.0]; ylab = ["67.0","66.5","66.0","65.5","65.0"]
    h1 = plot(xlabel="Likelihood",ylabel="Time [Ma]", xtickfontsize = 7, yticks=(yloc,ylab), ytickfontsize = 7)
    timespan1 = 65:0.001:67

    #Likelihood of LIP Eruption within interval
    plot!(h1, reverse(within_interval.(65.590,0.0135,66.358,0.0185,timespan1)), timespan1, label="Deccan Traps", color=:blue)

    #Boundary Gausian PDF
    boundary_pdf = zeros(length(timespan1))
    i = 0
    for x = timespan1
        global i += 1
        boundary_pdf[i] = normpdf(66.016,0.025,x)
    end
    plot!(h1, reverse(boundary_pdf), timespan1, label="K-Pg Boundary PDF", color=:green)

## --- Tr-J boundary and CAMP
    yloc = [200.0,200.5,201.0,201.5,202.0]; ylab = ["202.0","201.5","201.0","200.5","200.0"]
    h2 = plot(xtickfontsize = 7, yticks=(yloc,ylab), ytickfontsize = 7)
    timespan2 = 200:0.001:202

    #Likelihood of LIP Eruption within interval
    plot!(h2, reverse(within_interval.(200.916,0.032,201.639,0.0145,timespan2)), timespan2, label="CAMP", color=:blue)

    #Boundary Gausian PDF
    boundary_pdf = zeros(length(timespan2))
    i = 0
    for x = timespan2
        global i += 1
        boundary_pdf[i] = normpdf(201.36,0.085,x)
    end
    plot!(h2, reverse(boundary_pdf), timespan2, label="Tr-J Boundary PDF", color=:green)

## --- P-T boundary and Siberian Traps
    yloc = [250.5,251.0,251.5,252.0,252.5]; ylab = ["252.5","252.0","251.5","251.0","250.5"]
    h3 = plot(xtickfontsize = 7, yticks=(yloc,ylab), ytickfontsize = 7)
    timespan3 = 250.5:0.0001:252.5

    #Likelihood of LIP Eruption within interval
    plot!(h3, reverse(within_interval.(251.354,0.044,252.24,0.05,timespan3)),timespan3, label="Siberian Traps", color=:blue)

    #Boundary Gausian PDF
    boundary_pdf = zeros(length(timespan3))
    i = 0
    for x = timespan3
        global i += 1
        boundary_pdf[i] = normpdf(251.902,0.024,x)
    end
    plot!(h3, reverse(boundary_pdf), timespan3, label="PT Boundary PDF", color=:green)

## --- End Guadilupian and Emeishan
    yloc = [259.0,259.5,260.0,260.5,261.0]; ylab = ["261.0","260.5","260.0","259.5","259.0"]
    h4 = plot(xtickfontsize = 7, yticks=(yloc,ylab), ytickfontsize = 7)
    timespan4 = 259:0.001:261

    #Likelihood of LIP Eruption within interval
    plot!(h4, reverse(within_interval.(259.51,0.095,259.9,0.4,timespan4)), timespan4, label="Emeishan", color=:blue)

    #Boundary Gausian PDF
    boundary_pdf = zeros(length(timespan4))
    i = 0
    for x = timespan4
        global i += 1
        boundary_pdf[i] = normpdf(259.8,0.15,x)
    end

    plot!(h4, reverse(boundary_pdf), timespan4, label="End-Guadilupian Boundary PDF", color=:green)

## --- Frasnian-Famennian boundary and Viluy
    yloc = [360,365,370,375,380]; ylab = ["380","375","370","365","360"]
    h5 = plot(xtickfontsize = 7, yticks=(yloc,ylab), ytickfontsize = 7)
    timespan5 = 360:0.001:380

    #Likelihood of LIP Eruption within interval
    plot!(h5, reverse(within_interval.(364.4,0.85,376.7,0.85,timespan5)), timespan5, label="Viluy", color=:blue)

    #Boundary Gausian PDF
    boundary_pdf = zeros(length(timespan5))
    i = 0
    for x = timespan5
        global i += 1
        boundary_pdf[i] = normpdf(371.855,0.075,x)
    end
    plot!(h5, reverse(boundary_pdf), timespan5, label="Famennian Boundary PDF", color=:green)

## --- Combined Figure
#h6 = plot(h1,h2,h3,h4,h5,layout=(1,5),legend=:best, legendfontsize=3)
h6 = plot(h1,h2,h3,layout=(1,3),legend=:best, legendfontsize=8)
savefig(h6, "CombinedTimeLikelihood.pdf")
display(h6)


## --- End of File
