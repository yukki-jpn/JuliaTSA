module Simulation

# Write your package code here.
    # install necessary packages (open source)
    using Statistics, Plots, Distributions, FFTW, LinearAlgebra
    
    # install necessary packages (original)
    using .FundamentalStats

    # function name
    #OK
    #export SpectraPlot
    export BarcodePlot
    export BarcodePlot!
    export BinningPlot
    export BinningPlot!
    export Correlogram
    export Correlogram!
    export VariogramPlot
    export VariogramPlot!
    export RCDFPlot
    export RCDFPlot!
    export ECDFPlot
    export ECDFPlot!
end