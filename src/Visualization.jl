module Visualization

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

    function LinearFilterPlot()
    end

    function BarcodePlot(t::Vector=[nothing], c="black", d::Dict = Dict(); data::Vector = t,color=c,kwargs::Dict = d)
        r = plot([data[1],data[1]],[0,1],linestyle=:solid,
            color=color,
            size=(3000,300),
            ylabel=nothing,
            ytick=nothing,
            xlim = [0,t[end]]
            ;kwargs)
        for i in data[2:end]
            plot!(r,[i,i],[0,1],linestyle=:solid,
            label=nothing,
            color=color,
            size=(3000,300),
            ylabel=nothing,
            ytick=nothing,
            xlim = [0,t[end]]
            ;kwargs
            )
        end
        return r
    end

    function BarcodePlot!(t::Vector=[nothing], c="black", d::Dict = Dict(); data::Vector = t,color=c,kwargs::Dict = d)
        r = plot!([data[1],data[1]],[0,1],linestyle=:solid,
            color=color,
            size=(3000,300),
            ylabel=nothing,
            ytick=nothing,
            xlim = [0,t[end]]
            ;kwargs)
        for i in data[2:end]
            plot!(r,[i,i],[0,1],linestyle=:solid,
            label=nothing,
            color=color,
            size=(3000,300),
            ylabel=nothing,
            ytick=nothing,
            xlim = [0,t[end]]
            ;kwargs
            )
        end
    end

    function Binningplot(X::Vector=[nothing], Y::Vector=[nothing], N = sqrt(length(X)), k::Dict = Dict() ;
        x::Vector=X, y::Vector=Y, num = N, kwargs::Dict = k )
        a,b,c = Binning(x = x, y = y, num = N)
        return plot(a,b,yerror = c;kwargs)
    end

    function Binningplot!(X::Vector=[nothing], Y::Vector=[nothing], N = sqrt(length(X)), k::Dict = Dict() ;
        x::Vector=X, y::Vector=Y, num = N, kwargs::Dict = k )
        a,b,c = Binning(x = x, y = y, num = N)
        plot!(a,b,yerror = c;kwargs)
    end
    
    function Correlogram(x::Vector = [nothing], l::Int64 = 10, k::Dict = Dict(); 
                            data::Vector = x, lag::Int64 = l, kwargs::Dict = k)
        x::Vector{Float64} = collect(0:lag)
        y::Vector{Float64} = AutoCorrelation(data = data,lag = lag+1)
        return plot(x, y; kwargs)
    end

    function Correlogram!(x::Vector = [nothing], l::Int64 = 10, k::Dict = Dict(); 
        data::Vector = x, lag::Int64 = l, kwargs::Dict = k)
        x::Vector{Float64} = collect(0:lag)
        y::Vector{Float64} = AutoCorrelation(data = data,lag = lag+1)
        plot!(x, y; kwargs)
    end

    function VariogramPlot(x::Vector = [nothing], l::Int64 = 10, k::Dict = Dict(); 
                            data::Vector=x, lag::Int64 = l, kwargs::Dict = k)
        x::Vector{Float64},y::Vector{Float64} = Variogram(data = data,lag = lag)
        return plot(x, y; kwargs)
    end

    function VariogramPlot!(x::Vector = [nothing], l::Int64 = 10,k::Dict = Dict();
                                 data::Vector=x,lag::Int64 = l, kwargs::Dict = k)
        x::Vector, y::Vector = Variogram(data = data,lag = lag)
        plot!(x, y; kwargs)
    end

    function RCDFPlot(x::Vector = [nothing], b::Bool = true, k::Dict = Dict(); 
                                data::Vector = x,complementary::Bool = b, kwargs::Dict= k)
        x::Vector, y::Vector = RCDF(data = data,complementary = complementary)
        return plot(x, y; kwargs)
    end

    function RCDFPlot!(x::Vector = [nothing], b::Bool = true, k::Dict = Dict(); 
                                data::Vector = x,complementary::Bool = b, kwargs::Dict = k)
        x::Vector, y::Vector = RCDF(data = data, complementary = complementary)
        plot!(x, y; kwargs)
    end

    function ECDFPlot(x::Vector = [nothing], b::Bool = true, k::Dict = Dict();
                        data::Vector = x, complementary::Bool = b, kwargs::Dict = k)
        x::Vector, y::Vector = ECDF(data = data, complementary = complementary)
        return plot(x, y; kwargs)
    end
    
    function ECDFPlot!(x::Vector = [nothing], b::Bool = true, k::Dict = Dict();
                        data::Vector=x, complementary::Bool = b, kwargs::Dict = k)
        x::Vector, y::Vector = ECDF(data = data, complementary = complementary)
        plot!(x, y; kwargs)
    end
end
