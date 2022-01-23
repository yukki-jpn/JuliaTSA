module PointProcess

# Write your package code here.
    # install necessary packages (open source)
    using Statistics, Plots, Distributions , LinearAlgebra
    
    # install necessary packages (original)
    using .FundamentalStats
    export EmpiricalFireRate
    export MakeWindow
    export FanoFactor
    function MakeWindow(t::Vector=[nothing],s = 0., f = nothing , delta = 1.;
                    data::Vector = t, start = s, finish = f , dt = delta)::Taple{Vector{Float64},Vector{Int64}}
        if finish == nothing
            finish = t[end]
        end
        T::Vector{Float64} = collect(start:dt:finish)
        N::Vector{Int64} = [sum( (data.>= i) .& (data.< i+dt) ) for i in T]
        return T, N 
    end

    function EmpiricalFireRate(t::Vector=[nothing],s = 0., f = nothing , delta = 1.;
                    data::Vector = t, start = s, finish = f , dt = delta)::Taple{Vector{Float64},Vector{Float64}}
        X::Vector{Float64},Y::Vector{Float64} = MakeWindow(data = data, start = start, finish = finish, dt = dt)
        return X, Y./dt
    end
    function FanoFactor(t::Vector = [nothing], s = 0., f = nothing, delta = 1.;  
                            data::Vector = t, start = s, finish = f , dt = delta)::Float64
        _, Y::Vector{Float64} = MakeWindow(data = data, start = start, finish = finish, dt = dt)
        return var(Y)/mean(Y)
    end
    
end
