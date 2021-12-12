module StochasticPointProcess

# Write your package code here.
    using Statistics
    # about poisson process
    export SimulateStationaryPoisson
    export EstimateStationaryPoisson

    # basic statistical tools
    ```
    For analyses on the stochastic process, explanatory data analysis do have crucial role. We thus provide some statistical tools, such as estimation of empirical cumulative distributions and ranked cumulative distributions.
    ```
    function RCDF(x::Vector)::Tuple{Vector, Vector{Float64}} 
        return sort(x), [i for i in 1.0:1.0:length(x)]
    end
    
    function ECDF(x::Vector,complementary::Bool = true)::Tuple{Vector{Float64}, Vector{Float64}} 
        n::Int64 = length(x)
        a::Vector = diff(sort(x)) .!= 0
        y::Vector = [i for i in 1:length(x)-1]
        num::Vector = y[a]
        if complementary
            pushfirst!(num,0.0)
            num = n.-num
            return num./n,sort(unique(x))
        else
            push!(num,length(x))
            return num./n,sort(unique(x))
        end    
    end

    function EstimateStationaryPoisson(t::Vector, interval::Bool = false)::Float64
        if interval
            return mean(t)
        else
            return mean(diff(t))
        end
    end
    
    function EstimateStationaryPoisson(t::Array, interval::Bool = false)::Float64
        if interval
            return mean(t)
        else
            return mean(diff(t))
        end
    end
    
    function NonStationarPoiss(λ, t::Float64)::Vector{Float64}
        dt::Float64 = 1/λ(0)/100.
        T::Float64 = 0.
        r::Float64 = 0.
        result::Vector{Float64} = zeros(Float64,10^8)
        num::Int64 = 1
        while T < t
            r = rand()
            if r < λ(T) * dt
                result[num] = T
                num+=1
            end
            T += dt
            dt = 1/λ(T)/100
        end
        print(T)
        return result[1:num-1]
    end

    function TimeReshape(λ, t::Float64)::Vector{Float64}
        dt::Float64 = 1/λ(0)/100.
        T::Float64 = 0.
        r::Float64 = 0.
        result::Vector{Float64} = zeros(Float64,10^8)
        num::Int64 = 1
        while T < t
            r = rand()
            if r < λ(T) * dt
                result[num] = T
                num+=1
            end
            T += dt
            dt = 1/λ(T)/100
        end
        print(T)
        return result[1:num-1]
    end

    
    function Gaussian(t::Vector,σ)::Vector
        return (1/(2*π*σ)) .* exp.(-(t.^2)./(2 *σ^2) )
    end
    function Gaussian(t,σ)::Float64
        return (1/(2*π*σ)) * exp(-(t^2)/(2 *σ^2) )
    end
end
