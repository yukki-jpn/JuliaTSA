module FundamentalStats

# Write your package code here.
    using Statistics, FFTW, Distributions
    # about poisson process
    export RCDF
    export ECDF
    export PowerSpectrum
    export AutoCorrelation
    export Variogram
    export Binning
    export Smooth
    export Q_Quantile
    export GaussFilter
    export StepFilter
        
    function StepFilter(t::Vector=[nothing],s = 1.; 
                            data::Vector = t, stepsize = s)
        
        x2 = [i for i in 0:(1/stepsize):t[end]]
        y2 = [sum( (data.<i+stepsize/2) .& (data.>i-stepsize/2) )/stepsize for i in x2]
        return x2,y2
    end
    

    function Gaussian(t::Vector=[nothing],σ = 1.; 
                            data::Vector = t, sigma = σ)::Vector{Float64}
        return (1/(sqrt(2*π)*sigma) ) .* exp.(-(data.^2)./(2 *sigma^2))
    end

    function GaussFilter(t::Vector = [nothing],  σ = 1.; 
                                    data::Vector = t, sigma = σ)
        result = zeros(length(data))
        for i in 1:1:length(data)
            result[i] = sum(Gaussian(data = data.-data[i],σ))
        end
        return result
    end


    function Q_Quantile(x = nothing,y = Normal(0,1), num = 100; data_x = x, data_y = y)::Taple{Vector{Float64},Vector{Float64}}
        t = collect(0:(1/num):1)
        Q1 = quantile(data_x,t)
        Q2 = quantile(data_y,t)
        return Q1, Q2
    end

    function MakeWindowFilter(t::Vector,start = 0., finish = nothing,dt = 1., filter=Gaussian, params::Vector = [1.])
        if finish == nothing
            finish = t[end]
        end
        T::Vector{Float64} = [i for i in start:dt:finish]
        N::Vector{Float64} = [ sum(filter(t.-i, params...)) for i in T]
        return T, N
    end


    function Binning(X::Vector=[nothing], Y::Vector=[nothing], N = sqrt(length(X)) ;
         x::Vector=X, y::Vector=Y, num = N )::Tuple{Vector{Float64},Vector{Float64},Vector{Float64} }
        p::Vector{Int64} = sortperm(x)
        x_::Vector{Float64} = x[p]
        y_::Vector{Float64} = y[p]
        step::Float64 = (x_[end] - x_[1])/(num)
        start::Float64 = x_[1] + step
        x_save::Vector{Float64} = zeros(Float64,num+1)
        y_save::Vector{Float64} = zeros(Float64,num+1)
        y_err::Vector{Float64}  = zeros(Float64,num+1)
        P::Vector{Bool} = zeros(Bool,num+1)
        n::Int64 = 1
        Y::Float64 = 0
        Y2::Float64 = 0
        N::Int64 = 0
        for i in 1:length(x)
            while x_[i] > start
                x_save[n] = start - step/2
                y_save[n] = Y/N
                y_err[n] = sqrt( Y2/(N)-(Y/N)^2 )/sqrt(N)
                P[n] = (N!=0)
                Y = 0
                Y2 = 0
                start += step
                n += 1
                N=0
            end
            Y2 += y_[i]^2
            Y += y_[i]
            N += 1
        end
        x_save[n] = start - step/2
        y_save[n] = Y/N
        y_err[n] = sqrt( Y2/(N)-(Y/N)^2 )/sqrt(N)
        P[n] = (N!=0)
        return x_save[P], y_save[P],y_err[P]
    end

    function PowerSpectrum(x::Vector = [nothing], f::Bool = false;
                                data::Vector = x, full::Bool = f)::Vector{Float64}
        X = zeros(length(data)*2-1)
        X[1:length(data)] .= (data .- mean(data))
        X = abs2.(fft(X))
        if full
            return X
        else
            return X[1:length(data)]
        end
    end

    function AutoCorrelation(x::Vector = [nothing], l::Int64=10; 
                                data::Vector = x,lag::Int64 = l)::Tuple{Vector{Float64}, Vector{Float64}}
        X = PowerSpectrum(data = data, full = true)
        X = real.(ifft(X))[1:lag]
        X = X./X[1]
        return X
    end
    
    function Variogram(x::Vector = [nothing], l::Int64=10; 
                                data::Vector = x,lag::Int64 = l)::Tuple{Vector{Float64}, Vector{Float64}}
        n::Float64 = length(data)
        X::Vector{Float64} = collect(range(1,lag, length=lag))
        Y::Vector{Float64} = [sum( (data[1:end-i]-data[i+1:end]).^2 )/(n-i) for i in 1:l]
        return X, Y
    end

    # basic statistical tools
    ```
    For analyses on the stochastic process, explanatory data analysis do have crucial role. We thus provide some statistical tools, such as estimation of empirical cumulative distributions and ranked cumulative distributions.
    ```

    function RCDF(x::Vector = [nothing],b::Bool = true;
                            data::Vector = x, complementary::Bool = b)::Tuple{Vector, Vector{Float64}}
        if complementary 
            return sort(data), [i for i in 1:1:length(data)]
        else
            return sort(data), [i for i in length(data):-1:1]
        end
    end
    
    function ECDF(x::Vector = [nothing], b::Bool = true;
                            data::Vector = x, complementary::Bool = b)::Tuple{Vector{Float64}, Vector{Float64}} 
        n::Int64 = length(data)
        a::Vector = diff(sort(data)) .!= 0
        y::Vector = [i for i in 1:length(data)-1]
        num::Vector = y[a]
        if complementary
            pushfirst!(num,0.0)
            num = n.-num
            return num./n,sort(unique(data))
        else
            push!(num,length(data))
            return num./n,sort(unique(data))
        end    
    end

end