"""
CV(data::AbstractVector)
CV(data::AbstractVector, TimeWindow::Int64, Start = data[1], End = data[end])
CV(data::AbstractVector, TimeWindow::Float64,Start = data[1], End = data[end])


CV function calculate coefficient of validation, which can be described with 
$$
CV = \frac{\sqrt{\langle \tau \rangle}}{\langle \tau \rangle}.
$$
CV can be underestimated if the dataset is right censored. 
"""
function CV(x::AbstractVector, tw::Float64 = nothing;
         data::AbstractVector = x, TimeWindow::Float64= tw, Start = data[1], End = data[end])
    if TimeWindow == nothing
        return std(data)/mean(data)
    else
        t::Vector{Float64} = [i for i in start:tw:End]
        CVs::Vector{Float64} = zeros(Float64,length(t))
        n::Int64 = 0
        for i in t
            a = data[ i .< data .<  i + tw ]
            if length(a) != 0
                CVs[n] = std(a)/mean(a)
            else
                CVs[n] = nothing
            end
        end
        return t,CVs
    end
end

function CV(x::AbstractVector, tw::Int64 = nothing; 
                data::AbstractVector = x, TimeWindow::Int64= tw, Start = data[1], End = data[end])
    if TimeWindow == nothing
        return std(data)/mean(data)
    else
        t::Vector{Float64} = [i for i in start:tw:End]
        CVs::Vector{Float64} = zeros(Float64,length(t))
        n::Int64 = 0
        for i in t
            a = data[ i .< data .<  i + tw ]
            if length(a) != 0
                CVs[n] = std(a)/mean(a)
            else
                CVs[n] = nothing
            end
        end
        return t,CVs
    end
end


function CV(x::AbstractVector, tw::Float64 = nothing, s,e;
                data::AbstractVector = x, TimeWindow::Float64= tw, Start = s, End = e)
    if TimeWindow == nothing
        return std(data)/mean(data)
    else
        t::Vector{Float64} = [i for i in start:tw:End]
        CVs::Vector{Float64} = zeros(Float64,length(t))
        n::Int64 = 0
        for i in t
            a = data[ i .< data .<  i + tw ]
            if length(a) != 0
                CVs[n] = std(a)/mean(a)
            else
                CVs[n] = nothing
            end
        end
        return t,CVs
    end
end

function CV(x::AbstractVector, tw::Int64 = nothing, s,e;
                data::AbstractVector = x, TimeWindow::Int64= tw, Start = s, End = e)
    if TimeWindow == nothing
        return std(data)/mean(data)
    else
        t::Vector{Float64} = [i for i in start:tw:End]
        CVs::Vector{Float64} = zeros(Float64,length(t))
        n::Int64 = 0
        for i in t
            a = data[ i .< data .<  i + tw ]
            if length(a) != 0
                CVs[n] = std(a)/mean(a)
            else
                CVs[n] = nothing
            end
        end
        return t,CVs
    end
end