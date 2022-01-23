"""
CV(data::AbstractVector, TimeWindow::Int64 = nothing)

CV function calculate coefficient of validation, which is defined by 
$$
CV = \frac{\sqrt{\langle \tau \rangle}}{\langle \tau \rangle}
$$
"""
function CV(x::AbstractVector,tw::Float64=nothing; data::AbstractVector = x, TimeWindow::Float64= tw)
    if TimeWindow==nothing
        return std(data)/mean(data)
    else
        D::Vector = []
        return std(D)/mean(D)
    end
end
