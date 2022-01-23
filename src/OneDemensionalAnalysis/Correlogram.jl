"""
Correlogram(data::AbstractVector, lag::Int64 = 20, fft::Bool = true)


Correlogram function calculates the autocorrelation function of given data. To use this software, you must input your dataset into a 'data' argument. You have two options to calculate autocorrelation, which use Fourier fast transformation (FFT) and calculate by its definition. We calculate the autocorrelation function by lag if you set 'lag' arguments. If you set 'fft' argument true, we calculate the autocorrelation function by FFT. Otherwise, we calculate them by definition.
"""
function Correlogram(d::AbstractVector, l::Int64 = 20, f::Bool = true ;data::AbstractVector = d, lag::Int64 = l, fft::Bool = f)
    n::Int64 = length(data)
    s::Float64 = zeros(Float64, 2*n-1 )
    s[1:n] .= data
    Fs = real.(ifft(abs2.(fft(s))))
    return Fs
end