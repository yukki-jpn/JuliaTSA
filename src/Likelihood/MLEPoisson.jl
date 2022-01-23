"""
MLEPoisson(data::AbstractVector)
MLEPoisson estimates intensity of homogenous poisson process. Homogenous poisson process has one independent prameter, intensity $\lambda$, and the result of maximum likelihood shows following relationship
$$
\lambda = E[t_{i}-t_{i-1}].
$$
"""

function MLEPoisson(x::AbstractVector;data::AbstractVector = x)
    mu = diff(x)
    return 1/mean(mu)
end