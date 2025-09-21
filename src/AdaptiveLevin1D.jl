module AdaptiveLevin1D

using LinearAlgebra
include("ChebyshevUtils.jl")
using .ChebyshevUtils

export adaptive_levin_1d, levin_1d


function levin_1d(f, g, a::Float64, b::Float64; k::Int=32, qr_tol_factor::Float64=1.0)
    # Chebyshev nodes + differentiation
    x_ref, D_ref = cheb_nodes_and_D(k)
    x_phys = (b - a) * (x_ref .+ 1.0) ./ 2.0 .+ a
    scale = 2.0 / (b - a)
    D = scale * D_ref
    eps0 = 1e-12

    fvals = f.(x_phys)
    gvals = g.(x_phys)
    gprime_hat = D * gvals

    # Build A and solve with pivoted QR
    A = D .+ 1im .* Diagonal(gprime_hat)
    rhs = fvals

    F = qr(A, Val(true))
    R = F.R
    piv = F.p
    diagR = abs.(diag(R))
    Anorm = maximum(diagR)
    tol_sing = Anorm * eps0 * qr_tol_factor

    pvals = zeros(ComplexF64, length(rhs))
    min_sigma = Inf
    if all(diagR .< tol_sing)
        pvals .= 0.0
    else
        l = count(≥(tol_sing), diagR)
        y = F.Q[:,1:l]' * rhs
        ytop = R[1:l,1:l] \ y
        colp = zeros(ComplexF64, length(rhs))
        colp[piv[1:l]] = ytop
        pvals .= colp
        if l > 0
            min_sigma = minimum(diagR[1:l])
        end
    end

    # boundary contribution
    idx_left = argmin(x_phys)
    idx_right = argmax(x_phys)
    g_left, g_right = gvals[idx_left], gvals[idx_right]

    shift_im = 0.5*(imag(g_left) + imag(g_right))
    shift = 1im * shift_im
    factor = exp(shift)

    term_right = pvals[idx_right] * exp(1im*(g_right - shift))
    term_left  = pvals[idx_left]  * exp(1im*(g_left  - shift))
    val = (term_right - term_left) * factor

    return - val
end

"""
    adaptive_levin_1d(f, g, a, b; tol=1e-8, k=32, qr_tol_factor=1.0, max_intervals=100000)

Compute ∫_a^b f(x) e^{i g(x)} dx adaptively using Levin collocation

- Solves local systems with pivoted QR.
- Subdivides intervals if subdivision gives a more accurate result exceeding the threshold. That is how adaptivity is achieved.
"""


function adaptive_levin_1d(f, g, a::Float64, b::Float64; tol::Float64=1e-8, k::Int=32, qr_tol_factor::Float64=1.0, max_intervals::Int=10000)

    # stack of intervals to process
    intervals = [(a, b)]
    total = 0.0 + 0.0im
    length_ab = b - a

    # safety guard
    processed = 0
    while !isempty(intervals)
        a0, b0 = pop!(intervals)
        sublength = b0 - a0
        processed += 1
        if processed > max_intervals
            println("Warning: $max_intervals intervals is not enough to $tol accuracy. The output may not be accurate.")
            return levin_1d(f, g, a, b; k=100, qr_tol_factor=qr_tol_factor)
        end

        # compute estimate on full interval
        val0 = levin_1d(f, g, a0, b0; k=k, qr_tol_factor=qr_tol_factor)

        # compute halves
        c0 = 0.5*(a0 + b0)
        valL = levin_1d(f, g, a0, c0; k=k, qr_tol_factor=qr_tol_factor)
        valR = levin_1d(f, g, c0, b0; k=k, qr_tol_factor=qr_tol_factor)

        # mixed absolute/relative criterion
        scale = max(1.0, abs(val0), abs(valL + valR))
        if abs(val0 - (valL + valR)) <= tol * scale * length_ab / sublength
            total += valL + valR
        else
            # subdivide
            push!(intervals, (a0, c0))
            push!(intervals, (c0, b0))
        end
    end

    return total
end

end