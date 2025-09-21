module AdaptiveLevin2D

using LinearAlgebra
include("ChebyshevUtils.jl")
using .ChebyshevUtils
include("AdaptiveLevin1D.jl")
using .AdaptiveLevin1D

export adaptive_levin_2d, levin_2d

# Solve one column using truncated QR
function solve_column(Dx::Matrix{}, gcol::Vector{}, fcol::Vector{}; qr_tol_factor::Float64=1.0)
    gprime_hat = Dx * gcol
    Aj = Dx .+ 1im .* Diagonal(gprime_hat)
    k = length(fcol)
    eps0 = 1e-12
    
    F = qr(Aj, Val(true))  #QR
    R, piv = F.R, F.p
    diagR = abs.(diag(R))
    tol = maximum(diagR) * eps0 * qr_tol_factor
    if all(diagR .< tol)
        return zeros(ComplexF64, k)
    else
        l = count(≥(tol), diagR)
        y = F.Q[:,1:l]' * fcol
        ytop = R[1:l,1:l] \ y
        colp = zeros(ComplexF64, k)
        colp[piv[1:l]] = ytop
        return colp
    end
end

# 2D Levin via delamination
function levin_2d_general(f, g, fdomain::Vector{Float64}, gdomain::Vector{Float64}; k::Int=32, tol::Float64=1e-8, qr_tol_factor::Float64=1.0)
    a, b = fdomain
    c, d = gdomain
    
    # Chebyshev grid and differentiation
    t_ref, D_ref = cheb_nodes_and_D(k)
    x_phys = (b - a)*(t_ref .+ 1)/2 .+ a
    y_phys = (d - c)*(t_ref .+ 1)/2 .+ c

    fvals = Array{ComplexF64}(undef, k, k)
    gvals = Array{ComplexF64}(undef, k, k)
    for j in 1:k
        fvals[:, j] = f.(x_phys, y_phys[j])
        gvals[:, j] = g.(x_phys, y_phys[j])
    end

    Dx = (2.0 / (b - a)) .* D_ref
    pb_mat = zeros(ComplexF64, k, k)
    for j in 1:k
        pb_mat[:, j] = solve_column(Dx, gvals[:, j], fvals[:, j]; qr_tol_factor=qr_tol_factor)
    end
    
    # Solve columns along x (delamination)
    pb_mat = zeros(ComplexF64, k, k)
    for j in 1:k
        fcol = f.(x_phys, y_phys[j])
        gcol = g.(x_phys, y_phys[j])
        pb_mat[:,j] = solve_column(Dx, gcol, fcol; qr_tol_factor=qr_tol_factor)
    end
    
    # Interpolate along y boundaries
    idx_left = argmin(x_phys)
    idx_right = argmax(x_phys)
    pb_right_vals, pb_left_vals = pb_mat[idx_right, :], pb_mat[idx_left, :]

    interp_pb_right = cheb_barycentric_interpolator(y_phys, pb_right_vals)
    interp_pb_left  = cheb_barycentric_interpolator(y_phys, pb_left_vals)

    # Build functions for 1D integrals
    f_right(y) = interp_pb_right(y)
    g_right(y) = g(b, y)
    f_left(y)  = interp_pb_left(y)
    g_left(y)  = g(a, y)
    
    # Adaptive 1D Levin along y
    I_right = adaptive_levin_1d(f_right, g_right, c, d; tol=tol, k=k, qr_tol_factor=qr_tol_factor)
    I_left  = adaptive_levin_1d(f_left,  g_left,  c, d; tol=tol, k=k, qr_tol_factor=qr_tol_factor)
    
    return - I_right + I_left
end

function levin_2d_separable(f, gx, gy, fdomain::Vector{Float64}, gdomain::Vector{Float64}; k::Int=32, tol::Float64=1e-8, qr_tol_factor::Float64=1.0)
    a, b = fdomain
    c, d = gdomain

    # Outer integrand in x
    function F_outer(x)
        f_inner(y) = f(x, y)
        g_inner(y) = gy(y)
        return adaptive_levin_1d(f_inner, g_inner, c, d; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
    end

    # Outer Levin in x
    return adaptive_levin_1d(F_outer, gx, a, b; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
end

function levin_2d(f, g, fdomain::Vector{Float64}, gdomain::Vector{Float64}; k::Int=32, tol::Float64=1e-8, qr_tol_factor::Float64=1.0)
    # --- Case 1: separable phase g = [gx, gy] ---
    if g isa Vector{Function} && length(g) == 2
        gx, gy = g
        return levin_2d_separable(f, gx, gy, fdomain, gdomain; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
    end

    # --- Case 2: general 2D phase g(x,y) ---
    if g isa Function
        return levin_2d_general(f, g, fdomain, gdomain; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
    end

    error("levin_2d: argument g must be either Function or Vector{Function} of length 2")
end

function adaptive_levin_2d(f, g, fdomain::Vector{Float64}, gdomain::Vector{Float64}; k::Int=32, tol::Float64=1e-8, qr_tol_factor::Float64=1.0, max_subdivisions::Int=10000)
    a0, b0 = fdomain
    c0, d0 = gdomain
    rectangles = [(a0,b0,c0,d0)]
    total = 0.0 + 0.0im
    area = (b0 - a0) * (d0 - c0)

    processed = 0
    while !isempty(rectangles)
        (a,b,c,d) = pop!(rectangles)
        subarea = (b - a) * (d - c)

        processed += 1
        if processed > max_subdivisions
            println("Warning: $max_subdivisions subdivisions is not enough to $tol accuracy. The output may not be accurate.")
            return levin_2d(f,g,[a0,b0],[c0,d0]; k=100, tol=tol, qr_tol_factor=qr_tol_factor)
        end
        # Whole rectangle estimate
        val0 = levin_2d(f,g,[a,b],[c,d]; k=k, tol=tol, qr_tol_factor=qr_tol_factor)

        # Subdivide into quadrants
        mx, my = 0.5*(a+b), 0.5*(c+d)
        R1 = levin_2d(f,g,[a,mx],[c,my]; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
        R2 = levin_2d(f,g,[mx,b],[c,my]; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
        R3 = levin_2d(f,g,[a,mx],[my,d]; k=k, tol=tol, qr_tol_factor=qr_tol_factor)
        R4 = levin_2d(f,g,[mx,b],[my,d]; k=k, tol=tol, qr_tol_factor=qr_tol_factor)

        scale = max(1.0, abs(val0), abs(R1+R2+R3+R4))
        if abs(val0 - (R1+R2+R3+R4)) < tol * scale * sqrt(area / subarea)
            total += R1+R2+R3+R4
        else
            push!(rectangles, (a,mx,c,my))
            push!(rectangles, (mx,b,c,my))
            push!(rectangles, (a,mx,my,d))
            push!(rectangles, (mx,b,my,d))
        end
    end
    return total
end

end