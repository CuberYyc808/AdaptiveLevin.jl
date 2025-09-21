module ChebyshevUtils

using LinearAlgebra

export cheb_nodes_and_D, cheb_barycentric_interpolator

"""
    cheb_nodes_and_D(k)

Generate Chebyshev-Lobatto nodes and the spectral differentiation matrix.

Arguments
---------
- `k::Int`: number of Chebyshev-Lobatto nodes (must be ≥ 2).

Returns
-------
- `x::Vector{Float64}`: Chebyshev-Lobatto nodes in [-1,1].
- `D::Matrix{Float64}`: Differentiation matrix such that D*v ≈ dv/dx at nodes.
"""

function cheb_nodes_and_D(k::Int)
    if k < 2
        error("k must be ≥ 2")
    end
    N = k - 1
    j = collect(0:N)
    x = cos.(π * j ./ N)     # Chebyshev-Lobatto nodes
    c = ones(Float64, N+1)
    c[1] = 2.0; c[end] = 2.0
    c = c .* ((-1) .^ j)

    X = repeat(x', N+1, 1)
    dX = X .- X'
    D = (c * (1.0 ./ c)') ./ (dX .+ Matrix(I, N+1, N+1))   # off-diagonal
    D = D .- diagm(0 => sum(D, dims=2)[:])                 # diagonal correction
    return x, D
end

"""
    cheb_barycentric_interpolator(x_nodes, vvals)

Build a barycentric Lagrange interpolator over Chebyshev-Lobatto nodes.

Arguments
---------
- `x_nodes::Vector{Float64}`: Chebyshev-Lobatto nodes.
- `vvals::Vector{ComplexF64}`: function values at nodes.

Returns
-------
- `eval_fn(y::Float64)`: closure that interpolates v(y).
"""

function cheb_barycentric_interpolator(x_nodes::Vector{}, vvals::Vector{})
    k = length(x_nodes)
    w = [ ((j==1 || j==k) ? 0.5 : 1.0) * (-1.0)^(j-1) for j in 1:k ]  # barycentric weights

    function eval_at(y)
        # Handle exact node matches
        for (jj, xv) in enumerate(x_nodes)
            if isapprox(y, xv; atol=0, rtol=0)
                return vvals[jj]
            end
        end
        num = 0.0
        den = 0.0
        for j in 1:k
            diff = y - x_nodes[j]
            bj = w[j] / diff
            num += bj * vvals[j]
            den += bj
        end
        return num / den
    end
    return eval_at
end

end 