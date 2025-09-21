# AdaptiveLevin.jl

A Julia package for **highly oscillatory integrals** using the **Levin collocation method** with Chebyshev nodes.  
Supports **1D** and **2D oscillatory integrals**, with both fixed and adaptive algorithms.

---

## Introduction

- Efficient evaluation of integrals of the form
  
  I = \int f(x) e^{i g(x)} dx
  
  and
  
  I = \iint f(x,y) exp\[i g(x,y)\] dx dy

- Chebyshev collocation + truncated QR for solving the Levin ODE system.
- Adaptive subdivision to handle varying oscillation frequency.
- Support for **separable phases** \( g(x,y) = g_x(x) + g_y(y) \).
- Clean Julia implementation, no dependencies beyond `LinearAlgebra`.

---

## Installation

Currently not registered in the General registry.  
You can install from GitHub:

```julia
using Pkg
Pkg.add(url="git@github.com:CuberYyc808/AdaptiveLevin.jl.git")
```
---

## Usage Examples

### 1. One-dimensional oscillatory integral

```julia
using AdaptiveLevin, QuadGK
```
# High-frequency case
```
ω = 1000
f(x) = 1 / (1 + x^2 + 1.23im * x)
g(x) = ω * (x^3 - 2 * x^2 + 5 * x)
integrand(x) = f(x) * exp(1im * g(x))
```
Run Gauss–Kronrod, Levin, and adaptive Levin:
```
@time val_ref = quadgk(integrand, -10.0, 10.0; rtol=1e-8)[1]
@time val_levin = levin_1d(f, g, -10.0, 10.0; k=500)
@time val_adapt = adaptive_levin_1d(f, g, -10.0, 10.0; tol=1e-15, k=12)
```
The three results and time are
- -5.779965334115397e-8 - 3.228457794398132e-11im, 2.284366 seconds
- -5.779965533861069e-8 - 3.229091669112004e-11im, 0.031174 seconds
- -5.779965526562382e-8 - 3.22910511958691e-11im, 0.000978 seconds

For highly oscillatory integrals, Levin method is way better than traditional Gauss–Kronrod method.
If the phase function \(g(x)\) is uneven, adaptive Levin method can significantly improve efficiency compared with that without adaptivity.

# Low-frequency case
ω = 0.001
f_low(x) = 1 / (1 + x^2)
g_low(x) = ω * (x^3 + 2x)
integrand_low(x) = f_low(x) * exp(1im * g_low(x))

println("\n=== Low-frequency (ω=0.001) ===")

@time val_ref_low = quadgk(integrand_low, -10.0, 10.0; rtol=1e-10)[1]
@time val_adapt_low = adaptive_levin_1d(f_low, g_low, -10.0, 10.0; tol=1e-10, k=8)

println("quadgk             = $val_ref_low")
println("adaptive_levin_1d  = $val_adapt_low")



