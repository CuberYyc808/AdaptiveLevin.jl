# AdaptiveLevin.jl

A Julia package for **highly oscillatory integrals** using the **Levin collocation method** with Chebyshev nodes.  
Supports **1D** and **2D oscillatory integrals**, with both fixed and adaptive algorithms.

---

## Introduction

- Efficient evaluation of integrals of the form
  
  $$I = \int f(x) \exp\[i g(x)\] dx$$
  
  and
  
  $$I = \iint f(x,y) \exp\[i g(x,y)\] dx dy$$

- Chebyshev collocation + truncated QR for solving the Levin ODE system.
- Adaptive subdivision to handle varying oscillation frequency.
- Support for **separable phases** $` g(x,y) = g_x(x) + g_y(y) `$.
- Clean Julia implementation, no dependencies beyond `LinearAlgebra`.

---

## Installation

You can install it directly by 

```julia
using Pkg
Pkg.add("AdaptiveLevin")
```
---

## Usage Examples

### 1. One-dimensional oscillatory integral

```julia
using AdaptiveLevin, QuadGK

ω = 10000
f(x) = 1 / (1 + x^2)
g(x) = ω * (x^2 + 2 * x)
integrand(x) = f(x) * exp(1im * g(x))
```
Run Gauss–Kronrod and adaptive Levin:
```julia
@time val_ref = quadgk(integrand, -10.0, 10.0; rtol=1e-8)[1]
@time val_adapt = adaptive_levin_1d(f, g, -10.0, 10.0; tol=1e-15, k=24)
```
The results and time are
```julia
-0.00788183605964596 - 0.004051885499350059im, 1.835586 seconds # Gauss–Kronrod
-0.00788183605964634 - 0.004051885499358002im, 0.002167 seconds # adaptive-Levin
```
- For highly oscillatory integrals, Levin method is way better than traditional Gauss–Kronrod method.
- If the phase function $`g(x)`$ is uneven, adaptive Levin method can significantly improve efficiency compared with that without adaptivity.
- Adaptivity can effectively avoid the pathological phenomena at stable points. For example, $`g'(-2)=0`$ in the above case.
- In low frequency cases, the adaptive Levin method is till fast and accurate. (Try ω=0.001 to see its performance.) 

### 2. Two-dimensional oscillatory integral
```julia
using AdaptiveLevin, Integrals

ω = 100
fxy(x,y) = 1.0 / (1.0 + x^2 + y^2)
gxy(x,y) = ω * (x^2 + y^3)
integrand(x, y) = fxy(x, y) * exp(1im * gxy(x, y))
```
Reference with HCubature:
```julia
f(u, p) = integrand(u[1], u[2])
domain = ([-1.0, -1.0], [1.0, 1.0])
prob = IntegralProblem(f, domain)
@time val_ref = solve(prob, HCubatureJL(); reltol = 1e-8).u
```
2D adaptive Levin method:
```julia
@time val_levin = adaptive_levin_2d(fxy, gxy, [-1.0, 1.0], [-1.0, 1.0]; tol=1e-8, k=16)
```
The results and time are
```julia
0.04087258215148527 + 0.03989088975628127im, 9.716638 seconds # HCubature
0.04087258215149002 + 0.03989088975628428im, 0.142308 seconds # adaptive-Levin
```
If the phase function is separable, namely $`g(x,y) = g_x(x) + g_y(y)`$, we also accept input as `g = [gx, gy]`:
```julia
gx(x) = ω * x^2
gy(y) = ω * y^3
@time val_sep = adaptive_levin_2d(fxy, [gx, gy], [-1.0, 1.0], [-1.0, 1.0]; tol=1e-8, k=12)
```
It gives
```julia
0.04087258215150031 + 0.03989088975690701im, 0.380358 seconds # adaptive-Levin separable phase
```


