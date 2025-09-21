# AdaptiveLevin.jl

A Julia package for **highly oscillatory integrals** using the **Levin collocation method** with Chebyshev nodes.  
Supports **1D** and **2D oscillatory integrals**, with both fixed and adaptive algorithms.

---

## ✨ Features

- Efficient evaluation of integrals of the form  
  \[
    I = \int f(x) e^{i g(x)} dx
  \]
  and  
  \[
    I = \iint f(x,y) e^{i g(x,y)} dx dy
  \]

- Chebyshev collocation + truncated QR for solving the Levin ODE system.
- Adaptive subdivision to handle varying oscillation frequency.
- Support for **separable phases** \( g(x,y) = g_x(x) + g_y(y) \).
- Clean Julia implementation, no dependencies beyond `LinearAlgebra`.

---

## 📦 Installation

Currently not registered in the General registry.  
You can install from GitHub:

```julia
using Pkg
Pkg.add(url="git@github.com:CuberYyc808/AdaptiveLevin.jl.git")
```
