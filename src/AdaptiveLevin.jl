module AdaptiveLevin

include("AdaptiveLevin1D.jl")
include("AdaptiveLevin2D.jl")
using .AdaptiveLevin1D
using .AdaptiveLevin2D

export adaptive_levin_1d, adaptive_levin_2d, levin_1d, levin_2d

end