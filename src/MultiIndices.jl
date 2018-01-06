module MultiIndices

import Base: sub2ind, size

abstract type MultiIndex{N,O} end

size(mi::MultiIndex, i::Integer) = size(mi)[i]

Base.getindex(mi::MultiIndex, args...) = sub2ind(mi, args...)

function sub2ind(mi::MultiIndex{N,O}, args...) where N where O
    sub2ind(tuple((size(mi, o) for o in O)...),
            (args[o] for o in O)...)
end

struct ProductSpace <: MultiIndex{2,(2,1)}

end

"""
    𝔄 ⊗ 𝔅

Construct the product space of which `𝔄` and `𝔅` are factors. `𝔅` will
be the "fast dimension".
"""
function ⊗(𝔄::MultiIndex, 𝔅::MultiIndex)
    println("Product space of $(𝔄) and $(𝔅)")
    𝔄
end

"""
    N ⊗ 𝔖

Construct the product space `N` copies of `𝔖`. `𝔖` will
be the "fast dimension".
"""
function ⊗(N::Integer, 𝔖::MultiIndex)
    𝔖
end

# Cartesian coordinate systems

struct Cartesian{N,O} <: MultiIndex{N,O}
    dims::Tuple
end

size(mi::Cartesian) = mi.dims

Cartesian(i::Integer) = Cartesian{1,(1,)}((i,))
Cartesian{N,O}(args...) where N where O = Cartesian{N,O}(args)

CartesianXY = Cartesian{2, (1,2)}
CartesianYX = Cartesian{2, (2,1)}

CartesianXYZ = Cartesian{3, (1,2,3)}
CartesianXZY = Cartesian{3, (1,3,2)}

CartesianYXZ = Cartesian{3, (2,1,3)}
CartesianYZX = Cartesian{3, (2,3,1)}

CartesianZXY = Cartesian{3, (3,1,2)}
CartesianZYX = Cartesian{3, (3,2,1)}

# Polar coordinates

struct Polar{O} <: MultiIndex{2,O}
    nr::Integer
    nϕ::Integer
end
size(mi::Polar) = (mi.nr, mi.nϕ)

# Spherical coordinates

struct Spherical{O} <: MultiIndex{3,O}
    nr::Integer
    nθ::Integer
    nϕ::Integer
end
size(mi::Spherical) = (mi.nr, mi.nθ, mi.nϕ)

# For expansions in spherical harmonics, given an ℓ, m ∈
# [-ℓ,ℓ]. However, it may be desired to truncate the expansion not
# only in ℓ, but also in m. The "fast dimension" should be m.

struct SphericalHarmonics <: MultiIndex{2,(2,1)}
    nℓ::Integer
    nm::Integer
end
SphericalHarmonics(nℓ) = SphericalHarmonics(nℓ, 2nℓ-1)

# It would be nice if ℓ was zero-based, and m from -ℓ–ℓ. Maybe
# "custom" indexing can be implemented like this?
function (mi::SphericalHarmonics)(ℓ::Integer, m::Integer)
    ℓ ∉ 0:mi.nℓ-1 && error("Invalid value of ℓ = $(ℓ) ∉ [0,$(mi.nℓ-1)]")
    abs(m) > ℓ && error("Invalid value of m = $(m) ∉ [-$(ℓ),$(ℓ)]")
    2abs(m) > mi.nm && error("Invalid value of |m| = $(abs(m)) > $((mi.nm-1)/1) (max m)")

end

# This function actually returns the *maximum* nm, that is nm for the
# largest ℓ. The expansion is not "rectangular", but rather "triangular":
#=

              ℓ

    0   1    2    3    4   ...

                      -4
                 -3   -3
            -2   -2   -2
       -1   -1   -1   -1
m   0   0    0    0    0
        1    1    1    1
             2    2    2
                  3    3
                       4
=#
size(mi::SphericalHarmonics) = (mi.nℓ, mi.nm)

# This coordinate system is actually proper 3d, but in cylindrically
# symmetric problems (e.g. linear polarization), the inital m quantum
# number is conserved, and the computational effort is 2d.

struct SphericalHarmonics2d <: MultiIndex{1, (1,)}
    nℓ::Integer
    m₀::Integer
end
SphericalHarmonics2d(nℓ) = SphericalHarmonics2d(nℓ,0)

# Spin angular momentum; given a spin magnitude s, the projection can
# take values in the interval mₛ ∈ [-s,s].

struct Spin{S} <: MultiIndex{1, (1,)} end
size(mi::Spin{S}) where S = (Integer(2S+1),)

# Exports

export MultiIndex,
    size, getindex, sub2ind, ⊗,
    Cartesian, CartesianXY, CartesianYX,
    CartesianXYZ, CartesianXZY, CartesianYXZ,
    CartesianYZX, CartesianZXY, CartesianZYX,
    Polar, Spherical,
    SphericalHarmonics, SphericalHarmonics2d,
    Spin

end # module
