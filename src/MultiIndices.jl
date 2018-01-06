module MultiIndices

import Base: sub2ind, size, show, getindex

# * Abstract spaces

abstract type Space{N,O} end

# size(mi::Space, i::Integer) = size(mi)[i]

dims(mi::Space{N,O}) where N where O = N

order(mi::Space{N,O}) where N where O = O
function order(𝔄::Space{N𝔄,O𝔄},𝔅::Space{N𝔅,O𝔅}) where N𝔄 where O𝔄 where N𝔅 where O𝔅
    mO𝔄 = maximum(O𝔄)
    (O𝔄..., (o+mO𝔄 for o in O𝔅)...)
end

Base.getindex(mi::Space, args...) = sub2ind(mi, args...)

function sub2ind(mi::Space, args...)
    s = size(mi)
    O = order(mi)
    sub2ind(tuple((s[o] for o in O)...),
            (args[o] for o in O)...)
end

# ** Sum space

struct SumSpace{S<:Space{N,O} where N where O} <: Space{1,1}
    𝔖::Vector{S}
end

dims(mi::SumSpace{S}) where S <: Space{N,O} where N where O = 1
order(mi::SumSpace{S}) where S <: Space{N,O} where N where O = 1

function size(mi::SumSpace)
    sizes = map(mi.𝔖) do s
        prod(size(s))
    end
    sum(sizes)
end

getindex(mi::MultiIndices.SumSpace, i::Integer, args...) = mi.𝔖[i][args...]
getindex(mi::MultiIndices.SumSpace, i::Integer) = mi.𝔖[i]

"""
    𝔄 ⊕ 𝔅

Construct the sum space of which `𝔄` and `𝔅` are terms.
"""
⊕(𝔄::S, 𝔅::S) where S <: Space =
    SumSpace([𝔄, 𝔅])

"""
    𝔖 ⊕ 𝔄

Construct the sum space of which `𝔖` and `𝔄` are terms.
"""
⊕(𝔖::SumSpace{S}, 𝔄::S) where S <: Space =
    SumSpace([𝔖.𝔖..., 𝔄])

# ** Product space

struct ProductSpace{N,O} <: Space{2,(1,2)} where N where O
    𝔄::Space
    𝔅::Space
end
function ProductSpace(𝔄::Space,𝔅::Space)
    N = dims(𝔄) + dims(𝔅)
    O = order(𝔄, 𝔅)
    ProductSpace{N,O}(𝔄, 𝔅)
end

dims(mi::ProductSpace{N,O}) where N where O = N
order(mi::ProductSpace{N,O}) where N where O = O

function size(mi::ProductSpace)
    (size(mi.𝔄)...,size(mi.𝔅)...)
end

"""
    𝔄 ⊗ 𝔅

Construct the product space of which `𝔄` and `𝔅` are factors. `𝔄` will
be the "fast dimension".
"""
⊗(𝔄::Space, 𝔅::Space) = ProductSpace(𝔄, 𝔅)

struct CopySpace{N,O} <: Space{2,(1,2)} where N where O
    𝔖::Space
    n::Integer
end
function CopySpace(𝔖::Space, n::Integer)
    N = dims(𝔖) + 1
    O = order(𝔖)
    O = (O..., length(O)+1)
    CopySpace{N,O}(𝔖, n)
end
dims(mi::CopySpace{N,O}) where N where O = N
order(mi::CopySpace{N,O}) where N where O = O
size(mi::CopySpace) = (size(mi.𝔖)..., mi.n)

# ** Copy space

"""
    𝔖 ⊗ n

Construct the product space from `N` copies of `𝔖`. `𝔖` will
be the "fast dimension".
"""
⊗(𝔖::Space, n::Integer) = CopySpace(𝔖, n)


# * Concrete spaces
# ** Cartesian coordinate systems

struct Cartesian{N,O} <: Space{N,O}
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

# ** Polar coordinates

struct Polar{O} <: Space{2,O}
    nr::Integer
    nϕ::Integer
end
size(mi::Polar) = (mi.nr, mi.nϕ)

# ** Spherical coordinates

struct Spherical{O} <: Space{3,O}
    nr::Integer
    nθ::Integer
    nϕ::Integer
end
size(mi::Spherical) = (mi.nr, mi.nθ, mi.nϕ)

# *** Spherical harmonics

# For expansions in spherical harmonics, given an ℓ, m ∈
# [-ℓ,ℓ]. However, it may be desired to truncate the expansion not
# only in ℓ, but also in m. The "fast dimension" should be m.

struct SphericalHarmonics <: Space{2,(2,1)}
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

struct RotationSpace <: Space{1, (1,)}
    ℓ::Integer
    mmax::Integer
    RotationSpace(ℓ::Integer, mmax::Integer) = new(ℓ, min(ℓ, mmax))
end
RotationSpace(ℓ::Integer) = RotationSpace(ℓ, ℓ)
size(mi::RotationSpace) = 2mi.mmax + 1

ℓs = "spdfghiklmnoqrtuvwxyz"
ℓm(ℓ::Integer, m::Integer) = (ℓ+1, ℓ+m+1)
ℓm(ℓ::Char, m::Integer) = ℓm(search(ℓs, ℓ)-1, m)

function show(io::IO, mi::RotationSpace)
    write(io, "Rot{ℓ = ")
    if mi.ℓ + 1 < length(ℓs)
        write(io, "$(ℓs[mi.ℓ+1])")
    else
        write(io, "$(mi.ℓ)")
    end
    if mi.mmax < mi.ℓ
        write(io, ", mmax=$(mi.mmax)")
    end
    write(io, "}")
end


# This coordinate system is actually proper 3d, but in cylindrically
# symmetric problems (e.g. linear polarization), the inital m quantum
# number is conserved, and the computational effort is 2d.

struct SphericalHarmonics2d <: Space{1, (1,)}
    nℓ::Integer
    m₀::Integer
end
SphericalHarmonics2d(nℓ) = SphericalHarmonics2d(nℓ,0)

# Spin angular momentum; given a spin magnitude s, the projection can
# take values in the interval mₛ ∈ [-s,s].

# ** Spin

struct Spin{S} <: Space{1, (1,)} end
size(mi::Spin{S}) where S = (Integer(2S+1),)

# * Exports

export Space,
    size, getindex, sub2ind, dims, order,
    ⊕, ⊗,
    Cartesian, CartesianXY, CartesianYX,
    CartesianXYZ, CartesianXZY, CartesianYXZ,
    CartesianYZX, CartesianZXY, CartesianZYX,
    Polar, Spherical,
    SphericalHarmonics, SphericalHarmonics2d,
    RotationSpace, ℓm,
    Spin

# * End

end # module
