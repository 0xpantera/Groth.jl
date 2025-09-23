"""
ProjectivePoint{Curve,F}

Generic Jacobian projective point storage shared across curve families. It holds
(X, Y, Z) coordinates in homogeneous form; concrete curves provide doubling and
addition formulas.
"""
struct ProjectivePoint{Curve,F} <: GroupElem{Curve}
    coords::SVector{3,F}

    function ProjectivePoint{Curve,F}(coords::SVector{3,F}) where {Curve,F}
        new{Curve,F}(coords)
    end

    function ProjectivePoint{Curve,F}(X::F, Y::F, Z::F) where {Curve,F}
        new{Curve,F}(SVector(X, Y, Z))
    end
end

# Coordinate helpers
Base.getindex(p::ProjectivePoint{Curve,F}, i::Int) where {Curve,F} = p.coords[i]

x_coord(p::ProjectivePoint) = p[1]
y_coord(p::ProjectivePoint) = p[2]
z_coord(p::ProjectivePoint) = p[3]

# Zero / infinity handling
Base.zero(::Type{ProjectivePoint{Curve,F}}) where {Curve,F} =
    ProjectivePoint{Curve,F}(one(F), one(F), zero(F))
Base.zero(p::ProjectivePoint{Curve,F}) where {Curve,F} = zero(ProjectivePoint{Curve,F})

Base.iszero(p::ProjectivePoint) = iszero(z_coord(p))

# Equality (converted to affine coordinates when needed)
function Base.:(==)(p::ProjectivePoint{Curve,F}, q::ProjectivePoint{Curve,F}) where {Curve,F}
    if iszero(p) && iszero(q)
        return true
    elseif iszero(p) || iszero(q)
        return false
    end
    xp, yp = to_affine(p)
    xq, yq = to_affine(q)
    return xp == xq && yp == yq
end

# to_affine is curve-specific because it depends on field operations; dispatch here
function to_affine end

# Utility to build from integers when field constructor is provided
function projective_point_from_integers(::Type{ProjectivePoint{Curve,F}}, fx, fy, fz) where {Curve,F}
    return ProjectivePoint{Curve,F}(fx, fy, fz)
end

export ProjectivePoint
