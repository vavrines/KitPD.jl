"""
sol: temperature field
flux: thermal flux
np: number of particles (except fictitious layers)
nph: particle number in the horizon
idx: particle index in the horizon
bcflag: boundary condition flag
x: particle position
u: displacement
v: velocity
kp: micro-conductivity
ℓ: Bonds' length at initial position
cp: micro-modulus
expan: thermal expansion coefficient
crt: thermal correction
dx: particle interval
isfail: bond-state materix {1: undamaged 0: broken}
ρ: density
cᵥ: specific heat capacity
dt: time step
"""
function update_temperature!(
    sol,
    flux,
    np,
    nph,
    idx,
    bcflag,
    x,
    u,
    v,
    kp,
    ℓ,
    cp,
    expan,
    crt,
    dx,
    isfail,
    ρ,
    cᵥ,
    dt,
)
    @inbounds for i = 1:np
        flux[i, 1] = 0
        for j = 1:nph[i, 1]
            m = idx[i, j]
            maa = begin
                if bcflag[m, 1] == 1
                    0
                else
                    1
                end
            end
            yx = x[m, 1] - x[i, 1] + u[m, 1] - u[i, 1]
            yy = x[m, 2] - x[i, 2] + u[m, 2] - u[i, 2]
            ts = sqrt(yx^2 + yy^2)
            ev = (yx * (v[m, 1] - v[i, 1]) + yy * (v[m, 2] - v[i, 2])) / ts
            ΔT = sol[m, 1] - sol[i, 1]
            flux[i, 1] +=
                (kp * ΔT / ℓ[i, j] - ev * 0.5 * cp * expan * maa) *
                crt[i, j] *
                (dx)^3 *
                isfail[i, j]
        end
    end
    @. sol += flux * dt / (ρ * cᵥ)

    return nothing
end


"""
u: displacement
v: velocity
T: temperature field
pforce: thermal force
bf: body force
np: number of particles (except fictitious layers)
nph: particle number in the horizon
idx: particle index in the horizon
bcflag: boundary condition flag
x: particle position
ℓ: Bonds' length at initial position
cp: micro-modulus
expan: thermal expansion coefficient
crt: mechanical correction
dx: particle interval
isfail: bond-state materix {1: undamaged 0: broken}
mv: 
dmg: damage
ρ: density
dt: time step
"""
function update_mechanics!(
    u,
    v,
    T,
    pforce,
    bf,
    np,
    nph,
    idx,
    bcflag,
    x,
    ℓ,
    cp,
    expan,
    crt,
    dx,
    isfail,
    sc,
    mv,
    dmg,
    ρ,
    dt,
)
    @inbounds for i = 1:np
        pforce[i, :] .= 0.0
        dmg1 = 0.0
        dmg2 = 0.0
        for j = 1:nph[i, 1]
            m = idx[i, j]
            maa = begin
                if bcflag[m, 1] == 1
                    0
                else
                    1
                end
            end
            yx = x[m, 1] - x[i, 1] + u[m, 1] - u[i, 1]
            yy = x[m, 2] - x[i, 2] + u[m, 2] - u[i, 2]
            ts = sqrt(yx^2 + yy^2)
            s = ts / ℓ[i, j] - 1.0
            pforce[i, 1] +=
                (cp * s - (T[m, 1] + T[i, 1]) * 0.5 * cp * expan * maa) * yx / ts *
                crt[i, j] *
                (dx)^3 *
                isfail[i, j]
            pforce[i, 2] +=
                (cp * s - (T[m, 1] + T[i, 1]) * 0.5 * cp * expan * maa) * yy / ts *
                crt[i, j] *
                (dx)^3 *
                isfail[i, j]
            # calculate the crack
            if abs(s - (T[m, 1] + T[i, 1]) * 0.5 * expan) > sc
                isfail[i, j] = 0
                bnd += 1
            end
            dmg1 += isfail[i, j] * mv[i, j]
            dmg2 += mv[i, j]
        end
        dmg[i, 3] = 1 - dmg1 / dmg2
    end
    @. v += (pforce + bf) * dt / ρ
    @. u += v * dt

    return nothing
end

"""
Same as above, except for inputing/returning bnd
"""
function update_mechanics1!(
    u,
    v,
    T,
    pforce,
    bf,
    np,
    nph,
    idx,
    bcflag,
    x,
    ℓ,
    cp,
    expan,
    crt,
    dx,
    isfail,
    sc,
    mv,
    dmg,
    ρ,
    dt,
    bnd,
)
    @inbounds for i = 1:np
        pforce[i, :] .= 0.0
        dmg1 = 0.0
        dmg2 = 0.0
        for j = 1:nph[i, 1]
            m = idx[i, j]
            maa = begin
                if bcflag[m, 1] == 1
                    0
                else
                    1
                end
            end
            yx = x[m, 1] - x[i, 1] + u[m, 1] - u[i, 1]
            yy = x[m, 2] - x[i, 2] + u[m, 2] - u[i, 2]
            ts = sqrt(yx^2 + yy^2)
            s = ts / ℓ[i, j] - 1.0
            pforce[i, 1] +=
                (cp * s - (T[m, 1] + T[i, 1]) * 0.5 * cp * expan * maa) * yx / ts *
                crt[i, j] *
                (dx)^3 *
                isfail[i, j]
            pforce[i, 2] +=
                (cp * s - (T[m, 1] + T[i, 1]) * 0.5 * cp * expan * maa) * yy / ts *
                crt[i, j] *
                (dx)^3 *
                isfail[i, j]
            # calculate the crack
            if abs(s - (T[m, 1] + T[i, 1]) * 0.5 * expan) > sc
                isfail[i, j] = 0
                bnd += 1
            end
            dmg1 += isfail[i, j] * mv[i, j]
            dmg2 += mv[i, j]
        end
        dmg[i, 3] = 1 - dmg1 / dmg2
    end
    @. v += (pforce + bf) * dt / ρ
    @. u += v * dt

    return bnd
end
