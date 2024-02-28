using Base: @kwdef

@kwdef mutable struct PDMater{T} <: KitBase.AbstractProperty
    δ::T
    rh::T = 3.015 * δ
    emod::T = 2e11 # modulus Pa
    ν::T = 0.333 # poison ratio
    kbk::T = emod / (2 - 2 * ν) # bulk modulus
    ksh::T = emod / (2 + 2 * ν) # shear modulus
    c::T = 9.0 * emod / (π * (rh^3) * δ) # micro-modulus in PD
    kc::T = 51.9 # thermal conductivity W/mK
    kp::T = 6.0 * kc / (π * (rh^3) * dx) # micro-conductivity in PD
    aph::T = 11.5e-6 # thermal expansion K-1
    cᵥ::T = 472.0 # specific heat capacity J/kgK
    ρ::T = 7870.0 # kg/m^3
    gc::T = 42320.0 # critial energy release rate
    sc::T = sqrt(gc / (rh * (ksh * 6.0 / π + 16.0 / (9.0 * (π^2)) * (kbk - 2.0 * ksh))))
end

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


function update_ghost!(ctr, ps, gas, ib)
    ghost_ids, xbis, nbis, xips, ip_cids, ip_nids, ip_bids =
        ib.idg, ib.xb, ib.nb, ib.xi, ib.idic, ib.idin, ib.idib

    for iter in eachindex(ip_nids)
        nid = ip_nids[iter]
        bid = ip_bids[iter]
        xf, yf = xips[iter]
        if length(nid) == 4
            pos = [1, xf, yf, xf * yf]

            # U
            w1 = [ctr[idx].prim[2] for idx in nid]
            w2 = zeros(length(bid))
            w = [w1; w2]

            C = KB.bilinear_coeffs(ps, xbis, nbis, nid, bid, w)
            U1 = C' * pos

            # V
            w1 = [ctr[idx].prim[3] for idx in nid]
            w2 = zeros(length(bid))
            w = [w1; w2]

            C = KB.bilinear_coeffs(ps, xbis, nbis, nid, bid, w)
            V1 = C' * pos

            # T
            w1 = [1 / ctr[idx].prim[4] for idx in nid]
            w2 = ones(length(bid))
            w = [w1; w2]

            C = KB.bilinear_coeffs(ps, xbis, nbis, nid, bid, w)
            T1 = C' * pos

            # p
            w1 = [0.5 * ctr[idx].prim[1] / ctr[idx].prim[4] for idx in nid]
            w = w1

            C = KB.bilinear_coeffs(ps, xbis, nbis, nid, bid, w)
            P1 = C' * pos

            ρ1 = 2 * P1 / T1
        else
            idx1 = ip_cids[iter]
            ρ1, U1, V1, λ1 = ctr[idx1].prim
            T1 = 1.0 / λ1
            P1 = 0.5 * ρ1 / λ1
        end
        
        T0 = 2 - T1 # here the temperature of solid wall is set as 2.
                    # we only need to update here when PD codes get involved.
        ρ0 = ρ1 * T1 / T0

        idx = ghost_ids[iter]
        ctr[idx].prim .= [ρ0, -U1, -V1, 1 / T0]
        ctr[idx].w .= prim_conserve(ctr[idx].prim, gas.γ)
    end
end

function update_field!(KS, ctr, a1face, a2face, flags, residual)
    nx, ny, dx, dy = KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for j ∈ 1:ny
        for i ∈ 1:nx
            if flags[i, j] == 1
                KB.step!(
                    ctr[i, j].w,
                    ctr[i, j].prim,
                    a1face[i, j].fw,
                    a1face[i+1, j].fw,
                    a2face[i, j].fw,
                    a2face[i, j+1].fw,
                    KS.gas.γ,
                    dx[i, j] * dy[i, j],
                    sumRes,
                    sumAvg,
                    :bgk,
                )
            end
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    KitBase.bc_extra!(ctr; dirc = :xr)
    KitBase.bc_extra!(ctr; dirc = :yr)
    KitBase.bc_mirror!(ctr; dirc = :yl)

    return nothing
end
