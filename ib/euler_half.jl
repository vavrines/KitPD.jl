using KitBase, Plots
using KitBase.JLD2
using Base.Threads: @threads
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)

set = Setup(
    case = "cylinder",
    space = "2d0f0v",
    boundary = ["fix", "extra", "mirror", "extra"],
    limiter = "vanleer",
    cfl = 0.3,
    maxTime = 2.0, # time
    flux = "hllc",
    hasForce = false,
)
ps = PSpace2D(0, 6, 120, 0, 2, 40, 1, 1)
gas = Gas(Kn = 1e-3, Ma = 0.8, K = 1.0)

prim0 = [1.0, 0.0, 0.0, 1.0]
prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
fw = (args...) -> prim_conserve(prim1, gas.γ)
bc = function (x, y, args...)
    if abs((x - 3)^2 + y^2 - 1) < 1e-3
        return prim0
    else
        return prim1
    end
end
ib0 = IB(fw, bc, NamedTuple())

ks = SolverSet(set, ps, nothing, gas, ib0)
ctr, a1face, a2face = init_fvm(ks; structarray = true)

radius = 1
flags = ones(Int, axes(ps.x))
for i in axes(flags, 1), j in axes(flags, 2)
    if (ps.x[i, j] - 3)^2 + ps.y[i, j]^2 < radius # (x-3)^2+y^2=1
        flags[i, j] = 0
    end
end
flags[0, :] .= -1
flags[ks.ps.nx+1, :] .= -1
flags[:, 0] .= -1
flags[:, ks.ps.ny+1] .= -1
KB.ghost_flag!(ps, flags)

ghost_ids = findall(flags .== -2)
xbis = [Vector{Float64}(undef, 2) for iter = 1:length(ghost_ids)]
nbis = zero.(xbis)
for iter in axes(xbis, 1)
    idx = ghost_ids[iter]

    θ = atan(ks.ps.y[idx], (ks.ps.x[idx] - 3))
    xbis[iter][1] = 3 + radius * cos(θ)
    xbis[iter][2] = radius * sin(θ)

    nbis[iter] .= [radius * cos(θ), radius * sin(θ)]
end

xips = KB.ip_location(ps, ghost_ids, xbis)
ip_cids, ip_nids, ip_bids = KB.ip_connectivity(ps, xips, flags)
ib = SharpIB{
    typeof(flags),
    typeof(ghost_ids),
    typeof(xbis),
    typeof(ip_cids),
    typeof(ip_nids),
    typeof(ip_bids),
}(
    flags,
    ghost_ids,
    xbis,
    nbis,
    xips,
    ip_cids,
    ip_nids,
    ip_bids,
)

wbis = [zeros(4) for i = 1:size(ib.xb, 1)]
primbis = zero.(wbis)
solb = Solution2D(wbis, primbis)
soli = Solution2D(deepcopy(wbis), deepcopy(primbis))

function update_ghost!(ctr, ps, gas, ib)
    ghost_ids, xbis, nbis, xips, ip_nids, ip_bids =
        ib.idg, ib.xb, ib.nb, ib.xi, ib.idin, ib.idib

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
        else
            idx = ip_cids[iter]
            
            U1 = ctr[idx].prim[2]
            V1 = ctr[idx].prim[3]
            T1 = 1 / ctr[idx].prim[4]
            P1 = 0.5 * ctr[idx].prim[1] / ctr[idx].prim[4]
        end

        ρ1 = 2 * P1 / T1
        T0 = 2 - T1
        ρ0 = ρ1 * T1 / T0

        idx = ghost_ids[iter]
        ctr[idx].prim .= [ρ0, -U1, -V1, 1 / T0]
        ctr[idx].w .= prim_conserve(ctr[idx].prim, gas.γ)
    end

    return nothing
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

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)

@showprogress for iter = 1:nt
    evolve!(ks, ctr, a1face, a2face, dt)
    update_ghost!(ctr, ks.ps, ks.gas, ib)
    update_field!(ks, ctr, a1face, a2face, flags, res)

    global t += dt
end

plot(ks, ctr)
