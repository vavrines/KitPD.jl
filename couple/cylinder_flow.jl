using KitBase, Plots
using KitBase.JLD2
using Base.Threads: @threads
using KitBase.ProgressMeter: @showprogress

cd(@__DIR__)

set = Setup(
    case = "square",
    space = "2d0f0v",
    boundary = ["fix", "extra", "mirror", "extra"],
    limiter = "vanleer",
    cfl = 0.3,
    maxTime = 2.0,
    flux = "hll",
    hasForce = false,
)
#ps = PSpace2D(-0.15, 0.05, 400, -0.05, 0.15, 400, 1, 1) # Hu's geometry
ps = PSpace2D(-0.15, 0.05, 40, -0.05, 0.15, 40, 1, 1) # coarse geometry
vs = nothing
gas = Gas(Kn = 1e-3, Ma = 0.9, K = 1.0)

prim0 = [1.0, 0.0, 0.0, 1.0]
prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
fw = function (x, y, args...)
    pr = bc(x, y, args...)
    prim_conserve(pr, gas.γ)
end
bc = function (x, y, args...)
    if -ps.x1 < x < ps.x1 && ps.y0 < y < -ps.y0
        return prim0
    else
        return prim1
    end
end
ib0 = IB(fw, bc, NamedTuple())

ks = SolverSet(set, ps, vs, gas, ib0)
ctr, a1face, a2face = init_fvm(ks; structarray = true)

flags = ones(Int, axes(ps.x))
for i in axes(flags, 1), j in axes(flags, 2)
    if -ps.x1 < ps.x[i, j] < ps.x1 && ps.y0 < ps.y[i, j] < -ps.y0
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

    if ks.ps.y[idx] < -ks.ps.y0 - ks.ps.dx[1] / 2
        xbis[iter][1] = ks.ps.x[idx] - ks.ps.dx[idx] / 2
        xbis[iter][2] = ks.ps.y[idx]
        nbis[iter] .= [-1.0, 0.0]
    else
        xbis[iter][1] = ks.ps.x[idx]
        xbis[iter][2] = ks.ps.y[idx] + ks.ps.dy[idx] / 2
        nbis[iter] .= [0.0, 1.0]
    end
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

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)

@showprogress for iter = 1:200
    evolve!(ks, ctr, a1face, a2face, dt)
    update_ghost!(ctr, ks.ps, ks.gas, ib)
    update_field!(ks, ctr, a1face, a2face, flags, res)

    global t += dt
end

plot(ks, ctr)
