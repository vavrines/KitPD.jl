using KitBase, KitPD
using DelimitedFiles, Plots
using Base.Threads: @threads
using ProgressMeter: @showprogress

cd(@__DIR__)
include("init.jl")
include("fns.jl")

#--- material ---#
begin
    len = 0.1 # meter
    nx = 40
    dx = len / nx
    mat = PDMater(δ = dx)
    tem0 = 273 # initial temperature
end

begin # boundary fictitious layers
    lb = 3 # left for temperature
    rb = 3 # right for displacement
    ub = 0 # up
    db = 0 # down
    hm = 50 # maximum number of points in horizon
    tn = (nx + lb + rb) * (nx + ub + db) # total points number (real + ficitious)
end

begin # initialization
    pin = zeros(tn, 2) # initial position
    bf = zeros(tn, 2) # initial body force
    fcs = zeros(tn, 1) # fictitious bc-sign
    fbs = zeros(tn, 1) # bf bc-sign
    u = zeros(tn, 2) # initial displacement
    v = zeros(tn, 2) # initial velocity
    tem = zeros(Float64, tn, 1) # initial tempareture changes
    pforce = zeros(tn, 2) # initial inner force
    pflux = zeros(tn, 1) # initial thermal flux
    fail = ones(tn, hm) # initial bond-state materix {1: undamaged 0: broken}
    dmg = zeros(tn, 3) # initial damage array {0:undamaged}

    isedge = Array{Bool}(undef, tn)
    isedge .= false

    ixy_index = Array{CartesianIndex}(undef, 1840)
    #tem .*= tem0
end

begin
    ccl = 0.02 # pre-existing crack length
    clo = [0, 0] # pre-existing crack location
    bnd = 0
end

# discrection
nnum = 0
for i = 1:nx
    for j = 1:nx
        cx = -0.5 * (len - dx) + (i - 1) * dx
        cy = -0.5 * (len - dx) + (j - 1) * dx
        nnum += 1
        pin[nnum, 1] = cx
        pin[nnum, 2] = cy
        if cx < dx - 0.5 * len
            fbs[nnum, 1] = 1
        end

        if cx < dx - 0.5 * len ||
           cx > 0.5 * len - dx ||
           cy < dx - 0.5 * len ||
           cy > 0.5 * len - dx
            isedge[nnum] = true
        end

        ixy_index[nnum] = CartesianIndex(i, j)
    end
end
mtn = nnum # material points number

# temperature ficitous layers
for i = 1:lb
    for j = 1:nx
        nnum = nnum + 1
        pin[nnum, 1] = -0.5 * (len + dx) - (i - 1) * dx # left
        pin[nnum, 2] = -0.5 * (len - dx) + (j - 1) * dx
        fcs[nnum, 1] = 1 # thermal

        ixy_index[nnum] = CartesianIndex(-i + 1, j)
    end
end

## displacement ficitous layers
for i = 1:rb
    for j = 1:nx
        nnum = nnum + 1
        pin[nnum, 1] = 0.5 * (len + dx) + (i - 1) * dx # right
        pin[nnum, 2] = -0.5 * (len - dx) + (j - 1) * dx
        fcs[nnum, 1] = 2 # mechanical

        ixy_index[nnum] = CartesianIndex(nx + i, j)
    end
end

hnt, hct, hnm, hcm, ds, fail = init.Horizon(pin, tn, mat.rh, hm, fcs, ccl, clo)
# hnt: material points number in each horizon for thermo
# hct: Index matrix of each points in horizon for thermo
# hnm: material points number in each horizon for mechanic
# hcm: Index matrix of each points in horizon for mechanic
# ds: Bonds' length at initial position 
# fail: pre-existing crack, 1:undamaged 0:broken

# Correction 
## thc: thermal corrcection factors
## mec: mechanics corrcection factors 
thc = init.modifyth(tn, hm, pin, hnt, hct, ds, mat.rh, dx, mat.kp, mat.kc)
mec, fac = init.modifyme(tn, hm, pin, hnm, hcm, ds, mat.rh, dx, mat.c, mat.emod)
vmv = fac[1:mtn, :] * dx^3

#--- fluid ---#
set = Setup(
    case = "square",
    space = "2d0f0v",
    boundary = ["fix", "extra", "mirror", "extra"],
    limiter = "vanleer",
    cfl = 0.2,
    maxTime = 2.0,
    flux = "hllc",
    hasForce = false,
)
ps = PSpace2D(-0.15, 0.05, nx * 2, -0.05, 0.15, nx * 2, 1, 1)
vs = nothing
gas = Gas(Kn = 1e-3, Ma = 1.2, K = 1.0)

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

dtf = timestep(ks, ctr, 0.0)
dt = dtf * 0.002965641457831028
nt = 200
res = zeros(4)

ibpd_index = Array{Int}(undef, length(ib.idic))
for i in eachindex(ibpd_index)
    idx = ib.idg[i][1] - ks.ps.nx ÷ 2
    idy = ib.idg[i][2]
    id = findall(x -> x == CartesianIndex(idx, idy), ixy_index)
    @assert length(id) == 1
    ibpd_index[i] = id[1]
end

#--- main loop ---#
#@showprogress for tt = 1:10#nt
for tt = 1:100#nt
    println("Step=", tt, ", f-bonds=", bnd)
    # BC
    for i = mtn+1:tn
        if fcs[i, 1] == 1 # left
            ix = ks.ps.nx ÷ 2
            iy = ixy_index[i][2]
            Tf = 1.0 / ctr[ix, iy].prim[end] * tem0
            tem[i, 1] = Tf - tem[i-mtn, 1]
        elseif fcs[i, 1] == 2 # right
            u[i, :] .= 0
        end
    end

    # thermal loop
    update_temperature!(
        tem,
        pflux,
        mtn,
        hnt,
        hct,
        fcs,
        pin,
        u,
        v,
        mat.kp,
        ds,
        mat.c,
        mat.aph,
        thc,
        dx,
        fail,
        mat.ρ,
        mat.cᵥ,
        dt,
    )

    # mechanical loop
    bnd = update_mechanics1!(
        u,
        v,
        tem,
        pforce,
        bf,
        mtn,
        hnm,
        hcm,
        fcs,
        pin,
        ds,
        mat.c,
        mat.aph,
        mec,
        dx,
        fail,
        mat.sc,
        vmv,
        dmg,
        mat.ρ,
        dt,
        bnd,
    )

    evolve!(ks, ctr, a1face, a2face, dtf)

    #=@threads for j = 1:ks.ps.ny
        for i = 1:ks.ps.nx+1
            KB.flux_hllc!(a1face[i, j].fw, ctr[i-1, j].w, ctr[i, j].w, 5/3, dtf, dx)
        end
    end
    @threads for j = 1:ks.ps.ny+1
        for i = 1:ks.ps.nx
            wL = local_frame(ctr[i, j-1].w, [0, 1])
            wR = local_frame(ctr[i, j].w, [0, 1])
            KB.flux_hllc!(a2face[i, j].fw, wL, wR, 5/3, dtf, dx)
            a2face[i, j].fw .= global_frame(a2face[i, j].fw, [0, 1])
        end
    end=#

    update_ghost2!(ctr, ks.ps, ks.gas, ib, ibpd_index, tem)
    update_field!(ks, ctr, a1face, a2face, flags, res)
end

# plot
plot(ks, ctr)

sol = extract_sol(ks, ctr)
plot(ks.ps.x[1:40, 1], sol[1:40, 20, 2])

# write
dmg[:, 1:2] = pin[1:tn, :]
writedlm("Displacement", u)
writedlm("Temparature changes", tem)
writedlm("Damage.txt", dmg)
