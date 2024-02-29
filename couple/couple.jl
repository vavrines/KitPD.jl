"""
Qs
---
Universal constants
k = 1.38065e-23

Argon gas
m = 66.3e-27 kg
R = 208.2428355957768

Reference state
L₀ = 0.1 m
T₀ = 273 K
U₀ = 337.1951782503631 m/s
t₀ = 0.0002965641457831028 s
"""

using Base.Threads: @threads
using ProgressMeter: @showprogress
using DelimitedFiles
using KitBase

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

    xyidx = Array{CartesianIndex}(undef, 1840)
    #tem .*= tem0
end

begin
    dt = 1e-7
    nt = 400#15000 # above ~10000 there will be crack
    t0 = 8000 # critical time step
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

        xyidx[nnum] = CartesianIndex(i, j)
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

        xyidx[nnum] = CartesianIndex(-i + 1, j)
    end
end

## displacement ficitous layers
for i = 1:rb
    for j = 1:nx
        nnum = nnum + 1
        pin[nnum, 1] = 0.5 * (len + dx) + (i - 1) * dx # right
        pin[nnum, 2] = -0.5 * (len - dx) + (j - 1) * dx
        fcs[nnum, 1] = 2 # mechanical

        xyidx[nnum] = CartesianIndex(nx + i, j)
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

#@showprogress for tt = 1:10#nt
for tt = 1:nt
    println("Step=", tt, ", f-bonds=", bnd)
    # BC
    for i = 1:tn
        if fcs[i, 1] == 1 # left
            tem[i, 1] = 273 - tem[i-mtn, 1]#2 * 200 - tem[i-mtn, 1] # temperature bc
        #tem[i, 1] = 273 # temperature bc
        elseif fcs[i, 1] == 2 # right
            u[i, :] .= 0 # fixed bc 
        end
        # body force bc
        if fbs[i, 1] == 1
            if tt < t0
                #bf[i, 1] = 1e14 * tt * dt / dx#1e14 * tt * dt / dx
            else
                #bf[i, 1] = 1e14 * t0 * dt / dx#1e14 * t0 * dt / dx
            end
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
end

dmg[:, 1:2] = pin[1:tn, :]
writedlm("Displacement", u)
writedlm("Temparature changes", tem)
writedlm("Damage.txt", dmg)

#--- fluid ---#
#=set = Setup(
    case = "square",
    space = "2d0f0v",
    boundary = ["fix", "extra", "mirror", "extra"],
    limiter = "vanleer",
    cfl = 0.5,
    maxTime = 2.0,
    flux = "hll",
    hasForce = false,
)
ps = PSpace2D(-0.15, 0.05, nx*2, -0.05, 0.15, nx*2, 1, 1)
vs = nothing
gas = Gas(Kn = 1e-3, Ma = 0.9, K = 1.0)

prim0 = [1.0, 0.0, 0.0, 1.0]
prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
fw = function(x, y, args...)
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

    if ks.ps.y[idx] < -ks.ps.y0 - ks.ps.dx[1]/2
        xbis[iter][1] = ks.ps.x[idx] - ks.ps.dx[idx]/2
        xbis[iter][2] = ks.ps.y[idx]
        nbis[iter] .= [-1.0, 0.0]
    else
        xbis[iter][1] = ks.ps.x[idx]
        xbis[iter][2] = ks.ps.y[idx] + ks.ps.dy[idx]/2
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
=#
