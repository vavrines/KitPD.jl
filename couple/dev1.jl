"""
Qs
---
1. mechanical horizon always smaller than thermal case?

"""

using Base.Threads: @threads
using ProgressMeter: @showprogress
using DelimitedFiles

cd(@__DIR__)
include("init.jl")
include("fns.jl")

function main()
    length = 0.1 # meter
    width = 0.1
    nx = 200
    ny = 200
    dx = length / nx
    dy = width / ny
    rh = 3.015 * dx

    emod = 2e11 # modulus Pa
    miu = 0.333 # poison ratio
    kbk = emod / (2 - 2 * miu) # bulk modulus
    ksh = emod / (2 + 2 * miu) # shear modulus
    c = 9 * emod / (π * (rh^3) * dx) # micro-modulus in PD

    kc = 51.9 # thermal conductivity W/mK
    kp = 6 * kc / (π * (rh^3) * dx) # micro-conductivity in PD
    aph = 11.5e-6 # thermal expansion K-1
    cv = 472 # specific heat capacity J/kgK
    dens = 7870 # kg/m^3
    tem0 = 100 # initial temperature

    # boundary fictitious layers 
    lb = 3 # left for temperature
    rb = 3 # right for displacement
    ub = 0 # up
    db = 0 # down
    hm = 50 # maximum number of points in horizon
    tn = (nx + lb + rb) * (ny + ub + db) # total points number (real + ficitious)

    # initialization
    pin = zeros(tn, 2) # initial position
    bf = zeros(tn, 2) # initial body force
    fcs = zeros(tn, 1) # fictitious bc-sign
    fbs = zeros(tn, 1) # bf bc-sign

    isedge = Array{Bool}(undef, tn)
    isedge .= false

    gc = 42320 # critial energy release rate
    sc = sqrt(gc / (rh * (ksh * 6 / pi + 16 / (9 * (pi^2)) * (kbk - 2 * ksh))))

    dt = 1e-9
    nt = 4000#15000 # above ~10000 there will be crack
    t0 = 8000 # critical time step
    u = zeros(tn, 2) # initial displacement
    v = zeros(tn, 2) # initial velocity
    tem = zeros(tn, 1) # initial tempareture changes
    pforce = zeros(tn, 2) # initial inner force
    pflux = zeros(tn, 1) # initial thermal flux
    fail = ones(tn, hm) # initial bond-state materix {1: undamaged 0: broken}
    dmg = zeros(tn, 3) # initial damage array {0:undamaged}

    ccl = 0.02 # pre-existing crack length
    clo = [0, 0] # pre-existing crack location
    bnd = 0

    # discrection
    nnum = 0
    for i = 1:nx
        for j = 1:ny
            cx = -0.5 * (length - dx) + (i - 1) * dx
            cy = -0.5 * (width - dx) + (j - 1) * dx
            nnum += 1
            pin[nnum, 1] = cx
            pin[nnum, 2] = cy
            if cx < dx - 0.5 * length
                fbs[nnum, 1] = 1
            end
            if cx < dx - 0.5 * length ||
            cx > 0.5 * length - dx ||
            cy < dy - 0.5 * width ||
            cy > 0.5 * width - dy
                isedge[nnum] = true
            end
        end
    end
    mtn = nnum # material points number

    # temperature ficitous layers
    for i = 1:lb
        for j = 1:ny
            nnum = nnum + 1
            pin[nnum, 1] = -0.5 * (length + dx) - (i - 1) * dx # left
            pin[nnum, 2] = -0.5 * (width - dx) + (j - 1) * dx
            fcs[nnum, 1] = 1 # thermal
        end
    end

    ## displacement ficitous layers
    for i = 1:rb
        for j = 1:ny
            nnum = nnum + 1
            pin[nnum, 1] = 0.5 * (length + dx) + (i - 1) * dx # right
            pin[nnum, 2] = -0.5 * (width - dx) + (j - 1) * dx
            fcs[nnum, 1] = 2 # mechanical
        end
    end

    hnt, hct, hnm, hcm, ds, fail = init.Horizon(pin, tn, rh, hm, fcs, ccl, clo)
    # hnt: material points number in each horizon for thermo
    # hct: Index matrix of each points in horizon for thermo
    # hnm: material points number in each horizon for mechanic
    # hcm: Index matrix of each points in horizon for mechanic
    # ds: Bonds' length at initial position 
    # fail: pre-existing crack, 1:undamaged 0:broken

    # Correction 
    ## thc: thermal corrcection factors
    ## mec: mechanics corrcection factors 
    thc = init.modifyth(tn, hm, pin, hnt, hct, ds, rh, dx, kp, kc)
    mec, fac = init.modifyme(tn, hm, pin, hnm, hcm, ds, rh, dx, c, emod)
    vmv = fac[1:mtn, :] * dx^3

    #@showprogress for tt = 1:10#nt
    for tt = 1:nt
        println("Step=", tt, ", f-bonds=", bnd)
        # BC
        for i = 1:tn
            if fcs[i, 1] == 1 # left
                tem[i, 1] = 10 * 200 - tem[i-mtn, 1]#2 * 200 - tem[i-mtn, 1] # temperature bc
                #tem[i, 1] = 200 # temperature bc
            elseif fcs[i, 1] == 2 # right
                u[i, :] .= 0 # fixed bc 
            end
            # body force bc
            if fbs[i, 1] == 1
                if tt < t0
                    bf[i, 1] = 1e16 * tt * dt / dx#1e14 * tt * dt / dx
                else
                    bf[i, 1] = 1e16 * t0 * dt / dx#1e14 * t0 * dt / dx
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
            kp,
            ds,
            c,
            aph,
            thc,
            dx,
            fail,
            dens,
            cv,
            dt,
        )

        # mechanical loop
        #=bnd = update_mechanics!(
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
            c,
            aph,
            mec,
            dx,
            fail,
            sc,
            vmv,
            dmg,
            dens,
            dt,
            bnd,
        )=#

        #=update_mechanics1!(
            mtn,
            pforce,
            hnm,
            hcm,
            fcs,
            pin,
            u,
            v,
            ds,
            c,
            aph,
            mec,
            fail,
            tem,
            vmv,
            dmg,
            bf,
            dens,
            dx,
            dt,
            sc,
        )=#

        for i = 1:mtn
            pforce[i, :] .= 0
            dmg1 = 0
            dmg2 = 0
            for j = 1:hnm[i, 1]
                m = hcm[i, j]
                maa = 1
                if fcs[m, 1] == 2
                    maa = 0
                end
                yx = pin[m, 1] - pin[i, 1] + u[m, 1] - u[i, 1]
                yy = pin[m, 2] - pin[i, 2] + u[m, 2] - u[i, 2]
                ts = sqrt(yx^2 + yy^2)
                s = ts / ds[i, j] - 1
                pforce[i, 1] +=
                    (c * s - (tem[m, 1] + tem[i, 1]) * 0.5 * c * aph * maa) * yx / ts *
                    mec[i, j] *
                    (dx)^3 *
                    fail[i, j]
                pforce[i, 2] +=
                    (c * s - (tem[m, 1] + tem[i, 1]) * 0.5 * c * aph * maa) * yy / ts *
                    mec[i, j] *
                    (dx)^3 *
                    fail[i, j]
                ### calculate the crack
                if abs(s - (tem[m, 1] + tem[i, 1]) * 0.5 * aph) > sc ##&& abs(pin[i,1]-clo[1]) < 0.25*length   
                    fail[i, j] = 0
                    bnd += 1
                end
                dmg1 += fail[i, j] * vmv[i, j]
                dmg2 += vmv[i, j]
            end
            dmg[i, 3] = 1 - dmg1 / dmg2
        end
        v += (pforce + bf) * dt / dens     ### update velocity
        u += v * dt


    end

    dmg[:, 1:2] = pin[1:tn, :]
    writedlm("Displacement", u)
    writedlm("Temparature changes", tem)
    writedlm("Damage.txt", dmg)
end

main()
