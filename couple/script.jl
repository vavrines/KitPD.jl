### Main Fully Coupled Thermomechanics
### Bonds-based PD; Transient

cd(@__DIR__)
include("init.jl")

using DelimitedFiles

function main()
    length = 0.1 ## meter
    width = 0.1  ## meter
    nx = 200
    ny = 200
    dx = length / nx
    rh = 3.015 * dx
    emod = 2e11                ## module Pa
    miu = 0.333                ## poison ritio
    kbk = emod / (2 - 2 * miu)       ## Bulk module
    ksh = emod / (2 + 2 * miu)       ## Shear module
    c = 9 * emod / (π * (rh^3) * dx)  ## micromodule in PD

    kc = 51.9                  ## thermal conductivity W/mK
    kp = 6 * kc / (π * (rh^3) * dx)    ## microconductivity in PD
    aph = 11.5e-6              ## thermal expansion K-1
    cv = 472                   ## specific heat capacity J/kgK
    dens = 7870                ## kg/m^3
    tem0 = 100                 ## initial temperature

    #bc ficitious layers 
    lb = 3                     ## lift for temperature
    rb = 3                     ## right for displacement
    ub = 0                     ## up
    db = 0                     ## down
    hm = 50                    ## maximum number of Points in horizon
    tn = (nx + lb + rb) * (ny + ub + db)

    ## initialization
    pin = zeros(tn, 2)          ## initial position
    bf = zeros(tn, 2)           ## initial body force
    fcs = zeros(tn, 1)          ## fictitious bc-sign
    fbs = zeros(tn, 1)          ## bf bc-sign
    nnum = 0

    gc = 42320             ## critial energy release rate
    sc = sqrt(gc / (rh * (ksh * 6 / pi + 16 / (9 * (pi^2)) * (kbk - 2 * ksh))))

    dt = 1e-9                  ## time steps
    nt = 300                 ## total time step
    t0 = 8000                  ## total time step for leapfrogging
    u = zeros(tn, 2)            ## initialization displacement
    v = zeros(tn, 2)            ## initialization velocity
    tem = zeros(tn, 1)          ## initialization tempareture changes
    pforce = zeros(tn, 2)       ## initialization iner force
    pflux = zeros(tn, 1)        ## initialization thermal flux
    fail = ones(tn, hm)         ## initialization bond-state materix  1: undamaged 0: broken
    dmg = zeros(tn, 3)          ## initialization damage array 0:undamaged

    ccl = 0.02                 ## pre-existing crack length
    clo = [0, 0]                ## pre-existing crack location
    bnd = 0
    maa = 1

    ##  discrection 
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
        end
    end
    mtn = nnum                  ## material-points number

    ## temperature ficitous layers
    for i = 1:lb
        for j = 1:ny
            nnum = nnum + 1
            pin[nnum, 1] = -0.5 * (length + dx) - (i - 1) * dx
            pin[nnum, 2] = -0.5 * (width - dx) + (j - 1) * dx
            fcs[nnum, 1] = 1
        end
    end
    tnt = nnum                  ## material-points number + temperature-layer-points number 

    ## displacement ficitous layers
    for i = 1:rb
        for j = 1:ny
            nnum = nnum + 1
            pin[nnum, 1] = 0.5 * (length + dx) + (i - 1) * dx
            pin[nnum, 2] = -0.5 * (width - dx) + (j - 1) * dx
            fcs[nnum, 1] = 2
        end
    end

    tn = nnum                 ## total number of points
    tnf = mtn + tn - tnt      ## material-points number + displacement-layer-points number

    ### Horizon
    hnt, hct, hnm, hcm, ds, fail = init.Horizon(pin, tn, rh, hm, fcs, ccl, clo)

    ### Correction 
    #### thc : thermal corrcection factors
    #### mec : mechanics corrcection factors 
    thc = init.modifyth(tn, hm, pin, hnt, hct, ds, rh, dx, kp, kc)
    mec, fac = init.modifyme(tn, hm, pin, hnm, hcm, ds, rh, dx, c, emod)
    vmv = fac[1:mtn, :] * dx^3


    #### explicit 
    for tt = 1:nt
        println("Step = ", tt, " f-bonds = ", bnd)
        ## Displacement & tem-BCs
        for i = 1:tn
            if fcs[i, 1] == 1
                tem[i, 1] = 2 * 200 - tem[i-mtn, 1]      ### temperature bc      
            elseif fcs[i, 1] == 2
                u[i, :] .= 0         ### fixed bc 
            end
            ## bodt force-BCs           
            if fbs[i, 1] == 1
                if tt < t0
                    bf[i, 1] = 1e14 * tt * dt / dx
                else
                    bf[i, 1] = 1e14 * t0 * dt / dx
                end
            end
        end

        ### thermo loops
        for i = 1:mtn
            pflux[i, 1] = 0
            for j = 1:hnt[i, 1]
                m = hct[i, j]
                maa = 1
                if fcs[m, 1] == 1
                    maa = 0
                end
                yx = pin[m, 1] - pin[i, 1] + u[m, 1] - u[i, 1]
                yy = pin[m, 2] - pin[i, 2] + u[m, 2] - u[i, 2]
                ts = sqrt(yx^2 + yy^2)
                ev = (yx * (v[m, 1] - v[i, 1]) + yy * (v[m, 2] - v[i, 2])) / ts
                temvar = tem[m, 1] - tem[i, 1]
                pflux[i, 1] +=
                    (kp * temvar / ds[i, j] - ev * 0.5 * c * aph * maa) *
                    thc[i, j] *
                    (dx)^3 *
                    fail[i, j]
            end
        end
        tem += pflux * dt / (dens * cv)     ### update temperature change

        ### mechanic loops
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
