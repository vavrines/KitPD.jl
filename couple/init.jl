#### Initialization function
#### Horizon: horizon and adjacent relationships
#### modifyth：correction factors for thermo 
#### modifyme：correction factors for mechanic

module init
export Horizon, modifyth, modifyme

## pin: initial position
## tn: Total number of points
## rh: horizon radius
## hm： maximum number of points in horizon
## fcs: differentiate the thermo or mechanic boundary with a sign, where 1 represents thermo and 2 represents mechanics.
## ccl: pre-existing crack length
## clo: pre-existing crack location
## hornt: Array of the accurate quantities of material points in each horizon for thermo
## horct: Index matrix of each points in horizon for thermo
## hornm: Array of the accurate quantities of material points in each horizon for mechanic
## horcm: Index matrix of each points in horizon for mechanic
## dsini: Bonds' length at initial position 
## fail: pre-existing crack, 1:undamaged 0:broken

function Horizon(pin, tn, rh, hm, fcs, ccl, clo)
    hornt = zeros(Int, tn, 1)
    horct = zeros(Int, tn, hm)
    hornm = zeros(Int, tn, 1)
    horcm = zeros(Int, tn, hm)
    dsini = zeros(tn, hm)
    fail = ones(tn, hm)

    for i = 1:tn
        nd = 0
        ndt = 0
        ndm = 0
        for j = 1:tn
            if j != i
                dst = sqrt((pin[j, 1] - pin[i, 1])^2 + (pin[j, 2] - pin[i, 2])^2)
                if dst < rh
                    nd += 1
                    dsini[i, nd] = dst
                    if fcs[i, 1] != 2 && fcs[j, 1] != 2  ## remove the interecation of mechanic bc points 
                        ### Here, I excluded the excess displacement boundary fictitious layer during thermal calculation. 
                        ### Mark 2 indicates the PD points located at the displacement boundary.
                        ### The array horct does not include points at the displacement boundary.
                        ### This array is used for computing the surface correction for heat transfer and thermal conduction calculations.
                        ndt += 1
                        horct[i, ndt] = j
                    end

                    if fcs[i, 1] != 1 && fcs[j, 1] != 1  ## remove the interecation of thermo bc points
                        ### Here, I excluded the excess temperature boundary fictitious layer during thermal calculation. 
                        ### Mark 1 indicates the PD points located at the displacement boundary.
                        ### The array horcm does not include points at the temperature boundary.
                        ### This array is used for computing the surface correction for deformation.
                        ndm += 1
                        horcm[i, ndm] = j
                        if (pin[i, 1] - clo[1]) * (pin[j, 1] - clo[1]) < 0
                            if abs(pin[i, 2] - clo[2]) <= 0.5 * ccl ||
                               abs(pin[j, 2] - clo[2]) <= 0.5 * ccl
                                fail[i, ndm] = 0
                            end
                        end
                    end
                end
                hornt[i, 1] = ndt
                hornm[i, 1] = ndm
            end
        end
    end
    return hornt, horct, hornm, horcm, dsini, fail
end

##### tnt: Total number of points using in thermo 
### hnt as hornt; hct as horct; ds as dsini
### dx: space size
### kp:  thermal microconductivity in PD
### Kc: thermal conducticity in classical theory
### tmf: volume correction * surface correction of thermo

function modifyth(tn, hm, pin, hnt, hct, ds, rh, dx, kp, kc)
    tem = zeros(tn, 1) # temperature initialization
    factm = zeros(tn, hm) # volume correction factor initialization in thermo 
    thf = zeros(tn, 1) # surface correction factor initialization of thermo
    tmen = zeros(tn, 1) # thermoal potential initialization
    tmf = zeros(tn, hm) # total correction in thermo
        ### Here we perform surface and volume corrections separately, as the domains for heat and force are different. This function is for thermal。 

    tem[:, 1] = 0.001 * (pin[:, 1] + pin[:, 2])

    for i = 1:tn
        thg = 0
        for j = 1:hnt[i, 1]
            k = hct[i, j]
            ### volume modify 
            if ds[i, j] < (rh - 0.5 * dx)
                factm[i, j] = 1
            elseif ds[i, j] < (rh + 0.5 * dx)
                factm[i, j] = (rh + 0.5 * dx - ds[i, j]) / dx
            else
                factm[i, j] = 0
            end

            ### surface modify in x 
            thg += 0.25 * kp * (tem[k, 1] - tem[i, 1])^2 / ds[i, j] * factm[i, j] * (dx)^3 #thermo potential in PD
        end

        tmen[i, 1] = thg
        thf[i, 1] = (kc * 1e-6) / tmen[i, 1]
    end

    ### surface direction
    for i = 1:tn
        for j = 1:hnt[i, 1]
            k = hct[i, j]
            tmf[i, j] = 0.5 * (thf[i, 1] + thf[k, 1]) * factm[i, j]
        end
    end
    return tmf
end

###### tnf: Total number of points using in motion 
### hnm as hornm; hcm as horcm; ds as dsini
### dx: space size
### c:  micromudule in PD
### emod: module in classical theory
### fmf: volume correction * surface correction of motion

function modifyme(tn, hm, pin, hnm, hcm, ds, rh, dx, c, emod)
    factf = zeros(tn, hm) # volume correction factor initialization in motion
    dfmx = zeros(tn, 2) # deformation in x--initialization
    dfmy = zeros(tn, 2) # deformation in y--initialization
    mof = zeros(tn, 2) # surface correction factor initialization of motion
    sten = zeros(tn, 2) # strain energy density initialization
    fmf = zeros(tn, hm) # total correction in motion
     ### Here we perform surface and volume corrections separately, as the domains for heat and force are different. This function is for motion。 


    dfmx[:, 1] = 1.001 * pin[:, 1]
    dfmx[:, 2] = 1.000 * pin[:, 2]
    dfmy[:, 2] = 1.001 * pin[:, 2]
    dfmy[:, 1] = 1.000 * pin[:, 1]

    for i = 1:tn
        eng1 = 0
        eng2 = 0
        for j = 1:hnm[i, 1]
            k = hcm[i, j]
            ### volume modify 
            if ds[i, j] < (rh - 0.5 * dx)
                factf[i, j] = 1
            elseif ds[i, j] < (rh + 0.5 * dx)
                factf[i, j] = (rh + 0.5 * dx - ds[i, j]) / dx
            else
                factf[i, j] = 0
            end

            ### surface modify in x 
            tes1 = sqrt((dfmx[k, 1] - dfmx[i, 1])^2 + (dfmx[k, 2] - dfmx[i, 2])^2)
            eng1 = eng1 + 0.25 * c * (tes1 - ds[i, j])^2 / ds[i, j] * factf[i, j] * (dx)^3 #strain energy density in PD
            tes2 = sqrt((dfmy[k, 1] - dfmy[i, 1])^2 + (dfmy[k, 2] - dfmy[i, 2])^2)
            eng2 = eng2 + 0.25 * c * (tes2 - ds[i, j])^2 / ds[i, j] * factf[i, j] * (dx)^3 #strain energy density in PD
        end

        sten[i, 1] = eng1
        sten[i, 2] = eng2
        mof[i, 1] = (9 * emod * 1e-6) / (16 * sten[i, 1])
        mof[i, 2] = (9 * emod * 1e-6) / (16 * sten[i, 2])
    end

    ### surface direction
    scr = zeros(tn, hm)
    for i = 1:tn
        for j = 1:hnm[i, 1]
            k = hcm[i, j]
            if pin[k, 2] ≈ pin[i, 2]
                theta = 0
            elseif pin[k, 1] ≈ pin[i, 1]
                theta = π / 2
            else
                theta = atan(abs((pin[k, 2] - pin[i, 2]) / (pin[k, 1] - pin[i, 1])))
            end
            sx = (mof[k, 1] + mof[i, 1]) / 2
            sy = (mof[k, 2] + mof[i, 2]) / 2
            scr[i, j] = (cos(theta) / sx)^2 + (sin(theta) / sy)^2
            scr[i, j] = sqrt(1 / scr[i, j])
            fmf[i, j] = scr[i, j] * factf[i, j]
        end
    end
    return fmf, factf
end

end
