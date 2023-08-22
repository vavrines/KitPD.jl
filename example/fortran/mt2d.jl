"""
13.6.3 Plate Subjected to a Shock of Pressure and Temperature

Reminder:
- hs: heat source
- ta & ma: Dimensionless coefficients in Eq. (13.65)

Question:

"""

using KitBase.ProgressMeter: @showprogress
using Base.Threads: @threads

begin
    ndivx = 200
    ndivy = 200
    nbnd = 3
    totnode = ndivy * (ndivx + 2 * nbnd)
    nt = 3000
    maxfam = 100
    N = maxfam * totnode * 2
end

begin
    coord = zeros(totnode, 2)
    numfam = zeros(Int, totnode, 1)
    pointfam = zeros(Int, totnode, 1)
    nodefam = zeros(Int, N, 1)
    pflux = zeros(totnode, 1)
    pforce = zeros(totnode, 2)
    bforce = zeros(totnode, 2)
    tempddens = zeros(totnode, 1)
    stendens = zeros(totnode, 2)
    fncsth = zeros(totnode, 1)
    fncstm = zeros(totnode, 2)
    tem = zeros(totnode, 1)
    disp = zeros(totnode, 2)
    vel = zeros(totnode, 2)
    acc = zeros(totnode, 2)
    hs = zeros(totnode, 1)
end

begin
    length = 10.0
    width = 10.0
    dx = length / ndivx
    thick = dx
    delta = 3.015 * dx
    area = dx * dx
    ta = 6.0 / (delta^3) / pi / thick
    ma = ta
    cc = 0.10
    tb = cc / 2.0
    vol = area * dx
    tempcload1 = 1e-6
    stenload1 = 9.0 / 16.0 * 1 * 1e-6
    stenload2 = 9.0 / 16.0 * 1 * 1e-6
    dt = 1e-3
    totime = nt * dt
    ctime = 0.0
    idist = 0.0
    fac = 0.0
    radij = dx / 2.0
    nnum = 0
end

for i = 1:nbnd
    for j = 1:ndivy
        coordx = -1.0 / 2.0 * dx - (i - 1) * dx
        coordy = -1.0 / 2.0 * width + (dx / 2.0) + (j - 1) * dx
        nnum = nnum + 1
        coord[nnum, 1] = coordx
        coord[nnum, 2] = coordy
    end
end
np = nnum

for i = 1:ndivx
    for j = 1:ndivy
        coordx = 1.0 / 2.0 * dx + (i - 1) * dx
        coordy = -1.0 / 2.0 * width + (dx / 2.0) + (j - 1) * dx
        nnum = nnum + 1
        coord[nnum, 1] = coordx
        coord[nnum, 2] = coordy
    end
end
mp = nnum

for i = 1:nbnd
    for j = 1:ndivy
        coordx = length + (dx / 2.0) + (i - 1) * dx
        coordy = -1.0 / 2.0 * width + (dx / 2.0) + (j - 1) * dx
        nnum = nnum + 1
        coord[nnum, 1] = coordx
        coord[nnum, 2] = coordy
    end
end

function sort_family!(numfam, pointfam, coord, delta)
    totnode = size(pointfam, 1)

    for i = 1:totnode
        if i == 1
            pointfam[i, 1] = 1
        else
            pointfam[i, 1] = pointfam[i-1, 1] + numfam[i-1, 1]
        end

        @threads for j = 1:totnode
            idist = sqrt((coord[j, 1] - coord[i, 1])^2 + (coord[j, 2] - coord[i, 2])^2)
            if i != j && idist <= delta
                numfam[i, 1] += 1
                nodefam[pointfam[i, 1]+numfam[i, 1]-1, 1] = j
            end
        end
    end

    return nothing
end

sort_family!(numfam, pointfam, coord, delta)

for i = 1:totnode
    tem[i, 1] = 0.001 * (coord[i, 1] + coord[i, 2])
end

@showprogress for i = 1:totnode
    tempddens[i, 1] = 0.0
    for j = 1:numfam[i, 1]
        cnode = nodefam[pointfam[i, 1]+j-1, 1]
        idist = sqrt((coord[cnode, 1] - coord[i, 1])^2 + (coord[cnode, 2] - coord[i, 2])^2)
        vartem = tem[cnode, 1] - tem[i, 1]
        if idist <= delta - radij
            fac = 1.0
        elseif idist <= delta + radij
            fac = (delta + radij - idist) / (2.0 * radij)
        else
            fac = 0.0
        end
        tempddens[i, 1] += 0.25 * ta * vartem^2 / idist * vol * fac

    end
    fncsth[i, 1] = tempcload1 / tempddens[i, 1]
end

for i = 1:totnode
    disp[i, 1] = 0.001 * coord[i, 1]
    disp[i, 2] = 0.0
end

@showprogress for i = 1:totnode
    stendens[i, 1] = 0.0
    for j = 1:numfam[i, 1]
        cnode = nodefam[pointfam[i, 1]+j-1, 1]
        idist = sqrt((coord[cnode, 1] - coord[i, 1])^2 + (coord[cnode, 2] - coord[i, 2])^2)
        is = sqrt(
            (coord[cnode, 1] + disp[cnode, 1] - coord[i, 1] - disp[i, 1])^2 +
            (coord[cnode, 2] + disp[cnode, 2] - coord[i, 2] - disp[i, 2])^2,
        )
        if idist <= delta - radij
            fac = 1.0
        elseif idist <= delta + radij
            fac = (delta + radij - idist) / (2.0 * radij)
        else
            fac = 0.0
        end
        stendens[i, 1] += 3.0 / 8.0 * ma * ((is - idist) / idist)^2 * idist * vol * fac
    end

    fncstm[i, 1] = stenload1 / stendens[i, 1]
end

for i = 1:totnode
    disp[i, 1] = 0.0
    disp[i, 2] = 0.001 * coord[i, 2]
end

@showprogress for i = 1:totnode
    stendens[i, 2] = 0.0
    for j = 1:numfam[i, 1]
        cnode = nodefam[pointfam[i, 1]+j-1, 1]
        idist = sqrt((coord[cnode, 1] - coord[i, 1])^2 + (coord[cnode, 2] - coord[i, 2])^2)
        is = sqrt(
            (coord[cnode, 1] + disp[cnode, 1] - coord[i, 1] - disp[i, 1])^2 +
            (coord[cnode, 2] + disp[cnode, 2] - coord[i, 2] - disp[i, 2])^2,
        )
        if (idist <= delta - radij)
            fac = 1.0
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij)
        else
            fac = 0.0
        end
        stendens[i, 2] += 3.0 / 8.0 * ma * ((is - idist) / idist)^2 * idist * vol * fac
    end

    fncstm[i, 2] = stenload2 / stendens[i, 2]
end

for i = 1:totnode
    tem[i, 1] = 0.0
    disp[i, 1] = 0.0
    disp[i, 2] = 0.0
end

@showprogress for tt = 1:5#nt
    #--- temperature ---#
    for i = 1:np
        tem[i, 1] = 5.0 * dt * tt * exp(-2.0 * dt * tt)
    end

    for i = 1:totnode
        pflux[i, 1] = 0.0
        for j = 1:numfam[i, 1]
            cnode = nodefam[pointfam[i, 1]+j-1, 1]
            idist =
                sqrt((coord[cnode, 1] - coord[i, 1])^2 + (coord[cnode, 2] - coord[i, 2])^2)
            vartem = tem[cnode, 1] - tem[i, 1]
            is = sqrt(
                (coord[cnode, 1] + disp[cnode, 1] - coord[i, 1] - disp[i, 1])^2 +
                (coord[cnode, 2] + disp[cnode, 2] - coord[i, 2] - disp[i, 2])^2,
            )
            ev =
                (coord[cnode, 1] + disp[cnode, 1] - coord[i, 1] - disp[i, 1]) / is *
                (vel[cnode, 1] - vel[i, 1]) +
                (coord[cnode, 2] + disp[cnode, 2] - coord[i, 2] - disp[i, 2]) / is *
                (vel[cnode, 2] - vel[i, 2])

            if idist <= delta - radij
                fac = 1.0
            elseif idist <= delta + radij
                fac = (delta + radij - idist) / (2.0 * radij)
            else
                fac = 0.0
            end

            scrt = (fncsth[i, 1] + fncsth[cnode, 1]) / 2.0
            dflux1 = ta * (vartem / idist - 0.5 * cc * ev) * vol * scrt * fac

            pflux[i, 1] += dflux1
        end
    end

    for i = 1:totnode
        tem[i, 1] += dt * (pflux[i, 1] + hs[i, 1])
    end

    #--- mechanics ---#
    for i = (mp+1):totnode
        disp[i, 1] = 0.0
        disp[i, 2] = 0.0
        vel[i, 1] = 0.0
        vel[i, 2] = 0.0
    end

    for i = 1:np
        bforce[i, 1] = 5.0 * dt * tt * exp(-2.0 * dt * tt) / (dx * nbnd)
        bforce[i, 2] = 0.0
    end

    for i = 1:totnode
        pforce[i, 1] = 0.0
        pforce[i, 2] = 0.0
        for j = 1:numfam[i, 1]
            cnode = nodefam[pointfam[i, 1]+j-1, 1]
            idist =
                sqrt((coord[cnode, 1] - coord[i, 1])^2 + (coord[cnode, 2] - coord[i, 2])^2)
            is = sqrt(
                (coord[cnode, 1] + disp[cnode, 1] - coord[i, 1] - disp[i, 1])^2 +
                (coord[cnode, 2] + disp[cnode, 2] - coord[i, 2] - disp[i, 2])^2,
            )

            if idist <= delta - radij
                fac = 1.0
            elseif idist <= delta + radij
                fac = (delta + radij - idist) / (2.0 * radij)
            else
                fac = 0.0
            end

            if (abs(coord[cnode, 2] - coord[i, 2]) <= 1e-10)
                theta = 0.0
            elseif (abs(coord[cnode, 1] - coord[i, 1]) <= 1e-10)
                theta = 90.0 * pi / 180.0
            else
                theta = atan(
                    abs(coord[cnode, 2] - coord[i, 2]) / abs(coord[cnode, 1] - coord[i, 1]),
                )
            end

            scx = (fncstm[i, 1] + fncstm[cnode, 1]) / 2.0
            scy = (fncstm[i, 2] + fncstm[cnode, 2]) / 2.0
            scr = 1.0 / (((cos(theta))^2.0 / (scx)^2) + ((sin(theta))^2.0 / (scy)^2))
            scrm = sqrt(scr)

            dforce1 =
                ma * (coord[cnode, 1] + disp[cnode, 1] - coord[i, 1] - disp[i, 1]) / is *
                vol *
                scrm *
                fac *
                ((is - idist) / idist * 4.0 / 3.0 - (tem[cnode, 1] + tem[i, 1]) / 2)
            pforce[i, 1] += dforce1

            dforce2 =
                ma * (coord[cnode, 2] + disp[cnode, 2] - coord[i, 2] - disp[i, 2]) / is *
                vol *
                scrm *
                fac *
                ((is - idist) / idist * 4.0 / 3.0 - (tem[cnode, 1] + tem[i, 1]) / 2)
            pforce[i, 2] += dforce2
        end
    end

    for i = 1:totnode
        vel[i, 1] += dt * (pforce[i, 1] + bforce[i, 1])
        disp[i, 1] += dt * vel[i, 1]
        vel[i, 2] += dt * (pforce[i, 2] + bforce[i, 2])
        disp[i, 2] += dt * vel[i, 2]
    end
end
