using LinearAlgebra
using ProgressMeter: @showprogress

function loop1(m, n)
    global M
    global icount
    if (m > 0) && (m <= M)
        for i = 0:n
            loop1(m + 1, n - i)
            if m == M
                icount += 1
            end
        end
    end
end

function loop2(m, n, num)
    global icount
    global p
    global M
    num1 = zeros(10, 1)
    for ii = 1:10
        num1[ii] = num[ii]
    end

    if (m > 0) && (m <= M)
        for i = 0:n
            num1[m] = i
            loop2(m + 1, n - i, num1)
            if m == M
                icount += 1
                for mm = 1:M
                    p[icount, mm] = num1[mm]
                end
            end
        end
    end
end

function loop3(m, num)
    global icount
    global ndiv
    global coord
    global M
    global dx
    num1 = zeros(10, 1)
    for ii = 1:10
        num1[ii] = num[ii]
    end

    if (m > 0) && (m <= M)
        for i = 1:ndiv[m]
            num1[m] = i
            loop3(m + 1, num1)
            if m == M
                icount += 1
                for ii = 1:M
                    coord[icount, ii] = dx[ii] / 2.0 + (num1[ii] - 1) * dx[ii]
                end
            end
        end
    end
end

global icount
global p
global ndiv
global coord
global dx
global M
global N

M = 1
N = 3

len = zeros(M, 1)
ndiv = zeros(Int, M, 1)
porder = zeros(M, 1)
for i = 1:M
    len[i, 1] = 1
    ndiv[i, 1] = 100
    porder[i, 1] = 1
end

totnode = prod(ndiv)
coord = zeros(totnode, M)
icount = 0
loop1(1, N)
nsize = icount

p = zeros(Int, nsize, M)
num = zeros(10, 1)
icount = 0
loop2(1, N, num)
bb = zeros(nsize, 1)
for ii = 1:nsize
    bb[ii] = 1
    for mm = 1:M
        bb[ii] = bb[ii] * factorial(p[ii, mm])
    end
end

dmag = 0
dEntity = 1
dx = zeros(M, 1)
delta = zeros(M, 1)
for ii = 1:M
    dx[ii] = len[ii] / ndiv[ii]
    delta[ii] = dx[ii] * (N + 1)
    dmag = dmag + delta[ii] * delta[ii]
    dEntity = dEntity * dx[ii]
end
dmag = sqrt(dmag)
num = zeros(10, 1)
icount = 0

loop3(1, num)

fvec = zeros(totnode, 1)
for k = 1:totnode
    fvec[k] = coord[k, 1]^2 #+ coord[k,2]^2# + coord[k,3]^2 + coord[k,4]^2
end

nodefam = zeros(Int, 10000, 1)
idist = zeros(M, 1)
nmax = 0
#pvec = zeros(nmax,1)
#weight = zeros(nmax,1)

pvec = zeros(nsize)
weight = zeros(nsize)

dfvec = zeros(totnode, 1)
@showprogress for k = 1:totnode
    xsi = zeros(M)

    numfam = 1
    nodefam[1] = k
    for j = 1:totnode
        if j == k
            continue
        end
        for mm = 1:M
            idist[mm] = abs(coord[j, mm] - coord[k, mm])
        end
        inside = true
        for mm = 1:M
            if idist[mm] > delta[mm]
                inside = false
                break
            end
        end
        if inside
            numfam = numfam + 1
            nodefam[numfam] = j
        end
    end

    println("k = $k , numfam = $numfam")

    if (numfam > nmax)
        nmax = numfam
    end

    Amat = zeros(nsize, nsize)
    bvec = zeros(nsize, 1)
    for kk = 1:numfam
        j = nodefam[kk]
        xsimag = 0
        for mm = 1:M
            xsi[mm] = coord[j, mm] - coord[k, mm]
            xsimag = xsimag + xsi[mm] * xsi[mm]
        end
        xsimag = sqrt(xsimag)
        for ii = 1:nsize
            pvec[ii] = 1.0
            for mm = 1:M
                pvec[ii] = pvec[ii] * xsi[mm]^p[ii, mm]
            end
            weight[ii] = exp(-4 * (xsimag / dmag)^2)
        end
        for ii = 1:nsize
            for jj = 1:nsize
                Amat[ii, jj] = Amat[ii, jj] + weight[ii] * pvec[ii] * pvec[jj] * dEntity
            end
        end
    end

    AmatInv = inv(Amat)

    for ii = 1:nsize
        imatch = true
        for mm = 1:M
            if p[ii, mm] != porder[mm]
                imatch = false
                break
            end
        end
        if imatch
            bvec[ii] = bb[ii]
            break
        end
    end

    avec = AmatInv * bvec

    dfval = 0
    for kk = 1:numfam
        j = nodefam[kk]
        ff = fvec[j]
        xsimag = 0
        for mm = 1:M
            xsi[mm] = coord[j, mm] - coord[k, mm]
            xsimag = xsimag + xsi[mm] * xsi[mm]
        end
        xsimag = sqrt(xsimag)

        gfun = 0
        for ii = 1:nsize
            pvec[ii] = 1
            for mm = 1:M
                pvec[ii] = pvec[ii] * xsi[mm]^p[ii, mm]
            end
            weight[ii] = exp(-4 * (xsimag / dmag)^2)
            gfun = gfun + avec[ii] * weight[ii] * pvec[ii]
        end

        dfval = dfval + fvec[j] * gfun * dEntity
    end
    dfvec[k] = dfval
end
