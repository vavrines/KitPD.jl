"""
unfinished 2D work
"""

using LinearAlgebra

function num_ts(m, n, M, icount = 0)
    if m > 0 && m <= M
        for i = 0:n
            icount = num_ts(m + 1, n - i, M, icount)
            if m == M
                icount += 1
            end
        end
    end

    return icount
end

function assign_ts!(p, m, n, num, M, icount = 0)
    num1 = deepcopy(num)

    if m > 0 && m <= M
        for i = 0:n
            num1[m] = i
            icount = assign_ts!(p, m + 1, n - i, num1, M, icount)
            if m == M
                icount += 1
                for mm = 1:M
                    p[icount, mm] = num1[mm]
                end
            end
        end
    end

    return icount
end

function locate_grid!(coord, m, ndiv, dx, num, M, icount = 0)
    num1 = deepcopy(num)

    if m > 0 && m <= M
        for i = 1:ndiv[m]
            num1[m] = i
            locate_grid!(coord, m + 1, ndiv, dx, num1, M, icount)
            if m == M
                icount += 1
                for ii = 1:M
                    coord[icount, ii] = 0.5 * dx[ii] + (num1[ii] - 1) * dx[ii]
                end
            end
        end
    end

    return icount
end

# configuration
cfg = (
    np = 3, # maximum polynomial order
    ndim = 1, # number of dimensions
    lx = 1.0, # length of the domain
    nx = [100], # number of points
    ndiff = 1, # differentiation order
)

totnode = prod(cfg.nx)
x = zeros(totnode, cfg.ndim)
nts = num_ts(1, cfg.np, cfg.ndim)

num_ts(1, 3, 1)


tsorders = zeros(Int, nts, cfg.ndim)
assign_ts!(tsorders, 1, cfg.np, zeros(10), cfg.ndim)

bb = zeros(nts)
for ii = 1:nts
    bb[ii] = 1
    for mm = 1:cfg.ndim
        bb[ii] = bb[ii] * factorial(tsorders[ii, mm])
    end
end

dmag = 0
dEntity = 1
dx = zeros(cfg.ndim)
δ = zeros(cfg.ndim)
for ii = 1:cfg.ndim
    dx[ii] = len[ii] / cfg.nx[ii]
    δ[ii] = dx[ii] * (cfg.np + 1) # horizon size
    dmag += δ[ii]^2
    dEntity *= dx[ii]
end
dmag = sqrt(dmag)

locate_grid!(x, 1, cfg.nx, dx, zeros(10), cfg.ndim)

fvec = zeros(totnode, 1)
for k = 1:totnode
    fvec[k] = x[k, 1]^2
end

nodefam = zeros(10000)
idist = zeros(cfg.ndim)
nmax = 0
pvec = zeros(nts)
weight = zeros(nts)
