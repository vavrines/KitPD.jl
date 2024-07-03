using LinearAlgebra

# configuration
cfg = (
    np = 5, # maximum polynomial order
    lx = 1.0, # length of the domain
    nx = 100, # number of points
    ndiff = 1, # differentiation order
)

totnode = cfg.nx
nts = cfg.np + 1
tsorders = collect(0:cfg.np)

bbb = zeros(nts)
for ii = 1:nts
    bbb[ii] = factorial(tsorders[ii])
end

dx = cfg.lx / cfg.nx
δ = dx * (cfg.np + 1) # horizon size
dmag = δ
dEntity = dx
x = collect(dx/2:dx:cfg.lx-dx/2)

fi_smooth(x) = x^3
fi_shock(x) = begin
    if x < cfg.lx / 2
        1.0
    else
        0.0
    end
end
dfi_shock(x) = begin
    if cfg.lx / 2 - 1e-2 < x < cfg.lx / 2 + 1e-2
        -60.0
    else
        0.0
    end
end

fvec = zeros(totnode)
for k = 1:totnode
    #fvec[k] = fi_smooth(x[k]) # smooth
    fvec[k] = fi_shock(x[k]) # discontinuity
end

nodefam = zeros(Int, 10000)
nmax = 0
pvec = zeros(nts)
weight = zeros(nts)
dfvec = zeros(totnode)

for k = 1:totnode
    numfam = 1
    nodefam[1] = k
    for j = 1:totnode
        if j == k
            continue
        end

        idist = abs(x[j] - x[k])
        inside = true
        if idist > δ
            inside = false
        end
        if inside
            numfam += 1
            nodefam[numfam] = j
        end
    end

    if numfam > nmax
        global nmax = numfam
    end

    A = zeros(nts, nts)
    b = zeros(nts, 1)
    for kk = 1:numfam
        j = nodefam[kk]

        ξ = x[j] - x[k]
        ξmag = abs(ξ)

        for ii = 1:nts
            pvec[ii] = ξ^tsorders[ii]
            weight[ii] = exp(-4.0 * (ξmag / dmag)^2)
        end
        for ii = 1:nts
            for jj = 1:nts
                A[ii, jj] += weight[ii] * pvec[ii] * pvec[jj] * dEntity
            end
        end
    end

    for ii = 1:nts
        imatch = true
        if tsorders[ii] != cfg.ndiff
            imatch = false
        end

        if imatch
            b[ii] = bbb[ii]
            break
        end
    end

    a = A \ b

    dfval = 0.0
    for kk = 1:numfam
        j = nodefam[kk]
        ff = fvec[j]

        ξ = x[j] - x[k]
        ξmag = abs(ξ)

        gfun = 0
        for ii = 1:nts
            pvec[ii] = ξ^tsorders[ii]
            weight[ii] = exp(-4.0 * (ξmag / dmag)^2)
            gfun += a[ii] * weight[ii] * pvec[ii]
        end

        dfval += fvec[j] * gfun * dEntity
    end
    dfvec[k] = dfval
end

dfdm = zero(x)
for i in eachindex(dfdm)
    if i == 1
        dfdm[i] = (fvec[i+1] - fvec[i]) / dx
    elseif i == totnode
        dfdm[i] = (fvec[i] - fvec[i-1]) / dx
    else
        dfdm[i] = (fvec[i+1] - fvec[i-1]) / (2 * dx)
    end
end

using Plots
begin
    plot(x, dfvec, lw = 1.5, label = "PDDO", xlabel = "x", ylabel = "derivative")
    plot!(x, dfdm, lw = 1.5, label = "CDM", line = :dash)
    scatter!(x, dfi_shock.(x), alpha = 0.7, label = "exact")
    #scatter!(x, 3 .* x .^ 2, alpha = 0.7, label = "exact")
end
#savefig("smooth.png")
#savefig("sharp.png")

using BenchmarkTools

@btime begin
    for k = 1:totnode
        numfam = 1
        nodefam[1] = k
        for j = 1:totnode
            if j == k
                continue
            end
    
            idist = abs(x[j] - x[k])
            inside = true
            if idist > δ
                inside = false
            end
            if inside
                numfam += 1
                nodefam[numfam] = j
            end
        end
    
        if numfam > nmax
            global nmax = numfam
        end
    
        A = zeros(nts, nts)
        b = zeros(nts, 1)
        for kk = 1:numfam
            j = nodefam[kk]
    
            ξ = x[j] - x[k]
            ξmag = abs(ξ)
    
            for ii = 1:nts
                pvec[ii] = ξ^tsorders[ii]
                weight[ii] = exp(-4.0 * (ξmag / dmag)^2)
            end
            for ii = 1:nts
                for jj = 1:nts
                    A[ii, jj] += weight[ii] * pvec[ii] * pvec[jj] * dEntity
                end
            end
        end
    
        for ii = 1:nts
            imatch = true
            if tsorders[ii] != cfg.ndiff
                imatch = false
            end
    
            if imatch
                b[ii] = bbb[ii]
                break
            end
        end
    
        a = A \ b
    
        dfval = 0.0
        for kk = 1:numfam
            j = nodefam[kk]
            ff = fvec[j]
    
            ξ = x[j] - x[k]
            ξmag = abs(ξ)
    
            gfun = 0
            for ii = 1:nts
                pvec[ii] = ξ^tsorders[ii]
                weight[ii] = exp(-4.0 * (ξmag / dmag)^2)
                gfun += a[ii] * weight[ii] * pvec[ii]
            end
    
            dfval += fvec[j] * gfun * dEntity
        end
        dfvec[k] = dfval
    end
end

@btime begin
    for i in eachindex(dfdm)
        if i == 1
            dfdm[i] = (fvec[i+1] - fvec[i]) / dx
        elseif i == totnode
            dfdm[i] = (fvec[i] - fvec[i-1]) / dx
        else
            dfdm[i] = (fvec[i+1] - fvec[i-1]) / (2 * dx)
        end
    end
end
