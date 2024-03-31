"""
Test of the equivalence of original and modified recursive functions\
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
            icount = locate_grid!(coord, m + 1, ndiv, dx, num1, M, icount)
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

M = 2
N = 3

# test num_ts
icount = 0
loop1(1, N)
icount == num_ts(1, N, M)

# test assign_ts!
p = zeros(Int, icount, M)
num = zeros(10, 1)
icount = 0
loop2(1, N, num)

p1 = zero(p)
icount == assign_ts!(p1, 1, N, num, M)
p == p1

# test locate_grid!
ndiv = zeros(Int, M, 1)
dx = zeros(M, 1)
for i = 1:M
    ndiv[i, 1] = 100
    dx[i] = 1 / ndiv[i]
end
totnode = prod(ndiv)
coord = zeros(totnode, M)

icount = 0
loop3(1, num)

coord1 = zero(coord)
locate_grid!(coord1, 1, ndiv, dx, num, M)
coord == coord1
