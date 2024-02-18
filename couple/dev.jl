cd(@__DIR__)
include("init.jl")

length = 0.1 # meter
width = 0.1  # meter
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
tn = (nx + lb + rb) * (ny + ub + db)

# initialization
pin = zeros(tn, 2) # initial position
bf = zeros(tn, 2) # initial body force
fcs = zeros(tn, 1) # fictitious bc-sign
fbs = zeros(tn, 1) # bf bc-sign

isedge = Array{Bool}(undef, tn)
isedge .= false

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
mtn = nnum                  ## material-points number

ids = findall(x -> x == true, isedge)
for id in ids
    cx = pin[nnum, 1]
    cy = pin[nnum, 2]
    @assert cx < dx - 0.5 * length ||
    cx > 0.5 * length - dx ||
    cy < dy - 0.5 * width ||
    cy > 0.5 * width - dy
end

findall(x -> x == 1, fbs)
