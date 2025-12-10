"""
9.2 Plate with a Pre-existing Crack Under Velocity Boundary Conditions
Translated from Fortran source codes
"""

# Import necessary libraries
using LinearAlgebra
using Printf

begin
    # Define constants (parameters)
    ndivx = 100#500          # Number of divisions in x direction (excluding boundary)
    ndivy = 100#500          # Number of divisions in y direction (excluding boundary)
    nbnd = 3             # Number of boundary region divisions
    totnode = ndivx * (ndivy + 2 * nbnd)  # Total material points
    nt = 1250            # Total time steps
    maxfam = 100         # Max family members per material point

    # Material properties and simulation parameters
    length = 0.05        # Total plate length (m)
    width = 0.05         # Total plate width (m)
    dx = length / ndivx  # Spacing between material points
    delta = 3.015 * dx   # Horizon (interaction radius)
    thick = dx           # Plate thickness
    dens = 8000.0        # Density (kg/m³)
    emod = 192.0e9       # Elastic modulus (Pa)
    pratio = 1.0 / 3.0   # Poisson's ratio
    area = dx * dx       # Cross-sectional area
    vol = area * dx      # Volume per material point
    # Bond constant (micromodulus)
    bc = 9.0 * emod / (pi * thick * delta^3)  
    # Strain energy density for load cases
    sedload1 = 9.0 / 16.0 * emod * 1.0e-6
    sedload2 = 9.0 / 16.0 * emod * 1.0e-6
    # Time step (stable based on CFL condition)
    dt = 0.8 * sqrt(2.0 * dens * dx / (pi * delta^2 * dx * bc))
    totime = nt * dt     # Total simulation time
    crlength = 0.01      # Crack length (m)
    scr0 = 0.04472       # Critical stretch for bond failure
end

begin
    # Preallocate arrays
    coord = zeros(Float64, totnode, 2)        # Coordinates [x, y]
    pforce = zeros(Float64, totnode, 2)       # Peridynamic force [x, y]
    pforceold = zeros(Float64, totnode, 2)    # Previous force [x, y]
    bforce = zeros(Float64, totnode, 2)       # Body force [x, y]
    stendens = zeros(Float64, totnode, 2)     # Strain energy density [load1, load2]
    fncst = ones(Float64, totnode, 2)         # Surface correction factors [x, y]
    disp = zeros(Float64, totnode, 2)         # Displacement [x, y]
    vel = zeros(Float64, totnode, 2)          # Velocity [x, y]
    velhalfold = zeros(Float64, totnode, 2)   # Half-step velocity (old) [x, y]
    velhalf = zeros(Float64, totnode, 2)      # Half-step velocity [x, y]
    acc = zeros(Float64, totnode, 2)          # Acceleration [x, y]
    massvec = zeros(Float64, totnode, 2)      # Mass vector [x, y] (for adaptive dynamics)
    dmg = zeros(Float64, totnode)             # Damage parameter
    numfam = zeros(Int, totnode)              # Number of family members
    pointfam = zeros(Int, totnode)            # Index pointer for family list
    nodefam = zeros(Int, 10000000)            # Family member list (large preallocation)
    fail = ones(Int, totnode, maxfam)         # Bond failure flags (1 = intact, 0 = broken)

    # Initialize arrays
    ctime = 0.0  # Current time
    enddisp = zeros(Float64, nt)  # Displacement at end nodes
    endtime = zeros(Float64, nt)  # Time at each step

    # Initialize counters for region tracking
    nnum = 0      # Total node counter
    totint = 0    # Last internal node index
    totbottom = 0 # Last bottom boundary node index
    tottop = 0    # Last top boundary node index
end

# ----- Material Point Setup -----
# Internal region
for i in 1:ndivy
    for j in 1:ndivx
        nnum += 1
        coord[nnum, 1] = -length/2 + dx/2 + (j-1)*dx  # x-coordinate
        coord[nnum, 2] = -width/2 + dx/2 + (i-1)*dx   # y-coordinate
    end
end
totint = nnum  # Last internal node index

# Bottom boundary region
for i in 1:nbnd
    for j in 1:ndivx
        nnum += 1
        coord[nnum, 1] = -length/2 + dx/2 + (j-1)*dx   # x-coordinate
        coord[nnum, 2] = -width/2 - dx/2 - (i-1)*dx    # y-coordinate (below plate)
    end
end
totbottom = nnum  # Last bottom boundary node index

# Top boundary region
for i in 1:nbnd
    for j in 1:ndivx
        nnum += 1
        coord[nnum, 1] = -length/2 + dx/2 + (j-1)*dx   # x-coordinate
        coord[nnum, 2] = width/2 + dx/2 + (i-1)*dx     # y-coordinate (above plate)
    end
end
tottop = nnum  # Last top boundary node index

# ----- Family Identification -----
# Find all material points within horizon delta for each point
pointfam[1] = 1  # Starting index for first point
for i in 1:totnode
    println("Finding family for point $i of $totnode")
    if i > 1
        # Set starting index based on previous point
        pointfam[i] = pointfam[i-1] + numfam[i-1]
    end
    
    # Check all other points
    for j in 1:totnode
        if i != j
            # Calculate distance between points i and j
            idist = norm(coord[j, :] - coord[i, :])
            if idist <= delta
                numfam[i] += 1
                # Store neighbor index in nodefam
                nodefam[pointfam[i] + numfam[i] - 1] = j
            end
        end
    end
end

# ----- Pre-existing Crack Setup -----
# Break bonds crossing the crack line (x=0 between y=-crlength/2 to y=crlength/2)
for i in 1:totnode
    println("Setting crack for point $i of $totnode")
    for j in 1:numfam[i]
        cnode = nodefam[pointfam[i] + j - 1]  # Current neighbor
        
        # Check if bond crosses the crack plane (y=0)
        if (coord[cnode, 2] > 0 && coord[i, 2] < 0) ||
            (coord[i, 2] > 0 && coord[cnode, 2] < 0)
            
            # Check if either point is near the crack tip
            if (abs(coord[i, 1]) - crlength/2 <= 1e-10) ||
                (abs(coord[cnode, 1]) - crlength/2 <= 1e-10)
                fail[i, j] = 0  # Break the bond
            end
        end
    end
end

# ----- Surface Correction Factor Calculation -----
# Loading Case 1: Uniform x-strain (0.1%)
for i in 1:totnode
    disp[i, 1] = 0.001 * coord[i, 1]
end
for i in 1:totnode
    stendens[i, 1] = 0.0
    for j in 1:numfam[i]
        cnode = nodefam[pointfam[i] + j - 1]
        # Reference distance
        idist = norm(coord[cnode, :] - coord[i, :])
        # Deformed distance
        nlength = norm(coord[cnode, :] + disp[cnode, :] - coord[i, :] - disp[i, :])
        
        # Volume correction factor (partial volume treatment)
        if idist <= delta - dx/2
            fac = 1.0
        elseif idist <= delta + dx/2
            fac = (delta + dx/2 - idist) / dx
        else
            fac = 0.0
        end
        
        # Strain energy density contribution
        stendens[i, 1] += 0.5 * 0.5 * bc * ((nlength - idist)/idist)^2 * idist * vol * fac
    end
    # Surface correction factor for x-direction
    fncst[i, 1] = sedload1 / stendens[i, 1]
end

# Loading Case 2: Uniform y-strain (0.1%)
for i in 1:totnode
    disp[i, 1] = 0.0
    disp[i, 2] = 0.001 * coord[i, 2]
end
for i in 1:totnode
    stendens[i, 2] = 0.0
    for j in 1:numfam[i]
        cnode = nodefam[pointfam[i] + j - 1]
        idist = norm(coord[cnode, :] - coord[i, :])
        nlength = norm(coord[cnode, :] + disp[cnode, :] - coord[i, :] - disp[i, :])
        
        if idist <= delta - dx/2
            fac = 1.0
        elseif idist <= delta + dx/2
            fac = (delta + dx/2 - idist) / dx
        else
            fac = 0.0
        end
        
        stendens[i, 2] += 0.5 * 0.5 * bc * ((nlength - idist)/idist)^2 * idist * vol * fac
    end
    # Surface correction factor for y-direction
    fncst[i, 2] = sedload2 / stendens[i, 2]
end

# Reset displacements and velocities
disp .= 0.0
vel .= 0.0

# ----- Time Integration Loop -----
for tt in 1:nt
    ctime = tt * dt
    @printf("tt = %d\n", tt)

    # Apply velocity boundary conditions
    # Bottom boundary: constant downward velocity
    for i in (totint+1):totbottom
        vel[i, 2] = -20.0
        disp[i, 2] = -20.0 * ctime
    end
    # Top boundary: constant upward velocity
    for i in (totbottom+1):tottop
        vel[i, 2] = 20.0
        disp[i, 2] = 20.0 * ctime
    end

    # Reset forces and damage parameters
    pforce .= 0.0
    for i in 1:totnode
        dmgpar1 = 0.0  # Intact bond measure
        dmgpar2 = 0.0  # Total bond measure
        
        # Process all family members
        for j in 1:numfam[i]
            cnode = nodefam[pointfam[i] + j - 1]
            # Reference position vector
            ref_vec = coord[cnode, :] - coord[i, :]
            idist = norm(ref_vec)
            # Deformed position vector
            def_vec = (coord[cnode, :] + disp[cnode, :]) - (coord[i, :] + disp[i, :])
            nlength = norm(def_vec)
            
            # Volume correction factor
            if idist <= delta - dx/2
                fac = 1.0
            elseif idist <= delta + dx/2
                fac = (delta + dx/2 - idist) / dx
            else
                fac = 0.0
            end
            
            # Bond angle for anisotropic correction
            dx_diff = abs(coord[cnode, 1] - coord[i, 1])
            dy_diff = abs(coord[cnode, 2] - coord[i, 2])
            if dy_diff < 1e-10
                theta = 0.0
            elseif dx_diff < 1e-10
                theta = pi/2
            else
                theta = atan(dy_diff, dx_diff)
            end
            
            # Average surface correction and bond direction correction
            scx = (fncst[i, 1] + fncst[cnode, 1]) / 2.0
            scy = (fncst[i, 2] + fncst[cnode, 2]) / 2.0
            scr = 1.0 / sqrt((cos(theta)^2 / scx^2) + (sin(theta)^2 / scy^2))
            
            # Calculate peridynamic force if bond is intact
            if fail[i, j] == 1
                stretch = (nlength - idist) / idist
                # Force magnitude
                force_mag = bc * stretch * vol * scr * fac / nlength
                # Force components
                dforce1 = force_mag * def_vec[1]
                dforce2 = force_mag * def_vec[2]
            else
                dforce1 = 0.0
                dforce2 = 0.0
            end
            
            # Add to total force
            pforce[i, 1] += dforce1
            pforce[i, 2] += dforce2
            
            # Check bond failure condition (outside no-fail zone)
            stretch_ratio = abs(nlength - idist) / idist
            if stretch_ratio > scr0
                # Only break bonds outside the central region
                if abs(coord[i, 2]) > length/4
                    fail[i, j] = 0
                end
            end
            
            # Update damage parameters
            dmgpar1 += fail[i, j] * vol * fac
            dmgpar2 += vol * fac
        end
        
        # Calculate damage value (0 = intact, 1 = fully broken)
        dmg[i] = 1.0 - dmgpar1 / dmgpar2
    end

    # Update internal points (explicit Euler integration)
    for i in 1:totint
        # Acceleration = Force / Density (per volume)
        acc[i, 1] = pforce[i, 1] / dens
        acc[i, 2] = pforce[i, 2] / dens
        
        # Velocity update
        vel[i, 1] += acc[i, 1] * dt
        vel[i, 2] += acc[i, 2] * dt
        
        # Displacement update
        disp[i, 1] += vel[i, 1] * dt
        disp[i, 2] += vel[i, 2] * dt
    end

    # Update boundary regions (x-direction only)
    for i in (totint+1):totbottom
        acc[i, 1] = pforce[i, 1] / dens
        vel[i, 1] += acc[i, 1] * dt
        disp[i, 1] += vel[i, 1] * dt
    end
    for i in (totbottom+1):tottop
        acc[i, 1] = pforce[i, 1] / dens
        vel[i, 1] += acc[i, 1] * dt
        disp[i, 1] += vel[i, 1] * dt
    end

    # Record time
    endtime[tt] = ctime

    # Output results at specified time steps
    if tt in [750, 1000, 1250]
        filename = "coord_disp_pd_$(tt)_pwc_v20.txt"
        open(filename, "w") do f
            for i in 1:totint
                @printf(f, "%.5e   %.5e   %.5e   %.5e   %.5e\n",
                        coord[i, 1], coord[i, 2], 
                        disp[i, 1], disp[i, 2], 
                        dmg[i])
            end
        end
    end
end
