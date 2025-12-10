program main

implicit none

integer ndivx, ndivy, totnode, nt, maxfam, nnum, cnode, i, j, tt, nbnd, totint, totbottom, tottop
!ndivx: Number of divisions in x direction - except boundary region
parameter(ndivx = 500)
!ndivy: Number of divisions in y direction - except boundary region
parameter(ndivy = 500)
!nbnd: Number of divisions in the boundary region
parameter(nbnd = 3)
!totnode: Total number of material points
parameter (totnode = ndivx*(ndivy + 2 * nbnd)) 
!nt: Total number of time step
parameter(nt = 1250)
!maxfam: Maximum number of material points inside a horizon of a material point
parameter(maxfam = 100)

real *8 length,width, dx, delta, thick, dens, emod, pratio, area, vol, bc  
real *8 sedload1, sedload2, dt, totime, ctime, idist, fac, radij, nlength, dforce1, dforce2 
real *8 crlength, scr0, pi, tmpdx, tmpvol, tmpcx, tmpcy, tmpux, tmpuy, dmgpar1, dmgpar2, theta 
real *8 scx, scy, scr

real *8 coord(totnode,2), pforce(totnode,2), pforceold(totnode,2), bforce(totnode,2), stendens(totnode,2)
real *8 fncst(totnode,2), disp(totnode,2), vel(totnode,2), velhalfold(totnode,2), velhalf(totnode,2)
real *8 acc(totnode,2), massvec(totnode,2), enddisp(nt,1), endtime(nt,1), dmg(totnode,1)
integer numfam(totnode,1), pointfam(totnode,1), nodefam(10000000,1), fail(totnode,maxfam)

pi = dacos(-1.0d0)

do i = 1, totnode 
    !coord: Material point locations, 1:x-coord, 2:y-coord
	coord(i,1) = 0.0d0
	coord(i,2) = 0.0d0
    !numfam: Number of family members of each material point
	numfam(i,1) = 0
    !pointfam: index array to find the family members in nodefam array
	pointfam(i,1) = 0
    !pforce: total peridynamic force acting on a material point, 1:x-coord, 2:y-coord
	pforce(i,1) = 0.0d0
	pforce(i,2) = 0.0d0
    !pforceold: total peridynamic force acting on a material point in the previous time step
    !1:x-coord, 2:y-coord
	pforceold(i,1) = 0.0d0
	pforceold(i,2) = 0.0d0
    !bforce: body load acting on a material point, 1:x-coord, 2:y-coord
	bforce(i,1) = 0.0d0
	bforce(i,2) = 0.0d0
    !stendens: strain energy of a material point, 1:loading 1, 2:loading 2
	stendens(i,1) = 0.0d0
	stendens(i,2) = 0.0d0
    !fncst: surface correction factors of a material point, 1:loading 1, 2:loading 2
	fncst(i,1) = 1.0d0 
	fncst(i,2) = 1.0d0  
    !disp: displacement of a material point, 1:x-coord, 2:y-coord
	disp(i,1) = 0.0d0
	disp(i,2) = 0.0d0
    !vel: velocity of a material point, 1:x-coord, 2:y-coord
	vel(i,1) = 0.0d0
	vel(i,2) = 0.0d0
	velhalfold(i,1) = 0.0d0
	velhalfold(i,2) = 0.0d0
	velhalf(i,1) = 0.0d0
	velhalf(i,2) = 0.0d0
    !acc: acceleration of a material point, 1:x-coord, 2:y-coord 
	acc(i,1) = 0.0d0
	acc(i,2) = 0.0d0
    !massvec: massvector for adaptive dynamic relaxation, 1:x-coord, 2:y-coord
	massvec(i,1) = 0.0d0
	massvec(i,2) = 0.0d0
    !fail: Failure array
	do j = 1, maxfam
		fail(i,j) = 0
    enddo
    !dmg: Damage of a material point
	dmg(i,1) = 0.0d0
enddo

do i = 1, 1000000
    !nodefam: array containing family members of all material points
	nodefam(i,1) = 0
enddo

!length: Total length of the plate
length = 0.05d0
!width: Total width of the plate
width = 0.05d0
!dx: Spacing between material points
dx = length / ndivx
!delta: Horizon
delta = 3.015 * dx
!thick: Thickness of the plate
thick = dx
!dens: Density
dens = 8000.0d0
!emod: Elastic modulus
emod = 192.0d9
!pratio12 = Poisson's ratio
pratio = 1.0d0 / 3.0d0
!area: Cross-sectional area
area = dx * dx
!vol: Volume of a material point
vol = area * dx
!bc: Bond constant 
bc = 9.0d0 * emod / (pi * thick * (delta**3))
!sedload1: Strain energy density for the first loading
sedload1 = 9.0d0 / 16.0d0 * emod * 1.0d-6   
!sedload2: Strain energy density for the second loading
sedload2 = 9.0d0 / 16.0d0 * emod * 1.0d-6
!dt: Time interval
dt = 0.8d0 * dsqrt(2.0d0*dens*dx/(pi*delta**2*dx*bc))
!totime: Total time
totime = nt * dt
!ctime: Current time
ctime = 0.0d0
!idist: Initial distance
idist = 0.0d0
do i = 1, nt
	enddisp(i,1) = 0.0d0
	endtime(i,1) = 0.0d0
enddo
!fac: Volume correction factor
fac = 0.0d0
!radij: Material point radius
radij = dx / 2.0d0
!nnum: Material point number
nnum = 0
!cnode: Current material point
cnode = 0
!Length of deformed bond
nlength  = 0.0d0
!dforce1: x component of the PD force between two material points
dforce1 = 0.0d0
!dforce1: y component of the PD force between two material points
dforce2 = 0.0d0
!crlength: Crack length
crlength = 0.01d0
!scr0: Critical stretch
scr0 = 0.04472d0

!Initialization of fail flag array
!1 means no failure, 0 means failure of the PD bond
do i = 1,totnode
	do j = 1,maxfam
		fail(i,j) = 1
    enddo
enddo

!Specification of the locations of material points
!Material points of the internal region
do i = 1,ndivy
    do j = 1,ndivx
        nnum = nnum + 1
        coord(nnum,1) = (-1.0d0 * length / 2.0d0) + (dx / 2.0d0) + (j-1) * dx
        coord(nnum,2) = (-1.0d0 * width / 2.0d0) + (dx / 2.0d0) + (i-1) * dx
    enddo
enddo

totint = nnum

!Material points of the boundary region - bottom
do i = 1,nbnd
    do j = 1,ndivx
        nnum = nnum + 1
        coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
        coord(nnum,2) = -1.0d0 /2.0d0 * width - (dx / 2.0d0) - (i - 1) * dx
    enddo
enddo

totbottom = nnum

!Material points of the boundary region - top
do i = 1,nbnd
    do j = 1,ndivx
        nnum = nnum + 1
        coord(nnum,1) = -1.0d0 /2.0d0 * length + (dx / 2.0d0) + (j - 1) * dx
        coord(nnum,2) = 1.0d0 /2.0d0 * width + (dx / 2.0d0) + (i - 1) * dx
    enddo
enddo

tottop = nnum

!Determination of material points inside the horizon of each material point
do i = 1,totnode
    if (i.eq.1) then 
        pointfam(i,1) = 1
    else
        pointfam(i,1) = pointfam(i-1,1) + numfam(i-1,1)
    endif
    do j = 1,totnode
        idist = dsqrt((coord(j,1) - coord(i,1))**2 + (coord(j,2) - coord(i,2))**2)
        if (i.ne.j) then
            if(idist.le.delta) then
                numfam(i,1) = numfam(i,1) + 1
                nodefam(pointfam(i,1)+numfam(i,1)-1,1) = j
            endif
        endif
    enddo
enddo

!Definition of the crack surface
!PD bonds penetrating through the crack surface are broken
do i = 1,totnode
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
        if ((coord(cnode,2) > 0.0d0).and.(coord(i,2) < 0.0d0)) then
            if ((dabs(coord(i,1)) - (crlength / 2.0d0)).le.1.0d-10) then
				fail(i,j) = 0
            elseif ((dabs(coord(cnode,1)) - (crlength / 2.0d0)).le.1.0d-10) then
				fail(i,j) = 0
            endif
        elseif ((coord(i,2) > 0.0d0).and.(coord(cnode,2) < 0.0d0)) then
            if((dabs(coord(i,1)) - (crlength / 2.0d0)).le.1.0d-10) then 
				fail(i,j) = 0
			elseif((dabs(coord(cnode,1)) - (crlength / 2.0d0)).le.1.0e-10) then
				fail(i,j) = 0
            endif
        endif        
    enddo
enddo

!Loading 1
do i = 1,totnode
    disp(i,1) = 0.001d0 * coord(i,1)
    disp(i,2) = 0.0d0
enddo

do i = 1,totnode
    stendens(i,1) = 0.0d0
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        nlength = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif
                       
        stendens(i,1) = stendens(i,1) + 0.5d0 * 0.5d0 * bc * ((nlength - idist) / idist)**2 * idist * vol * fac  
    enddo
    !Calculation of surface correction factor in x direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncst(i,1) = sedload1 / stendens(i,1)
enddo
    
!Loading 2
do i = 1,totnode
    disp(i,1) = 0.0d0
    disp(i,2) = 0.001d0 * coord(i,2)
enddo

do i = 1,totnode
    stendens(i,2) = 0.0d0
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        nlength = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif   
                      
        stendens(i,2) = stendens(i,2) + 0.5d0 * 0.5d0 * bc * ((nlength - idist) / idist)**2 * idist * vol * fac 
    enddo
    !Calculation of surface correction factor in y direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncst(i,2) = sedload2 / stendens(i,2)
enddo
    
!Initialization of displacements and velocities
do i = 1,totnode
    vel(i,1) = 0.0d0
    vel(i,2) = 0.0d0
    disp(i,1) = 0.0d0
    disp(i,2) = 0.0d0         
enddo

!Time integration
do tt = 1,nt
    write(*,*) 'tt = ', tt
	ctime = tt * dt

    !Application of boundary conditions at the top and bottom edges
    do i = (totint+1), totbottom
        vel(i,2) = -20.0d0
        disp(i,2) = -20.0d0 * tt * dt
    enddo

    do i = (totbottom+1), tottop
        vel(i,2) = 20.0d0
        disp(i,2) = 20.0d0 * tt * dt
    enddo    
    
    do i = 1,totnode
        dmgpar1 = 0.0d0
        dmgpar2 = 0.0d0
        pforce(i,1) = 0.0d0
        pforce(i,2) = 0.0d0
        do j = 1,numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1)
                idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
                nlength = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
                !Volume correction
                if (idist.le.delta-radij) then
                    fac = 1.0d0
                elseif (idist.le.delta+radij) then
                    fac = (delta+radij-idist)/(2.0d0*radij)
                else
                    fac = 0.0d0
                endif
                if (dabs(coord(cnode,2) - coord(i,2)).le.1.0d-10) then 
                    theta = 0.0d0
                elseif (dabs(coord(cnode,1) - coord(i,1)).le.1.0d-10) then
                    theta = 90.0d0 * pi / 180.0d0
                else
                    theta = datan(dabs(coord(cnode,2) - coord(i,2)) / dabs(coord(cnode,1) - coord(i,1)))
                endif
                !Determination of the surface correction between two material points
                scx = (fncst(i,1) + fncst(cnode,1)) / 2.0d0
                scy = (fncst(i,2) + fncst(cnode,2)) / 2.0d0
                scr = 1.0d0 / (((dcos(theta))**2 / (scx)**2) + ((dsin(theta))**2 / (scy)**2))
                scr = dsqrt(scr)
                
                if (fail(i,j).eq.1) then
                    !Calculation of the peridynamic force in x and y directions 
                    !acting on a material point i due to a material point j
                    dforce1 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) / nlength             
                    dforce2 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) / nlength             
                else
                	dforce1 = 0.0d0
                    dforce2 = 0.0d0
                endif 
                pforce(i,1) = pforce(i,1) + dforce1             
                pforce(i,2) = pforce(i,2) + dforce2             
                
                !Definition of a no-fail zone             
                if (dabs((nlength - idist) / idist) > scr0) then
                    if (dabs(coord(i,2)).le.(length/4.0d0)) then
						fail(i,j) = 0 
					endif
                endif 					 
                            
                dmgpar1 = dmgpar1 + fail(i,j) * vol * fac
                dmgpar2 = dmgpar2 + vol * fac             
        enddo
        !Calculation of the damage parameter
        dmg(i,1) = 1.0d0 - dmgpar1 / dmgpar2
    enddo
    
    do i = 1,totint
        !Calculation of acceleration of material point i
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        acc(i,2) = (pforce(i,2) + bforce(i,2)) / dens
        !Calculation of velocity of material point i
        !by integrating the acceleration of material point i
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        vel(i,2) = vel(i,2) + acc(i,2) * dt
        !Calculation of displacement of material point i
        !by integrating the velocity of material point i
        disp(i,1) = disp(i,1) + vel(i,1) * dt
        disp(i,2) = disp(i,2) + vel(i,2) * dt
    enddo
    
    do i = (totint+1), totbottom
        acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
        vel(i,1) = vel(i,1) + acc(i,1) * dt
        disp(i,1) = disp(i,1) + vel(i,1) * dt
    enddo

    do i = (totbottom+1), tottop
            acc(i,1) = (pforce(i,1) + bforce(i,1)) / dens
            vel(i,1) = vel(i,1) + acc(i,1) * dt
            disp(i,1) = disp(i,1) + vel(i,1) * dt
    enddo
               
    endtime(tt,1) = ctime

	if (tt.eq.750) then
        !printing results to an output file
		open(26,file = 'coord_disp_pd_750_pwc_v20.txt')

		do i = 1, totint
			write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
		enddo

		close(26)
	elseif (tt.eq.1000) then
		open(26,file = 'coord_disp_pd_1000_pwc_v20.txt')

		do i = 1, totint
			write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
		enddo

		close(26)
	elseif (tt.eq.1250) then
		open(26,file = 'coord_disp_pd_1250_pwc_v20.txt')

		do i = 1, totint
			write(26,111) coord(i,1), coord(i,2), disp(i,1), disp(i,2), dmg(i,1)
		enddo

		close(26)
	endif

enddo

111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)

end program main