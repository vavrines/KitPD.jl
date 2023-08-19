! ------------------------------------------------------------
! 13.6.3 Plate Subjected to a Shock of Pressure and Temperature, 
! and Their Combination
! ------------------------------------------------------------

program main

implicit none

integer ndivx, ndivy, totnode, nt, maxfam, nnum, cnode, i, j, tt, nbnd, N, np, mp
!ndivx: Number of divisions in x direction - except boundary region
parameter(ndivx = 200)
parameter(ndivy = 200)
!nbnd: Number of divisions in the boundary region
parameter(nbnd = 3)
!totnode: Total number of material points
parameter (totnode = ndivy * (ndivx + 2 * nbnd)) 
!nt: Total number of time step
parameter(nt = 3000)
!maxfam: Maximum number of material points inside a horizon of a material point
parameter(maxfam = 100)
!N: tot array 
parameter(N = maxfam * totnode * 2)

real *8 vartem, dx, length, delta, area, vol, is, stenload1, stenload2, thick, theta, width
real *8 tempcload1, dt, totime, ctime, idist, fac, radij, dflux1, dforce1, dforce2
real *8 pi, ta, tb, ma, cn, cn1, cn2, cc, ev
real *8 scrt, scrm, scr, scx, scy, coordx, coordy

real *8 coord(totnode,2), pflux(totnode,1), tempddens(totnode,1), stendens(totnode,2), bforce(totnode,2)
real *8 fncsth(totnode,1), fncstm(totnode,2), tem(totnode,1), disp(totnode,2), pdtem(nt,1)
real *8 vel(totnode,2), acc(totnode,2), pforce(totnode,2), hs(totnode,1)
integer numfam(totnode,1), pointfam(totnode,1), nodefam(N,1)

pi = dacos(-1.0d0)

!coord: Material point locations
do i = 1, totnode 
    !coord: Material point locations
	coord(i,1) = 0.0d0
	coord(i,2) = 0.0d0
    !numfam: Number of family members of each material point
	numfam(i,1) = 0
    !pointfam: index array to find the family members in nodefam array
	pointfam(i,1) = 0
    !pflux: total peridynamic force acting on a material point
	pflux(i,1) = 0.0d0
	pforce(i,1) = 0.0d0
	bforce(i,1) = 0.0d0
	pforce(i,2) = 0.0d0
	bforce(i,2) = 0.0d0
    !tempddens: strain energy of a material point
	tempddens(i,1) = 0.0d0
	stendens(i,1) = 0.0d0
	stendens(i,2) = 0.0d0
    !fncst: surface correction factor of a material point
	fncsth(i,1) = 1.0d0
	fncstm(i,1) = 1.0d0
	fncstm(i,2) = 1.0d0
    !tem: tem of a material point    
	tem(i,1) = 0.0d0
	disp(i,1) = 0.0d0
	vel(i,1) = 0.0d0
	acc(i,1) = 0.0d0
	disp(i,2) = 0.0d0
	vel(i,2) = 0.0d0
	acc(i,2) = 0.0d0
       hs(i,1) = 0.0d0
enddo

do i = 1, N
    !nodefam: array containing family members of all material points
	nodefam(i,1) = 0
enddo

!length: Total length of the bar
length = 10.0d0
width = 10.0d0
!dx: Spacing between material points
dx = length / ndivx
thick = dx
!delta: Horizon
delta = 3.015d0 * dx
!area: Cross-sectional area
area = dx * dx
!ta
ta = 6.0d0 / (delta**3) / pi / thick
ma = ta
cc = 0.10d0
tb = cc / 2.0d0
!vol: Volume of a material point
vol = area * dx
!tempcload1: thermal potential of a material point for the first thermal condition
!based on classical continuum mechanics
tempcload1 =  1 * 1.0d-6  
stenload1 = 9.0d0 / 16.0d0 * 1 * 1.0d-6
stenload2 = 9.0d0 / 16.0d0 * 1 * 1.0d-6
!dt: Time interval
dt = 1.0d-3
!totime: Total time
totime = nt * dt
!ctime: Current time
ctime = 0.0d0
!idist: Initial distance
idist = 0.0d0
!fac: Volume correction factor
fac = 0.0d0
!radij: Material point radius
radij = dx / 2.0d0
!nnum: Material point number
nnum = 0




!thm boundaries of the plate
do i = 1,nbnd
    do j = 1,ndivy
        coordx = -1.0d0 /2.0d0 * dx - (i - 1) * dx
        coordy = -1.0d0 /2.0d0 * width + (dx / 2.0d0) + (j - 1) * dx
        nnum = nnum + 1
        coord(nnum,1) = coordx
        coord(nnum,2) = coordy
    enddo
enddo
np = nnum

!Specification of the locations of material points
!Material points of the plate
do i = 1,ndivx
    do j = 1,ndivy
        coordx = 1.0d0 /2.0d0 * dx + (i - 1) * dx
        coordy = -1.0d0 /2.0d0 * width + (dx / 2.0d0) + (j - 1) * dx
        nnum = nnum + 1
        coord(nnum,1) = coordx
        coord(nnum,2) = coordy
    enddo
enddo
mp = nnum

!dsp boundaries of the plate
do i = 1,nbnd
    do j = 1,ndivy
        coordx = length + (dx / 2.0d0) + (i - 1) * dx
        coordy = -1.0d0 /2.0d0 * width + (dx / 2.0d0) + (j - 1) * dx
        nnum = nnum + 1
        coord(nnum,1) = coordx
        coord(nnum,2) = coordy
    enddo
enddo

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

!Determination of surface correction factors 
! Tem Loading  0.
do i = 1,totnode
    tem(i,1) = 0.001d0 * (coord(i,1) + coord(i,2))
enddo

do i = 1,totnode
     tempddens(i,1) = 0.0d0
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        vartem = tem(cnode,1) - tem(i,1)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif            
        tempddens(i,1) =  tempddens(i,1) + 0.25d0 * ta * vartem**2 / idist * vol * fac 

    enddo
    !Calculation of surface correction factor in x direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncsth(i,1) = tempcload1 / tempddens(i,1)
enddo

! mechanism Loading -x 
do i = 1,totnode
    disp(i,1) = 0.001d0 * coord(i,1)
    disp(i,2) = 0.0d0
enddo

do i = 1,totnode
    stendens(i,1) = 0.0d0
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        is = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
(coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif            
        stendens(i,1) = stendens(i,1) + 3.0d0 / 8.0d0 * ma * ((is - idist) / idist)**2 * idist * vol * fac  
    enddo
    !Calculation of surface correction factor in x direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncstm(i,1) = stenload1 / stendens(i,1)
enddo

! mechanism Loading -yy 
do i = 1,totnode
    disp(i,1) = 0.0d0
    disp(i,2) = 0.001d0 * coord(i,2)
enddo

do i = 1,totnode
    stendens(i,2) = 0.0d0
    do j = 1,numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1)
        idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
        is = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
(coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)
        if (idist.le.delta-radij) then
            fac = 1.0d0
        elseif (idist.le.delta+radij) then
            fac = (delta+radij-idist)/(2.0d0*radij)
        else
            fac = 0.0d0
        endif            
        stendens(i,2) = stendens(i,2) + 3.0d0 / 8.0d0 * ma * ((is - idist) / idist)**2 * idist * vol * fac  
    enddo
    !Calculation of surface correction factor in x direction 
    !by finding the ratio of the analytical strain energy density value
    !to the strain energy density value obtained from PD Theory
    fncstm(i,2) = stenload2 / stendens(i,2)
enddo

open(34,file = 'fncst.txt')
do i = 1,totnode
    !Store therem and disp information for the material point- all
    write(34,112) coord(i,1), coord(i,2), fncsth(i,1), fncstm(i,1), fncstm(i,2)
    112 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
enddo
close(34)

!Initialization of tem & disp 
do i = 1,totnode
   tem(i,1) = 0.0d0 
   disp(i,1) = 0.0d0
   disp(i,2) = 0.0d0    
enddo

!Time integration
do tt = 1,nt
    write(*,*) 'tt = ', tt

   !Boundary condition  tem shock left
   do i = 1, np
      tem(i,1) = 5.0d0 * dt * tt * exp(-2.0d0 * dt * tt) 
   enddo

    do i = 1,totnode
        pflux(i,1) = 0.0d0
        do j = 1,numfam(i,1)            
                cnode = nodefam(pointfam(i,1)+j-1,1)
                idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
                vartem = tem(cnode,1) - tem(i,1)
                is = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
(coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2) 
                ev = (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) / is * (vel(cnode,1) - vel(i,1)) + &
(coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) / is * (vel(cnode,2) - vel(i,2)) 
                           
                !Volume correction
                if (idist.le.delta-radij) then
                    fac = 1.0d0
                elseif (idist.le.delta+radij) then
                    fac = (delta+radij-idist)/(2.0d0*radij)
                else
                    fac = 0.0d0
                endif

                !Determination of the surface correction between two material points
                scrt = (fncsth(i,1) + fncsth(cnode,1)) / 2.0d0
              
                !Calculation of the peridynamic force in x direction 
                !acting on a material point i due to a material point j
                dflux1 = ta * (vartem / idist - 0.5d0 * cc * ev) * vol * scrt * fac            
                pflux(i,1) = pflux(i,1) + dflux1                                                    				          
        enddo
    enddo

    do i = 1,totnode
        tem(i,1) = tem(i,1) + dt * (pflux(i,1) + hs(i,1))
   enddo

   !!!!!!!!!!~~~~~~~~loop--mechanics--disp

   !Boundary condition  constriction-right
   do i = (mp + 1), totnode
      disp(i,1) = 0.0d0
      disp(i,2) = 0.0d0
      vel(i,1) = 0.0d0
      vel(i,2) = 0.0d0
   enddo

   !Boundary condition  presu shock left
   do i = 1, np
      bforce(i,1) = 5.0d0 * dt * tt * exp(-2.0d0 * dt * tt) / (dx * nbnd)
      bforce(i,2) = 0.0d0 
   enddo

    do i = 1,totnode
        pforce(i,1) = 0.0d0
        pforce(i,2) = 0.0d0
        do j = 1,numfam(i,1)   
                cnode = nodefam(pointfam(i,1)+j-1,1)
                idist = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)
                is = dsqrt((coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1))**2 + &
(coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2))**2)

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
                scx = (fncstm(i,1) + fncstm(cnode,1)) / 2.0d0
                scy = (fncstm(i,2) + fncstm(cnode,2)) / 2.0d0
                scr = 1.0d0 / (((dcos(theta))**2.0d0 / (scx)**2.0d0) + ((dsin(theta))**2.0d0 / (scy)**2.0d0))
                scrm = dsqrt(scr)
                    
                !Calculation of the peridynamic force in x direction 
                !acting on a material point i due to a material point j
                dforce1 = ma * (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) / is *&
 vol * scrm * fac * ((is - idist) / idist * 4.0d0 / 3.0d0 - (tem(cnode,1) + tem(i,1)) / 2)            
                pforce(i,1) = pforce(i,1) + dforce1 

                dforce2 = ma * (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) / is *&
 vol * scrm * fac * ((is - idist) / idist * 4.0d0 / 3.0d0 - (tem(cnode,1) + tem(i,1)) / 2)            
                pforce(i,2) = pforce(i,2) + dforce2 
                                                   				          
        enddo
    enddo

    do i = 1,totnode
        vel(i,1) = vel(i,1) + dt * (pforce(i,1) + bforce(i,1))
        disp(i,1) = disp(i,1) + dt * vel(i,1)
        vel(i,2) = vel(i,2) + dt * (pforce(i,2) + bforce(i,2))
        disp(i,2) = disp(i,2) + dt * vel(i,2)
    enddo
enddo

open(35,file = '1dc.txt')
do i = (np + ndivy / 2),mp,ndivy
    !Store therem and disp information for the material point at the y = 0
    write(35,111) coord(i,1), tem(i,1), disp(i,1), disp(i,2)
    111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
enddo
close(35)

open(36,file = '1all.txt')
do i = (np + 1),mp
    !Store therem and disp information for the material point- all
    write(36,222) coord(i,1), tem(i,1), disp(i,1), disp(i,2)
    222 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
enddo
close(36)

end program main
