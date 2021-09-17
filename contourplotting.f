    subroutine Pancake

    include 'pinc/particles.f'
    include 'pinc/secondary.f'

    real pi
    integer uvbins i, j, k, a, b, BinIndex
    parameter ( pi = 3.141592653589793238462643, uvbins = 30 )
    real vp, vr, rad, phi, theta, radv, agelo, agehi, metlo, methi
    logical good, limits
    character*30 lbl

    common /ToContour/ values(uvbins,uvbins)
    real values, vmin, vmax, umin, umax, rmin, rmax
    real tr(6), lo, hi, disp(3), dis(3), vel(3)
    real ulow, uhigh, vlow, vhigh

c Asks user for new centre location and the radial velocity
    print*, "Please enter new centre Xc, Yc, Zc"
    read*, disp(1), disp(2), disp(3)
    print*, "Please enter the tangential velocity (m/s)"
    read*, radv 
    call AngleCalc (disp(1),disp(2),phi)
    
c splits the galactic velocity (related to the local standard of rest
c into relative x, y, z components based off the angle calculated above 
c from the subroutine (z = 0 always)
    vel(1) = radv*sin(phi)
    vel(2) = radv*cos(phi)
    vel(3) = 0

c translates the galactic centre and velocity centre
c based off the user inputs
    call ShiftOrigin (vel,.true.)
    call ShiftOrigin (disp,.false.)

c asks the user for other parameters of input, limits or filters on the range
c of data to be displayed, has the option to choose not to have some so it
c will auto place high range values for them that would cover any realistic simulation
    print*,"Enter Radial Range: rmin, rmax"
    read*, rmin, rmax
    print*,"Need more limits? (true/false)"
    read*, limits
    if ( .not. limits) then
        agelo = 0
        agehi = 100
        metlo = 0
        methi = 100
        ulow = -10000
        uhigh = 10000
        vlow = -10000
        vhigh = 10000
    else
        print*, "Enter Age Range: agemin, agemax"
        read*, agelo, agehi
        print*, "Enter Metallicity Range: metmin, metmax"
        read*, metlo, methi
        print*, "Enter U-Velocity Range: umin, umax"
        read*, ulow, uhigh
        print*, "Enter V-Velocity Range: vmin, vmax"
        read*, vlow, vhigh
    endif
    vmax = 0.
    vmin = 10000000.
    umax = 0.
    umin = 10000000.

c first loop through the star data to find the upper/ lower limits 
c of the data, which is effectively decided by the user inputs if 
c velocity limits were used
    do i = Ngas+Ndark+1, Ngas+Ndark+Nstar
        rad = sqrt( x(i)**2 + y(i)**2 )
        if( (rad .ge. rmin) .and. (rad .le. rmax) ) then
            if( (tform(i) .ge. agelo) .and. (tform(i) .le. agehi) ) then
                if ((metals(i) .ge. metlo) .and. (metals(i) .le. methi)) then
                    vr = (x(i)*vx(i) + y(i)*vy(i))/rad
                    vp = (-y(i)*vx(i) + x(i)*vy(i))/rad 
                    if ((vr .ge. ulow) .and. (vr .le. uhigh)) then
                        if ((vp .ge. vlow) .and. (vp .le. vhigh)) then 
                            vmax = max( vmax, vp )
                            vmin = min( vmin, vp )
                            umax = max( umax, vr )
                            umin = min( umin, vr )
                        end if
                    end if
                end if
            end if
        end if
    end do

c initialize the density array
    do i = 1, uvbins
        do j = 1, uvbins
            values(i,j) = 0.
        end do
    end do

c loop over all particles, k = counter
    k = 0
    do j = Ngas+Ndark+1, Ngas+Ndark+Nstar
        rad = sqrt( x(j)*x(j) + y(j)*y(j) )
        if( (rad .ge. rmin) .and. (rad .le. rmax) ) then
            if( (tform(j) .ge. agelo) .and. (tform(j) .le. agehi) )then
                if ((metals(j) .ge. metlo) .and. (metals(j) .le. methi))then
                    vr = (x(j)*vx(j) + y(j)*vy(j))/rad
                    vp = (-y(j)*vx(j) + x(j)*vy(j))/rad
                    if ((vr .ge. ulow) .and. (vr .le. uhigh)) then
                        if ((vp .ge. vlow) .and. (vp .le. vhigh)) then
                            a = BinIndex(vmin,vmax,uvbins,vp)
                            b = BinIndex(umin,umax,uvbins,vr)
                            if( (a .ne. 0) .and. (a .le. uvbins) .and. (b .ne. 0) .and. (b .le. uvbins) ) then
                                values(a,b) = values(a,b)
                                k = k + 1
                            end if
                        end if
                    end if
                end if
            end if
        end if
    end do
    print*,"Used ",k," particles."
    call RefineLims(vmin,vmax,umax,umin,values)
    
c re-initialize the density array
    do i = 1, uvbins
        do j = 1, uvbins
            values(i,j) = 0.
        end do
    end do

c loop again over all particles.
    do j = Ngas+Ndark+1, Ngas+Ndark+Nstar
        rad = sqrt( x(j)*x(j) + y(j)*y(j) )
        if( (rad .ge. rmin) .and. (rad .le. rmax) ) then
            if( (tform(j) .ge. agelo) .and. (tform(j) .le. agehi) ) then
                if ((metals(j) .ge. metlo) .and. (metals(j) .le. methi)) then
                    vr = (x(j)*vx(j) + y(j)*vy(j))/rad
                    vp = (-y(j)*vx(j) + x(j)*vy(j))/rad
                    if ((vr .ge. ulow) .and. (vr .le. uhigh)) then
                        if ((vp .ge. vlow) .and. (vp .le. vhigh)) then
                            a = BinIndex(vmin,vmax,uvbins,vp)
                            b = BinIndex(umin,umax,uvbins,vr)
                            if( (a .ne. 0) .and. (a .le. uvbins) .and. (b .ne. 0) .and. (b .le. uvbins) ) then
                                values(a,b) = values(a,b)
                            end if
                        end if
                    end if
                end if
            end if
        end if
    end do

c go to log space
    lo = 1000000000000.
    hi = -lo
    do b = 1, uvbins
        do a = 1, uvbins 
            if(values(a,b) .eq. 0.)then
                values(a,b) = -666.
            else
                values(a,b) = log10(values(a,b))
                lo = min( lo, values(a,b) )
                hi = max( hi, values(a,b) )
            end if
        end do
    end do 

c plotting routines
    call PGWINDOW( vmin, vmax, umin, umax )
    call GetContConsts(vmin,vmax,umin,umax,uvbins,uvbins,tr)
    call ContourPlot(tr,.false.,.false.,lo,hi)
    a = rmin*10
    b = -1
    call PGNUMB(a,b,1,lbl,i)
    lbl(i+1:30) = ' \\(2243) R \\(2243) '
    i = i + 20
    a = rmax*10
    b = -1
    call PGNUMB(a,b,1,lbl(i:30),j)
    call PGTEXT( vmin, umax*1.05, lbl )
    call PGBOX( 'BCNST', 0., 0, 'BCNST', 0., 0 )
    call PGLAB( 'V', 'U', ' ' )
    call LevelMeter(lo,hi,0.2,0.8,.false.,.false.) 
    return
    end
c End of subroutine Pancake

c co-ord and velocity translation subroutine depending on the if check
    subroutine ShiftOrigin (disp,vel)

    include 'pinc/particles.f'
    
    real disp(3)
    logical vel
    integer i

    do i= 1, Ntot
        if ( vel ) then 
            vx(i) = vx(i) - disp(1)
            vy(i) = vy(i) - disp(2)
            vz(i) = vz(i) - disp(3)
        else 
            x(i) = x(i) - disp(1)
            y(i) = y(i) - disp(2)
            z(i) = z(i) - disp(3)
        end if 
    end do
    return
    end
c end of subroutine ShiftOrigin

c calculates the angle based on an x-y co-ordinate relative to another
c position, with the x-axis = 0 degrees and going clockwise round from there
    subroutine AngleCalc(x,y,angle)

    parameter ( pi = 3.141592653589793238462643 )
    real x, y, angle 

    if (x .gt. 0) then
        if (y .gt. 0) then
            angle = 2*pi - atan(x/y)
        else
            if (y .eq. 0) then
                angle = 0 
            else
                angle = atan(x/y)
            endif
        endif
    else
        if (y .gt. 0) then
            if (x .eq. 0) then
                angle = 3*(pi/2)
            else
                angle = pi + atan(x/y)
            endif
        else
            if (y .eq. 0) then
                angle = pi
            else
                if (x .eq. 0) then
                    angle = pi/2
                else
                    angle = pi - atan(x/y)
                endif
            endif
        endif
    endif 
    return
    end
c end of subroutine AngleCalc 

c plotting subroutine
    subroutine ContourPlot (trans,vlike,dens,lo,hi)
    
    integer i, stp, uvbins
    parameter ( stp = 2, uvbins = 30 )
    real trans(6), lo, hi, h, c1, c2
    logical vlike, dens

    common /ToContour/ values(uvbins,uvbins)
    real values

    if( .not. dens )then
        do i = 1, 240, stp
            c1 = lo + real(i)*(hi-lo)/240.
            c2 = lo + real(i+stp)*(hi-lo)/240.
            if( vlike )then
                h = 120. + i
            else
                h = 360. - i
            end if
            call PGSHLS(12,h,0.5,1.0)
            call PGSCI(12)
            call PGCONF(values,uvbins,uvbins,1,uvbins,1,uvbins,c1,c2,trans)
        end do
    else
        call PGGRAY(values,uvbins,uvbins,1,uvbins,1,uvbins,lo,hi,trans)
    end if
    return
    end
c end of subroutine ContourPlot

c this is existing code at the project start provided by Victor 
c from here onwards, very little changes made here

c this subroutine simply plots up a level meter of two levels lo and hi 
c left is a logical that tells the code to use either the left hand side 
c or the right hand side.
    subroutine LevelMeter(lo,hi,y1,y2,vlike,left)

    real c1, c2, lo, hi, levels(2,240), h, tr(6), y1, y2
    logical vlike, left
    integer i

c draw the level meter
    call PGSHLS(12,0.0,0.0,0.0)
    call PGSCI(12)
    call PGSLS( 1 )
    call PGSLW( 1 )
    call PGSCH( 1.0 ) 
    if( left )then
        call PGSVP( -0.05, -0.01, y1, y2 )
    else
        call PGSVP( 0.82, 0.85, y1, y2 )
    end if
    call PGWINDOW( 0.3, 0.7, lo, hi )
    do i = 1, 240
        levels(1,i) = lo + real(i)*(hi-lo)/240.
        levels(2,i) = lo + real(i)*(hi-lo)/240.
    end do
    call GetContConsts(0.,1.,lo,hi,2,240,tr)
    do i = 1, 240
        c1 = lo + real(i)*(hi-lo)/240.
        c2 = lo + real(i+2)*(hi-lo)/240.
        if( vlike )then
            h = 120. + i
        else
            h = 360. - i
        end if
        call PGSHLS(12,h,0.5,1.0)
        call PGSCI(12)
        call PGCONB(levels,2,240,1,2,1,240,c1,1,tr,-666.)
        call PGCONF(levels,2,240,1,2,1,240,c1,c2,tr)
    end do
    call PGSHLS(12,0.0,0.0,0.0)
    call PGSCI(12)
    call PGSCH( 1.0 )
    if( left )then
        call PGBOX( 'BC', 0., 0, 'BCNTV', 0., 0 )
    else
        call PGBOX( 'BC', 0., 0, 'BCMTV', 0., 0 )
    end if
    return
    end
c end of subroutine LevelMeter
    
c change limits so I don't waste energy in regions where there are few particles
    subroutine RefineLims(vmin,vmax,umax,umin,values)

    integer i, j, vlo, vhi, uhi, uvbins
    parameter ( uvbins = 30 )
    real vmin, vmax, umax, umin, values(uvbins,uvbins)
    real vt1, vt2, mtot, mtmp

c squeeze out the excess space in the x-axis
    mtot = 0.
    do i = 1, uvbins
        do j = 1, uvbins
            mtot = mtot + values(i,j)
        end do
    end do
    mtmp = 0.
    vlo = 0
    do while( mtmp .lt. 0.001*mtot )
        vlo = vlo + 1
        do j = 1, uvbins
            mtmp = mtmp + values(vlo,j)
        end do
    end do
    if( vlo .gt. 1 )vlo = vlo - 1
        mtmp = 0.
        vhi = uvbins + 1
        do while( mtmp .lt. 0.001*mtot )
            vhi = vhi - 1
            do j = 1, uvbins
                mtmp = mtmp + values(vhi,j)
            end do
        end do
        if( vhi .lt. uvbins )vhi = vhi + 1
            vt1 = vmin + (vlo-1)*(vmax-vmin)/real(101)
            vt2 = vmax - (uvbins-vhi)*(vmax-vmin)/real(101)
            vmin = vt1
            vmax = vt2
c squeeze out the excess space in the y axis
            mtmp = 0.
            vlo = 0
            do while( mtmp .lt. 0.001*mtot )
                vlo = vlo + 1
                do j = 1, uvbins
                    mtmp = mtmp + values(j,vlo)
                end do
            end do
            if( vlo .gt. 1 )vlo = vlo - 1
                mtmp = 0.
                vhi = uvbins + 1
                do while( mtmp .lt. 0.001*mtot )
                    vhi = vhi - 1
                    do j = 1, uvbins
                        mtmp = mtmp + values(j,vhi)
                    end do
                end do
                if( vhi .lt. uvbins )vhi = vhi + 1
                    vt1 = umin + (vlo-1)*(umax-umin)/real(101)
                    vt2 = umax - (uvbins-vhi)*(umax-umin)/real(101)
                    umin = vt1
                    umax = vt2

                return
                end
c end of subroutine RefineLims