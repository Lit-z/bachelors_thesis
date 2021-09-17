    subroutine Pancake

    include 'pinc/particles.f'
    include 'pinc/secondary.f'

    real pi
    integer i, j, k, a, b 
    parameter ( pi = 3.141592653589793238462643, uvbins = 30 ) 
    real vp, vr, rad, phi, theta, radv, agelo, agehi, metlo, methi
    logical good
    character*30 lbl

    real values, vmin, vmax, umin, umax, rmin, rmax 
    real tr(6), lo, hi, disp(3), dis(3), vel(3) 
    real up(Nstar), ur(Nstar)


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
c of data to be displayed, doesn’t give as many choices as the contour plot
c as limiting the velocity won’t change a point plot in the same way it 
c changes a contour plot so if wanted to look at point plot between a certain
c velocity range zooming in would work just as effectively
    print*,"Enter Radial Range: rmin, rmax"
    read*, rmin, rmax
    print*,"Need more limits? (true/false)"
    read*, limits
    print*, "Enter Age Range: agemin, agemax"
    read*, agelo, agehi
    print*, "Enter Metallicity Range: metmin, metmax"
    read*, metlo, methi
    vmax = 0.
    vmin = 10000000.
    umax = 0.
    umin = 10000000.
    
c first loop through the star data to find the upper/ lower limits 
c of the data
    do i = Ngas+Ndark+1, Ngas+Ndark+Nstar
        rad = sqrt( x(i)**2 + y(i)**2 )
        if( (rad .ge. rmin) .and. (rad .le. rmax) )then
            if( (tform(i) .ge. agelo) .and. (tform(i) .le. agehi) )then
                if ((metals(i) .ge. metlo) .and. (metals(i) .le. methi))then
                    vr = (x(i)*vx(i) + y(i)*vy(i))/rad
                    vp = (-y(i)*vx(i) + x(i)*vy(i))/rad
                    vmax = max( vmax, vp )
                    vmin = min( vmin, vp )
                    umax = max( umax, vr )
                    umin = min( umin, vr )
                end if
            end if
        end if
    end do

c initialize the points array
    do i = 1, Nstar
        ur(i) = 0
        up(i) = 0
    end do

C loop over all particles, k = counter
    k = 0
    do j = Ngas+Ndark+1, Ngas+Ndark+Nstar
        rad = sqrt( x(j)*x(j) + y(j)*y(j) )
        if( (rad .ge. rmin) .and. (rad .le. rmax) )then
            if( (tform(j) .ge. agelo) .and. (tform(j) .le. agehi) )then
                if ((metals(j) .ge. metlo) .and. (metals(j) .le. methi))then
                    vr = (x(j)*vx(j) + y(j)*vy(j))/rad
                    vp = (-y(j)*vx(j) + x(j)*vy(j))/rad
                    k = k + 1
                    ur(k) = vr
                    up(k) = vp
                end if
            end if
        end if
    end do
    print*,"Used ",k," particles."


c plotting routines
    call PGWINDOW( vmin, vmax, umin, umax )
    do j = 1, k
        call PGPT1(up(j),ur(j),1)
    end do
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
c position, with the x-axis = 0 and going clockwise round from there
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