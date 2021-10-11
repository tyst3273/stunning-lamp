!---------------------------------------------------------------
program helium_hf
!---------------------------------------------------------------
  !
  ! Solve the Hartree-Fock equations for He in the ground state
  ! using self-consistency on the potential and the Numerov algorithm 
  ! to integrate the radial Schroedinger equation under an effective
  ! potential (assumed to have spherical symmetry)
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi
  logical :: is_converged
  integer :: max_scf_iter=100
  integer :: mesh
  integer :: n, l, i, iter
  real (dp) :: zeta=2.0_dp, zmesh, rmax, xmin, dx, del, beta, tol
  real (dp) :: e, ehx, evion, ekin, etot, etot1
  real (dp), allocatable :: r(:), sqr(:), r2(:), y(:), vpot(:), &
       rho(:), vscf(:), vhx(:), deltav(:)
  !
  ! read atomic charge and other parameters
  !
  write(*,*) 
  write(*,'(" Hartree-Fock calculation for Helium-like atom")')
  write(*,*) 
  ! write (*,'(a,f10.6)') 'Atomic Charge = ', zeta
  write (*,'(" Mixing parameter beta [0.0-1.0] > ")',advance='no')
  read(*,*) beta
  if ( beta < 0.0_dp .or. beta > 1.0_dp ) stop 'beta not in [0.0-1.0]'
  write (*,'(" SCF accuracy > ")',advance='no')
  read (*,*) tol
  if ( tol < 0.0_dp ) stop 'tol should be positive'
  !
  ! initialize logarithmic mesh
  !
  zmesh=  zeta 
  rmax = 100.0_dp
  xmin = -8.0_dp
  dx   =  0.01_dp
  !
  mesh = (log(zmesh*rmax) - xmin)/dx
  !
  allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh),&
       rho(0:mesh), vscf(0:mesh), vhx(0:mesh), deltav(0:mesh) )
  !
  call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
  !
  ! initialize the potential
  !
  call init_pot ( zeta, r, mesh, vpot )
  vscf = vpot
  !
  ! The GS configuration for helium-like atoms is (1s)**2
  !
  n = 1
  l = 0
  !
  ! SCF cycle
  !
  do iter = 1, max_scf_iter
     write(*,'("####################################################")')
     write(*,'(" SCF iteration # ",i3)') iter
     !
     ! solve the schroedinger equation in radial coordinates by Numerov method
     !
     call solve_sheq ( n, l, e, mesh, dx, r, sqr, r2, vscf, zeta, y )
     !
     ! calculate the charge density from the wfc
     !
     call rho_of_r ( mesh, r, r2, y, rho )
     !
     ! calculate the Hartree + Exchange potential and energy
     !
     call v_of_rho ( mesh, dx, r, r2, rho, vhx, ehx )
     !
     ! calculate the kinetic energy and the energy in the external potential
     ! 
     evion = 0.0_dp
     ekin = 2.0_dp * e
     do i=0,mesh
        evion = evion + vpot(i) * rho(i) * fpi * r2(i) * r(i) * dx
        ekin  = ekin - vscf(i) * rho(i) * fpi * r2(i) * r(i) * dx
     end do
     !
     deltav = vpot + vhx -vscf
     vscf = vscf + beta * deltav
     !
     del = 0.0_dp
     do i=0,mesh
        del = del + deltav(i) * rho(i) * fpi * r2(i) * r(i) * dx
     end do
     !
     ! write out the eigenvalue energy to be compared with the external potential
     !
     etot = 2 * e - ehx + del
     etot1= ekin + evion + ehx
     write (*,'(" eigenvalue   = ",f15.6)') e
     write (*,'(" Eigenvalue energy       ",f15.6)' )  2.0_dp * e 
     write (*,'(" Kinetic energy          ",f15.6)' )  ekin
     write (*,'(" External pot. energy    ",f15.6)' )  evion
     write (*,'(" Hartree+Exchange energy ",f15.6)' )  ehx
     write (*,'(" Variational correction  ",e15.6)' )  del
     write (*,'(" Total energy            ",2f15.6)')  etot, etot1
     write (*,'(" Virial check            ",f15.6)' ) -(evion+ehx)/ekin
     !
     ! The variational correction "del" is used to check for convergence
     !
     is_converged = ( abs(del) < tol )
     if ( is_converged )  exit
     !
  end do
  !
  write (*,*) 
  if ( is_converged ) then
     write (*,'(" SCF Convergence has been achieved")')
  else
     write (*,'(" WARNING: SCF Convergence not achieved after ",i4," iterations")') max_scf_iter
  end if
  write (*,*) 
  write (*,'(" compute additional single-particle states")')
  write (*,*) 
  !
  ! write potentials to file pot.out
  !
  open (7,file='pot.out',status='unknown',form='formatted')
  write(7,"('# ',12x,'r',20x,'Vcoul(r)',18x,'Vhar(r)',14x,'Vscf(r)')")
  do i =0,mesh
     write (7,*) r(i),vpot(i),vhx(i),vscf(i)
  end do
  close(7)
  !
  ! open wfc file
  !
  open (7,file='wfc.out',status='unknown',form='formatted')
  write(7,"('# ',12x,'r',24x,'R(r)',18x,'X(r)=rR(r)',18x,'Vscf(r)')")
  !
  do
     !
     ! read principal and angular quantum numbers
     !
     write (*,'(" n, l >" )', advance='no') 
     read (*,*) n, l
     if ( n < 1 ) stop
     if ( n < l+1 ) then
        write(*,*) 'error in main: n < l+1 -> wrong number of nodes '
        cycle
     end if
     call solve_sheq ( n, l, e, mesh, dx, r, sqr, r2, vscf, zeta, y )
     write (*,'(" eigenvalue ",2i2," = ",f15.6)') n,l,e
     do i=0,mesh
        write (7,*) r(i),y(i)/sqr(i), y(i)*sqr(i), vscf(i)
     end do
     write (7,'(/)')
  end do
  
end program helium_hf
!
!--------------------------------------------------------------------
subroutine rho_of_r ( mesh, r, r2, y, rho )
  !--------------------------------------------------------------------
  !
  ! compute the charge density of an He-like atom
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi, &
       nelec=2.0_dp
  
  integer, intent(in) :: mesh
  real (dp), intent(in) :: r(0:mesh), r2(0:mesh), y(0:mesh)
  real (dp), intent(out):: rho(0:mesh)
  
  integer :: i
  
  rho = nelec * y**2 * r / (fpi*r2)
  
  return
end subroutine rho_of_r
!
!--------------------------------------------------------------------
subroutine v_of_rho ( mesh, dx, r, r2, rho, vhx, ehx )
  !--------------------------------------------------------------------
  !
  ! compute the Hartree + Exchange potential for an He-like atom
  ! Exchange cancels exactly half of the Hartree result for both 
  ! potential and energy
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi
  real (dp), parameter :: e2=2.0_dp ! assuming atomic Ry units
  integer, intent(in) ::mesh
  real (dp), intent (in) :: dx, r(0:mesh), r2(0:mesh), rho(0:mesh)
  real (dp), intent (out):: vhx(0:mesh), ehx
  
  integer :: i
  real (dp) :: charge
  !
  ! calculate the Hartree potential and energy by integrating the 
  ! electric field generated by the electronic charge. 
  ! This is done in 2 steps
  !
  ! 1) calculate the charge inside a sphere and fill vhx with the 
  !    electric field generated by this charge
  !
  charge = 0.0_dp
  do i=0,mesh
     charge = charge + rho(i) * fpi * r2(i) * r(i) * dx
     vhx(i) = e2*charge/r2(i)
  end do
  !
  !   (the total charge is written in output as a check)
  ! 
  write (*,'(" Check: total charge = ",f15.6)') charge
  !
  ! 2) integrate the electric field from +\infty to r to get the potential
  !    and integrate V_{Hartree}*rho to get the energy
  !
  ehx = 0.d0
  vhx(mesh) = e2 * charge / r(mesh)
  do i = mesh-1, 0, -1
     vhx(i) = vhx(i+1) + vhx(i) * r(i) * dx
     ehx = ehx + vhx(i) * rho(i) * fpi * r2(i) * r(i) * dx
  end do
  ehx = ehx/2.0_dp
  !
  ! Exchange cancels exactly half of the Hartree result for both 
  ! potential and energy
  !
  ehx = 0.5_dp * ehx 
  vhx = 0.5_dp * vhx
  !
  return
end subroutine v_of_rho
!
!--------------------------------------------------------------------
subroutine do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: mesh
  real (dp), intent (in) :: zmesh, xmin, dx
  real (dp), intent (out) :: r(0:mesh), sqr(0:mesh), r2(0:mesh)
  !
  integer :: i
  real(dp) :: x
  !
  do i=0,mesh
     x = xmin + dx * i
     r(i)  = exp(x)/zmesh
     sqr(i)= sqrt(r(i))
     r2(i) = r(i) * r(i)
  end do
  write(*,'(/" radial grid information:")')
  write(*,'(" dx =",f10.6,", xmin =",f10.6,", zmesh =",f10.6)') &
       dx,xmin,zmesh
  write(*,'(" mesh =",i6,", r(0) =",f10.6,", r(mesh) =",f10.6)') &
       mesh, r(0), r(mesh)
  write(*,*) 
  !
  return
end subroutine do_mesh
!
!--------------------------------------------------------------------
subroutine init_pot ( zeta, r, mesh, vpot )
!--------------------------------------------------------------------
  !
  ! initialize potential
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: mesh
  real (dp), intent(in) :: zeta, r(0:mesh)
  real (dp), intent(out):: vpot(0:mesh)
  integer :: i

  open (7,file='pot.out',status='unknown',form='formatted')
  write(7,'("#       r             V(r)")')
  do i = 0,mesh
     vpot (i) = - 2.0_dp*zeta/r(i)
     write (7,*) r(i),vpot(i)
  end do
  close(7)
  
  return
end subroutine init_pot
!---------------------------------------------------------------------
subroutine solve_sheq ( n, l, e, mesh, dx, r, sqr, r2, vpot, zeta, y )
  !---------------------------------------------------------------------
  !
  ! solve the Schroedinger equation in radial coordinates on a 
  ! logarithmic grid by Numerov method - atomic (Ry) units
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200), maxiter=100
  real(dp), parameter :: eps=1.0D-10
  integer, intent(in) :: mesh, n,l
  real(dp), intent(in) :: dx, r(0:mesh), sqr(0:mesh), r2(0:mesh), & 
       vpot(0:mesh), zeta
  real(dp), intent(out) :: e, y(0:mesh)
  
  integer :: i, j, icl, nodes, ncross, kkk
  real (dp) :: ddx12, sqlhf, x2l2, ycusp, dfcusp, fac, norm, eup, elw, de
  real (dp), allocatable :: f(:)
  
  allocate ( f(0:mesh) )
  ddx12=dx*dx/12.0_dp
  sqlhf = (l+0.5_dp)**2
  x2l2 = 2*l+2
  !
  ! set (very rough) initial lower and upper bounds to the eigenvalue
  !
  eup = vpot(mesh)
  elw=minval (sqlhf/r2(:)+vpot(:))
  if (eup-elw < eps) then
     write (*,*) elw, eup
     write (*,*) 'solve_sheq: lower and upper bounds are equal'
     stop
  end if
  e = 0.5_dp * (elw + eup)

  do kkk = 1, maxiter
     !
     ! set up the f-function and determine the position of its last
     ! change of sign
     ! f < 0 (approximately) means classically allowed   region
     ! f > 0         "         "        "      forbidden   "
     !
     icl = -1
     f(0) = ddx12 *( sqlhf + r2(0) * (vpot(0)-e) )
     do i=1,mesh
        f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
        !
        ! beware: if f(i) is exactly zero the change of sign is not observed
        ! the following line is a trick to prevent missing a change of sign 
        ! in this unlikely but not impossible case:
        !
        if( f(i) == 0.0_dp ) f(i)=1.d-20
        if( f(i) /= sign(f(i),f(i-1)) ) icl=i
     end do

     if ( icl < 0 .or. icl >= mesh-2 ) then
        !
        ! classical turning point not found or too far away
        ! no panic: it may follow from a bad choice of eup in
        ! the first iterations. Update e and eup and re-try
        eup = e
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     !
     ! f function as required by numerov method
     !
     f = 1.0_dp - f
     y = 0
     !
     ! determination of the wave-function in the first two points 
     ! (asymptotic behaviour - second term depends upon the potential)
     !
     nodes = n - l - 1
     y(0) = r(0)**(l+1) *(1.0_dp - 2.0_dp*zeta*r(0)/x2l2) / sqr(0)
     y(1) = r(1)**(l+1) *(1.0_dp - 2.0_dp*zeta*r(1)/x2l2) / sqr(1)
     !
     ! outward integration, count number of crossings
     !
     ncross=0
     do i =1, icl-1
        y(i+1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
        if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
     end do
     fac = y(icl) 
     !
     ! check number of crossings
     !
     if ( ncross /= nodes ) then
        if ( ncross > nodes ) then
           eup = e
        else
           elw = e
        end if
        e = 0.5_dp * (eup+elw)
        cycle
     end if
     !
     ! determination of the wave-function in the last two points 
     ! assuming y(mesh+1) = 0 and y(mesh) = dx
     !
     y(mesh) = dx
     y(mesh-1) = (12.0_dp-10.0_dp*f(mesh))*y(mesh)/f(mesh-1) 
     !
     ! inward integration 
     !
     do i = mesh-1,icl+1,-1
        y(i-1)=((12.0_dp-10.0_dp*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)
        if (y(i-1) > 1.0d10) then
           do j=mesh,i-1,-1
              y(j) = y(j)/y(i-1)
           end do
        end if
     end do
     !
     ! rescale function to match at the classical turning point (icl)
     !
     fac = fac/y(icl)
     y (icl:) = y(icl:)*fac
     !
     ! normalization - note the change of variable:
     !  \int f(r)dr => \sum_i f_i r_i Delta x
     !
     norm = 0.d0
     do i=0,mesh
        norm = norm + y(i)*y(i) * r2(i) * dx
     end do
     norm = sqrt(norm)
     y = y / norm
     !
     ! find the value of the cusp at the matching point (icl)
     !
     i = icl
     ycusp = (y(i-1)*f(i-1)+f(i+1)*y(i+1)+10.0_dp*f(i)*y(i)) / 12.0_dp
     dfcusp = f(i)*(y(i)/ycusp - 1.0_dp)
     !
     ! eigenvalue update using perturbation theory
     !
     de = dfcusp/ddx12 * ycusp*ycusp * dx 
     if (de > 0.0_dp) elw = e
     if (de < 0.0_dp) eup = e
     !
     ! prevent e to go out of bounds, i.e. e > eup or e < elw 
     ! (might happen far from convergence)
     !
     e = max( min (e+de,eup),elw)
     !
     ! convergence check
     !
     if (abs(de) < eps) exit
     !
  end do
  !
  ! was convergence achieved ?
  !
  if ( abs(de) > 1.d-10 ) then
     if ( ncross /= nodes ) then
        write (*,*) e, elw, eup, ncross, nodes, icl
     else
        write (*,*) e, de
     end if
     write (*,*) ' error in solve_sheq: too many iterations'
     stop
  else 
  ! ---- convergence has been achieved -----
     write (*,'(" convergence achieved at iter #",i3," de = ", e10.4)') kkk,de
  end if
  return
end subroutine solve_sheq
