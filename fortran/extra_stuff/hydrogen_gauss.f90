!
!---------------------------------------------------------------
program hydrogen_gauss
  !---------------------------------------------------------------
  !
  !       Solution of Hydrogenic atoms via expansion on a gaussian
  !       basis set and diagonalization (using linear algebra methods)
  !
  !       Atomic (Ry) units. Requires lapack dsyev
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  !
  real (dp), parameter :: e2=2.0_dp, pi=3.14159265358979_dp, tpi=2.0_dp*pi
  integer :: ia,ib,n_alpha
  real (dp) :: zeta, aa
  real (dp), allocatable :: alpha(:), e(:), s(:,:), h(:,:), v(:,:)
  real (dp) :: r,f,prob,q,dr
  integer :: iunit, ios, i,j, nrx
  character(len=80) :: filename
  !
  !       Input data
  !
  write (*,'(" Atomic Charge Z > ")',advance='no') 
  read (*,*) zeta
  ! write (*,'(" Atomic charge = ",f10.6)') zeta
  if ( zeta < 1.0_dp) stop 'zeta should be >= 1'
  write (*,'(" Parameters of the Gaussians from file > ")',advance='no')
  read (*,'(A)') filename
  if ( trim(filename) == '') then
     !  read from terminal
     iunit = 5
  else
     write (*,*)
     !  read from specified file
     iunit = 10
     open (unit=iunit, file=trim(filename), status='old', &
           form='formatted', iostat=ios )
     if ( ios /= 0 ) stop 'file not found'
  end if 
  if ( iunit == 5 ) write (*,'(" Number of Gaussians > ")',advance='no')
  read (iunit,*) n_alpha
  allocate (alpha (n_alpha))
  do i=1,n_alpha
     if ( iunit == 5 ) then
        write (*,'(" Gaussian # ",i3," coefficient > ")',advance='no') i
        read (iunit,*) alpha(i)
     else
        read (iunit,*) alpha(i)
        write (*,'(" Gaussian # ",i3," coefficient > ",f10.6)') i, alpha(i)
     end if
     if ( alpha(i) <= 0.0_dp ) stop 'wrong coefficient'
  end do
  if ( iunit /= 5 ) close (unit = iunit, status='keep')
  !
  !       Allocate arrays
  !
  allocate (s(n_alpha, n_alpha), h(n_alpha,n_alpha), v(n_alpha,n_alpha), &
            e(n_alpha) )
  !
  !*******************************************************************
  !   S States
  !*******************************************************************
  !
  !       Assign values of overlap integrals and of matrix elements 
  !       of the hamiltonian on the S basis:
  !
  do ia=1,n_alpha
     do ib=1,n_alpha
        aa = alpha(ia)+alpha(ib)
        s(ia,ib) = (pi/aa)**1.5_dp
        h(ia,ib) = s(ia,ib)*6.0_dp*alpha(ia)*alpha(ib)/aa - e2*zeta*tpi/aa 
     end do
  end do
  !
  !       Solution [expansion coefficients are stored into v(j,i)
  !                 j=basis function index, i= eigenvalue index]
  !
  call diag ( n_alpha, n_alpha, h, s, e, v )
  !
  write (*,'(/" Lowest eigenvalues, S states:",3f12.8)') (e(ia),ia=1,3)
  !
  !
  !       Write to output file the coefficients of the lowest-energy state:
  !
  open (7,file='s-coeff.out',status='unknown',form='formatted')
  do j=1,n_alpha
     write(7,*) j, alpha(j), v(j,1)
  end do
  close(7)
  !
  !       Write to output file the lowest-energy state:
  !
  open (7,file='s-wfc.out',status='unknown',form='formatted')
  dr =0.1_dp/zeta
  nrx = 200
  q = 0.0_dp
  do i=0, nrx
     r = dr*i
     f = 0.d0
     do j=1,n_alpha
        f = f + exp(-alpha(j)*r*r)*v(j,1)
     end do
     prob = 4.0_dp*pi*(r*f)**2
     q   = q + prob*dr
     write(7,*) r, f, r*f, prob
  end do
  !       verify normalization (if desired):
  !      write (*,*) q
  close(7)
  
  !*******************************************************************
  !   P States 
  !*******************************************************************
  !
  !       Assign values of overlap integrals and of matrix elements 
  !       of the hamiltonian on the S basis:
  !
  do ia=1,n_alpha
     do ib=1,n_alpha
        aa = alpha(ia)+alpha(ib)
        s(ia,ib) = 0.5_dp/aa*(pi/aa)**1.5_dp
        h(ia,ib) = s(ia,ib)*10.0_dp*alpha(ia)*alpha(ib)/aa - &
                   e2*zeta*tpi/aa**2/3.0_dp
     end do
  end do
  !
  !       Solution [expansion coefficients are stored into v(j,i)
  !                 j=basis function index, i= eigenvalue index]
  !
  call diag ( n_alpha, n_alpha, h, s, e, v )
  !
  write (*,'( " Lowest eigenvalues, P states:",3f12.8)') (e(ia),ia=1,3)
  !
  !       Write to output file the coefficients of the lowest-energy state:
  !
  open (7,file='p-coeff.out',status='unknown',form='formatted')
  do j=1,n_alpha
     write(7,*) j, alpha(j), v(j,1)
  end do
  close(7)
  !
  !       Write to output file the lowest-energy state:
  !
  open (7,file='p-wfc.out',status='unknown',form='formatted')
  dr =0.1_dp/zeta
  nrx = 200
  q = 0.0_dp
  do i=0, nrx
     r = dr*i
     f = 0.0_dp
     do j=1,n_alpha
        f = f + r*exp(-alpha(j)*r**2)*v(j,1)
     end do
     prob = 4.0_dp/3.0_dp*pi*(r*f)**2
     q   = q + prob*dr
     write(7,*) r, f, r*f, prob
  end do
  !       verify normalization (if desired):
  !     write (*,*) q
  close(7)
  
end program hydrogen_gauss
