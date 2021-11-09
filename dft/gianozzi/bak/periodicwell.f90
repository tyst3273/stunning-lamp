!---------------------------------------------------------------
program periodicwell
  !---------------------------------------------------------------
  !
  !     Band structure of a 1d model periodic system (Kronig-Penney)
  !     Expansion on a plane-wave basis set and diagonalization
  !     Units: hbar^2/2m = 1
  !     Requires lapack dsyev

  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi=3.14159265358979_dp
  integer :: n, npw
  real(dp) :: v0, a, b, ecut, k
  real(dp), allocatable :: g(:), e(:), h(:,:), work (:)
  complex(dp) :: f
  integer :: i,j,m, lwork, info
  !
  !       Input data
  !
  !     Potential well: V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2
  !     Periodicity:    V(x+a)=V(x)
  !
  write (*,"('Parameters for potential well: V_0, a, b > ')",advance='no')
  read (*,*) v0, a, b
  if ( v0 <= 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp .or. a <= b ) &
     stop ' wrong input parameters '
  write (*,"('   V_0, a, b =',3f10.4)") v0, a, b
  !
  !     Plane waves basis set: G_n=n*2pi/a, \hbar^2/2m*G^2 < Ecut
  !
  write (*,"('Cutoff for plane waves: ecut > ')", advance='no')
  read (*,*) ecut
  if ( ecut <= 0.0_dp) stop ' wrong input parameter '
  !
  !     npw = number of plane waves with \hbar^2G^2/2m <= ecut
  !
  npw = nint ( sqrt ( ecut/(2.0_dp*pi/a)**2 ) )
  npw = 2*npw+1
  write (*,"('Ecut=',f8.2,'  # of PWs=',i4)") ecut,npw
  allocate (g(npw), e(npw), work(3*npw), h(npw,npw) )
  !
  !     Assign values of G_n: n=0,+1,-1,+2,-2, etc
  !
  g(1) = 0.0_dp
  do i=2,npw-1,2
     g(i  ) = (i/2)*2.0_dp*pi/a
     g(i+1) =-(i/2)*2.0_dp*pi/a
  end do
  !
  !     Loop on k-vectors: k runs from -pi/a to pi/a
  !
  open (7,file='bands.out',status='unknown',form='formatted')
  n = 20
  do m=-n,n
     !
     k = m*pi/n/a
     !       cleanup
     h(:,:) = 0.0_dp
     !
     !       Assign values of the matrix elements of the hamiltonian 
     !       on the plane wave basis
     !
     do i=1,npw
        do j=1,npw
           if ( i ==j ) then
              h(i,j) = (k+g(i))**2 - v0/a*b
           else
              h(i,j) = -v0/a * sin( (g(j)-g(i))*b/2.0_dp ) / (g(j)-g(i))*2.0_dp
           end if
           ! print  '(2i4,f12.6)', i,j, h(i,j)
        end do
     end do
     !
     !       Solution [expansion coefficients are stored into h(j,i)
     !                 j=basis function index, i= eigenvalue index]
     !
     lwork = 3*npw
     call dsyev ( 'V', 'U', npw, h, npw, e, work, lwork, info )
     if (info /= 0) stop 'H-matrix diagonalization failed '
     !
     write (*,"('k=',f10.4,'    Lowest eigenvalues:',3f12.6)") k,e(1),e(2),e(3)
     !
     !       Write to output file the band dispersion e(k)
     !
     write (7,'(4f12.6)') k, e(1), e(2), e(3)
  end do
  close(7)
  deallocate ( h, work, e, g)
  !
end program periodicwell
