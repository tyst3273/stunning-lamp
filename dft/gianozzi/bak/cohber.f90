!---------------------------------------------------------------
program cohber
  !---------------------------------------------------------------
  !
  !     Band structure of semiconductors in the zincblende and
  !     diamond structure - Cohen and Bergstresser, PRB 141, 789 (1966)
  !     Expansion on a plane-wave basis set and diagonalization
  !     Units: hbar^2/2m = 1
  !     Requires lapack dsyev

  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi=3.14159265358979_dp, tpi=2.0_dp*pi
  !
  !    a  = lattice parameter (a.u.)
  !    VxN= pseudopotential form factors, in Ry, as in the paper
  !         x = symmetric term, a = antisymmetric term
  !         The following numbers are good for Silicon
  !
  real(dp), parameter :: a=10.26, vs3=-0.21, vs8=0.04, vs11=0.08, &
                                  va3=-0.00, va4=0.00, va11=0.00
  !    the two atoms are at -tau and +tau (in units of a) 
  real(dp), parameter ::tau1=0.125, tau2=0.125, tau3=0.125
  !                                  
  integer :: n, npw, npwx, nk
  real(dp) :: ecut, kg2, g2, vag, vsg, kg0(3), gij(3), h1(3), h2(3), h3(3)
  real(dp), allocatable :: k(:,:), kg(:,:), e(:), h(:,:), work (:)
  integer :: i,j,m, nmax, n1, n2, n3, lwork, info
  !
  !     Plane waves basis set: G_n=n*2pi/L, \hbar^2/2m*G^2 < Ecut
  !
  write (*,"('Cutoff for plane waves: ecut (Ry) > ')", advance='no')
  read (*,*) ecut
  if ( ecut <= 0.0_dp) stop ' wrong cutoff '
  !
  !     Number and list of k-vectors
  !
  write (*,"('Number of k-vectors > ')",advance='no')
  read (*,*) nk
  if ( nk <= 0 ) stop ' wrong input parameter '
  allocate ( k(3,nk) )
  write (*,"('k (in 2pi/a units) > ')")
  do i=1, nk
     read (*,*, end=10, err=10)  k(1,i), k(2, i), k(3,i)
  end do
  !
  !     Basis vectors for reciprocal lattice in 2pi/a units
  !
  h1(1) = 1.0_dp;   h1(2) = 1.0_dp; h1(3) =-1.0_dp
  h2(1) = 1.0_dp;   h2(2) =-1.0_dp; h2(3) = 1.0_dp
  h3(1) =-1.0_dp;   h3(2) = 1.0_dp; h3(3) = 1.0_dp
  !
  !     Loop on k-vectors
  !
  !  open (7,file='bands.out',status='unknown',form='formatted')
  !
  do n=1, nk
     !
     !    Count plane waves such that (\hbar^2/2m)(k+G)^2 < Ecut
     !    nmax is an estimate of max index useful in the generation
     !    of PWs as G(n1,n2,n3) =  n1*h1 + n2*h2 + n3*h3
     !
     nmax = nint ( sqrt (ecut) / (tpi/a * sqrt(3.0_dp) ) + 0.5 ) + 1
     npw = 0
     do n1 = -nmax, nmax
        do n2 = -nmax, nmax
           do n3 = -nmax, nmax
              ! this is k+G in 2pi/a units
              kg0(:) = k(:,n) + n1*h1(:) + n2*h2(:) + n3*h3(:)
              kg2 = (tpi/a)**2 * ( kg0(1)**2 + kg0(2)**2 + kg0(3)**2 )
              if ( kg2 <= ecut ) npw = npw+1
           end do
        end do
     end do
     print *, 'Number of plane waves=',npw
     !
     allocate (kg(3,npw), e(npw), work(3*npw), h(npw,npw) )
     !
     !    now generate PWs (beware: in 2pi/a units)
     !
     i = 0
     do n1 = -nmax, nmax
        do n2 = -nmax, nmax
           do n3 = -nmax, nmax
              kg0(:) = k(:,n) + n1*h1(:)+n2*h2(:)+n3*h3(:)
              kg2 = (tpi/a)**2 * ( kg0(1)**2+ kg0(2)**2+ kg0(3)**2 )
              if ( kg2 <= ecut ) then
                 i = i + 1
                 kg(:,i) = kg0(:)
              end if
           end do
        end do
     end do
     if ( i /= npw ) stop ' Some PWs missing in action'
     !       cleanup
     h(:,:) = 0.0_dp
     !
     !       Assign values of the matrix elements of the hamiltonian 
     !       on the plane wave basis
     !
     do i=1,npw
        do j=1,npw
           gij(:) = kg(:,i) - kg(:,j)
           g2 = gij(1)**2 + gij(2)**2 + gij(3)**2
           if ( abs (g2-3.0_dp) < 1.0d-6 ) then
              vsg = vs3
              vag = va3
           else if ( abs (g2-4.0_dp) < 1.0d-6 ) then
              vsg = 0.0_dp
              vag = va4
           else if ( abs (g2-8.0_dp) < 1.0d-6 ) then
              vsg = vs8
              vag = 0.0_dp
           else if ( abs (g2-11.0_dp) < 1.0d-6 ) then
              vsg = vs11
              vag = va11
           else
              vsg =0.0_dp
              vag =0.0_dp
           end if
           if ( i ==j ) then
              h(i,j) =  (tpi/a)**2 * (kg(1,i)**2+kg(2,i)**2+kg(3,i)**2) + vsg
           else
              h(i,j) = vsg * cos(tpi*(gij(1)*tau1+gij(2)*tau2+gij(3)*tau3)) 
                     !  vag * sin(tpi*(gij(1)*tau1+gij(2)*tau2+gij(3)*tau3))
           end if
           !print  '(2i4,f12.6)', i,j, h(i,j)
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
     write (*,'("k=",3f10.4)') k(:,n)
     write (*,'(4f12.4)') e(1:8)*13.6058
     !
     !       Write to output file the band dispersion e(k)
     !
     deallocate ( h, work, e, kg)
     !
  end do
  close(7)
  !
  stop
10 stop ' Error while reading k-vectors'

end program cohber
