!---------------------------------------------------------------
program h2_hf_gauss
  !---------------------------------------------------------------
  !
  !    Hartree-Fock ground-state solution for the H_2 molecule
  !    S=0 state - equivalent to Hartree approximation
  !    expansion on a gaussian basis set and self-consistency
  !    with diagonalization of the Fock matrix
  !    Ry atomic units: \hbar^2/2m = 1
  !
  !    Requires lapack and blas
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, tpi=2.0_dp*pi, &
       e2=2.0_dp, zeta=1.0_dp
  logical :: do_scf =.true.
  ! .true.  for H2 calculation, .false. for H2+ calculation
  integer :: ia,ib,ic,id, i,j,n, n_alpha, ngaus, iter, maxter=100
  real (dp) :: aa, ab, dab, rab, d2, dmin, dmax, dd, eold, enew, enuc, e0
  real (dp), allocatable ::  alpha(:), d(:), e(:), c(:)
  real (dp), allocatable ::  h(:,:), f(:,:), s(:,:), v(:,:)
  real (dp), allocatable ::  q(:,:,:,:)
  ! real (dp), external :: erf
  character (len=80) :: filename
  integer :: iunit, ios
  !
  !       Input data
  !
  write (*,'(" Parameters of the Gaussians from file > ")',advance='no')
  read (*,'(a)') filename
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
  allocate (alpha (2*n_alpha))
  do i=1,n_alpha
     if ( iunit == 5 ) then
        write (*,'(" Gaussian # ",i3,",  coefficient > ")',advance='no') i
        read (iunit,*) alpha(i)
     else
        read (iunit,*) alpha(i)
        write (*,'(" Gaussian # ",i3,",   coefficient > ",f10.6)') i, alpha(i)
     end if
     if ( alpha(i) <= 0.0_dp ) stop 'wrong coefficient'
  end do
  if ( iunit /= 5 ) close (unit = iunit, status='keep')
  !
  open (unit = 7, file='h2.out', form='formatted', status='unknown')
  !
  write (*,'(" H2 distance (a.u.): Min, Max, step >> ")',advance='no')
  read (5,*) dmin, dmax, dd
  !
  ngaus = 2*n_alpha
  allocate ( d(ngaus) )
  allocate ( s(ngaus, ngaus), h(ngaus, ngaus), f (ngaus, ngaus), &
       v(ngaus, ngaus), c(ngaus), e (ngaus)  )
  allocate ( q(ngaus, ngaus, ngaus, ngaus) )
  !
  ! Starting solution (very lousy guess) for the first calculation
  ! Following calculations start from the solution at previous distance
  !
  do ia=1,ngaus
     c(ia) = 0.d0
  end do
  c(1) = 1.d0
  !
  ! loop over H-H distances (d2)
  !
  n = nint((dmax-dmin)/dd)
  do j = 0,n
     !
     d2 = dmin + j*dd
     !
     ! Basis set: b_i(r) = exp (-alpha_i (r-R_i)^2)
     ! (first n_alpha centered at -R/2, second n_alpha centered at +R/2)
     !
     do i=1, n_alpha
        alpha(n_alpha+i) = alpha(i)
        d(i        ) = -d2/2
        d(i+n_alpha) = +d2/2
     end do
     !
     !       Assign values of overlap integrals S and of matrix elements 
     !       of the one-electron hamiltonian H on the gaussian basis:
     !
     do ib=1,ngaus
        do ia=1,ngaus
           aa = alpha(ia)+alpha(ib)
           ab = alpha(ia)*alpha(ib)/aa
           dab= (d(ia)-d(ib))**2
           ! overlap
           s(ia,ib) = (pi/aa)**1.5_dp*exp(-ab*dab)
           ! kinetic
           h(ia,ib) = s(ia,ib)*ab*(6.0_dp-4.0_dp*ab*dab)
           ! add coulombian: beware the limit R_i-R_j=0
           do i=-1,1,2
              rab= abs ( (alpha(ia)*d(ia)+alpha(ib)*d(ib))/aa - i*d2/2 )
              if ( rab < 1.0E-8 ) then
                 h(ia,ib) = h(ia,ib) - e2*zeta*tpi/aa
              else
                 h(ia,ib) = h(ia,ib) - e2*zeta*s(ia,ib) * &
                      erf(sqrt(aa)*rab)/rab
              end if
           end do
        end do
     end do
     !
     ! Fill the Q matrix with the matrix elements of the e-e interaction:
     ! q(i,j,k,l) = \int b_i(r) b_j(r') 1/|r-r'| b_k(r) b_l(r') dr dr'
     ! beware the ordering of indices!
     !
     do id=1,ngaus
        do ic=1,ngaus
           do ib=1,ngaus
              do ia=1,ngaus
                 aa = alpha(ia)+alpha(ic)
                 dab= (alpha(ia)*d(ia)+alpha(ic)*d(ic))/aa
                 ab = alpha(ib)+alpha(id)
                 rab= (alpha(ib)*d(ib)+alpha(id)*d(id))/ab
                 if ( abs(dab-rab) < 1.0E-8) then
                    q(ia,ib,ic,id) = e2 * s(ia,ic) * s(ib,id) * &
                         2.0_dp / sqrt(pi) * sqrt(aa*ab/(aa+ab))
                 else
                    q(ia,ib,ic,id) = e2 * s(ia,ic) * s(ib,id) / abs(dab-rab) *&
                         erf ( sqrt(aa*ab/(aa+ab))* abs(dab-rab) )
                 end if
              end do
           end do
        end do
     end do
     !
     !       Self-consistency iteration
     !
     enew = 0.d0
     !
     scf: do iter=1,maxter
        !
        !       Fill the Fock matrix
        !
        do ia=1,ngaus
           do ib=1,ngaus
              f(ia,ib) = h(ia,ib)
              if ( do_scf ) then
                 do id=1,ngaus
                    do ic=1,ngaus
                       ! This is for Hartree ...
                       ! f(ia,ib) = f(ia,ib) + q(ia,ic,ib,id)*c(ic)*c(id)
                       ! ...and this is for Hartree-Fock ...
                       f(ia,ib) = f(ia,ib) + (2.0_dp*q(ia,ic,ib,id) - &
                                             q(ia,ib,ic,id) )*c(ic)*c(id)
                    end do
                 end do
              end if
           end do
        end do
        !
        !       Solution [expansion coefficients are stored into v(j,i)
        !                 j=basis function index, i= eigenvalue index]
        !
        call diag ( ngaus, ngaus, f, s,  e, v )
        !
        c(:) =  v(:,1)
        !
        eold = enew
        enew = 0.d0
        do ia=1,ngaus
           do ib=1,ngaus
              enew = enew + 2.0_dp*h(ia,ib)*c(ia)*c(ib)
              do ic=1,ngaus
                 do id=1,ngaus
                    ! This is for Hartree ...
                    !enew = enew + q(ia,ic,ib,id)*c(ia)*c(ib)*c(ic)*c(id)
                    ! ...and this is for Hartree-Fock ...
                    enew = enew + c(ia)*c(ib)*c(ic)*c(id) * &
                           ( 2.0_dp*q(ia,ic,ib,id)-q(ia,ib,ic,id) )
                 end do
              end do
           end do
        end do
        !write (*,'(" Iteration # ",i3,":  HF eigenvalue, energy = ",2f12.6)')&
        !     iter, e(1), enew
        !
        if ( abs (enew-eold) < 1.d-8 ) then
           !write (*,'(/" Convergence achieved, stopping")')
           if (j==0) write (*,'(/8x,"d      electron,   nuclear,   total (Ry) and binding (eV) energy")')
           if (j==0) write (7,'(/8x,"#  d   electron,   nuclear,   total (Ry) and binding (eV) energy")')
           if ( d2 > 1.0e-8) then
              enuc =  zeta**2*e2/d2
           else
              enuc = 0.0_dp
           end if
           if ( do_scf) then
              e0   =-2.0_dp
           else
              e0   =-1.0_dp
              enew = e(1)
           end if
           write (*,'(5f12.6)') d2,enew,enuc,enew+enuc,(enew+enuc-e0)*13.6058
           write (7,'(5f12.6)') d2,enew,enuc,enew+enuc,(enew+enuc-e0)*13.6058
           exit scf
        end if
        !
     end do scf
     !
  end do
  !
  close(unit=7)
  deallocate ( q, e, c, v, f, h, s, d, alpha )
  !
  stop 
  !
end program h2_hf_gauss
!
!-----------------------------------------------------------------------
subroutine diag ( n, ldh, h, s, e, v )
!     =================
!
!       Finds eigenvalues and eigenvectors of the generalized problem
!       Hv=eSv, where H=hermitian matrix, S=overlap matrix
!       On return S and H are unchanged.
!       Uses level-2 blas dgemm for matrix-matrix multiplication
!       Uses lapack dsyev for matrix diagonalization, checks for
!       too small eigenvalues of S, removes corresponding eigenvectors
!
!-----------------------------------------------------------------------
!
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)

  integer, intent(in) :: n, &! dimension of the matrix to be diagonalized
                         ldh ! leading dimension of h,s and v, as declared
                             ! in the calling pgm unit
  real (dp), intent(in) :: h(ldh,n), & ! matrix to be diagonalized
                           s(ldh,n)    ! overlap matrix
  real (dp), intent(out):: e(n), &     ! eigenvalues
                           v(ldh,n)    ! eigenvectors (column-wise)
  !
  ! local variables
  integer lwork, info, i, j, nn
  !
  real (dp), parameter :: small=1.d-10
  real (dp), allocatable :: work(:), aux(:,:), h1(:,:)
  !
  info = 0
  lwork=3*n
  allocate (work(lwork), aux(ldh,n))
  !
  !       Copy S into an auxiliary matrix because dsyev destroys the matrix
  !
  aux = s
  !
  !       Diagonalize S
  !
  call dsyev ( 'V', 'U', n, aux, ldh, e, work, lwork, info )
  if (info /= 0) stop 'S-matrix diagonalization failed '
  !
  !   Keep only linearly independent combinations (within a given threshold)
  !   store into matrix "aux" the eigenvectors of S divided by the squares
  !   of the eigenvectors
  !
  nn = 0
  do i=1,n
     if (e(i) > small) then
        nn = nn + 1
        aux(:,nn) = aux(:,i) / sqrt(e(i))
     end if
  end do
  if ( nn < n) write(*,*) " # of linearly independent vectors =", nn, n
  !
  !   Transform H using the "aux" matrix
  !
  !   V(i,j) = \sum_{k=1}^{n} H(i,k) aux(k,j),  i=1,n, j=1,nn
  !
  call dgemm ( 'N', 'N', n, nn, n, 1.0_dp, h, ldh, aux, ldh, 0.0_dp, v, ldh )
  !
  !    h1(i,j) = \sum_{k=1}^{n} aux(k,i) v(k,j),  i=1,nn, j=1,nn
  !    H' = transpose(aux)*H*aux
  !
  allocate (h1(nn,nn) )
  call dgemm ( 'T', 'N', nn, nn, n, 1.0_dp, aux, ldh, v, ldh, 0.0_dp, h1, nn )
  !
  !    Diagonalize transformed H
  !
  info = 0
  call dsyev ('V', 'U', nn, h1, nn, e, work, lwork, info)
  if (info /= 0) stop 'H-matrix diagonalization failed '
  !
  !    Back-transform eigenvectors
  !
  call dgemm ( 'N', 'N', n, nn, nn, 1.0_dp, aux, ldh, h1, nn, 0.0_dp, v, ldh )
  !
  deallocate (h1, aux, work)
  !
end subroutine diag

