!---------------------------------------------------------------
program ah
  !---------------------------------------------------------------
  !
  !     Self-Consistent calculation for Si - local pseudopotentials
  !        Joel A. Appelbaum and D. R. Hamann, PRB 8, 1777 (1973)
  !     Hartree potential + Slater exchange, plane-wave basis set
  !     Units: hbar^2/2m = 1, e2=2 (Ry a.u.)
  !
  !     Use Fast Fourier Transform for G <-> r conversion,
  !     conventional diagonalization of the H matrix
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi=3.14159265358979_dp, tpi=2.0_dp*pi
  ! a small quantity
  real(dp), parameter :: eps=1.0d-8
  ! Important variables
  !    a  = lattice parameter (a.u.)
  !    alpha = mixing parameter for scf
  !    threshold = convergence threshold for scf procedure
  real(dp) :: a=10.26, alpha = 0.5, threshold = 1.0d-6
  !    the two atoms are at -tau and +tau (in units of a) 
  real(dp) :: tau1=0.125, tau2=0.125, tau3=0.125
  !    nelec = number of electrons per unit cell
  !    nk    = number of k-points in the Brillouin Zone
  !    nbands= number of occupied bands
  !    max_iter = max number of scf iterations
  integer :: nelec=8, nk = 1, nbands = 4, max_iter = 50
  !    npwx= max number of plane waves across k-points
  !    ngm = number of G vectors with G^2 < 4*Ecut
  !    nm* = Max value of Miller indices
  !    nr* = real-space grid dimensions
  integer :: npwx, ngm, nm1, nm2, nm3, nr1, nr2, nr3
  !    npw  = number of plane waves per k-point
  !    mill = Miller indices of G vectors
  !    indg = index of G vector corresponding to Miller indices
  !    igk  = index of G vector in list of k+G such that (k+G)^2 < Ecut
  integer,  allocatable :: npw(:), mill(:,:), indg(:,:,:), igk(:,:)
  !    ecut = kinetic energy cutoff for plane waves, in Ry
  !    a1, a2, a3: direct, h1, h2, h3: reciprocal lattice vectors
  !    omega= volume of the unit cell
  real(dp) :: ecut, a1(3), a2(3), a3(3), h1(3), h2(3), h3(3), omega
  !    k    = k-points, in units 2\pi/a
  !    wk   = normalized weights of k-points
  !    e, h = eigenvalues, eigenvectors (overwritten to hamiltonian)
  !    g    = G vectors (in units 2pi/a) with G^2 < 4*Ecut
  !    g2   = G^2, in atomic units
  !    sg   = geometrical structure factor
  real(dp), allocatable :: k(:,:), wk(:), e(:,:), h(:,:), g(:,:), g2(:), sg(:)
  !    rhoin =  input charge density in a scf step, in real space
  !    rhoout= output charge density in a scf step, in real space
  !    vg    = exchange and coulomb potential, in G space
  real (dp), allocatable :: rhoin(:,:,:),  rhoout(:,:,:), vg(:)
  !    psi = valence bands (in plane waves). In general, complex,
  !          here, real because of inversion symmetry
  real (dp), allocatable :: psi(:,:)
  !
  ! Work variables
  !
  real(dp) :: vsg, gg, g0(3), kg2, kg0(3), r(3), drho2
  real(dp), external :: form_factor, mix_charge
  integer :: i, j, ik, jk, n, ng, nn, n1,n2,n3, iter
  !
  !     Plane waves basis set: G_n=n*2pi/L, \hbar^2/2m*G^2 < Ecut
  !
  write (*,"('Cutoff for plane waves: ecut (Ry) > ')", advance='no')
  read (*,*) ecut
  if ( ecut <= 0.0_dp) stop ' wrong cutoff '
  !
  !     Grid of k-vectors
  !
  !write (*,"('Number of k-vectors > ')", advance='no')
  !read (*,*) nk
  !if ( nk <= 0 ) stop ' wrong input parameter '
  allocate ( k(3,nk), wk(nk), npw(nk) )
  !write (*,"('k (in 2pi/a units) > ')")
  !do i=1, nk
  !   read (*,*, end=10, err=10)  k(1,i), k(2, i), k(3,i)
  !end do
  k(1,1) = 0.6223_dp ; k(2,1) = 0.2953_dp ; k(3,1) = 0.0_dp; wk(1)=1.0_dp
  !
  !     Basis vectors for fcc lattice in lattice parameter "a" units
  !
  a1(1) = 0.5_dp;   a1(2) = 0.5_dp; a1(3) = 0.0_dp
  a2(1) = 0.5_dp;   a2(2) = 0.0_dp; a2(3) = 0.5_dp
  a3(1) = 0.0_dp;   a3(2) = 0.5_dp; a3(3) = 0.5_dp
  !     In general: omega = a1 * ( a2 x a3 )
  omega = a**3 /4.0_dp 
  !
  !     Basis vectors for reciprocal lattice in 2pi/a units
  !
  h1(1) = 1.0_dp;   h1(2) = 1.0_dp; h1(3) =-1.0_dp
  h2(1) = 1.0_dp;   h2(2) =-1.0_dp; h2(3) = 1.0_dp
  h3(1) =-1.0_dp;   h3(2) = 1.0_dp; h3(3) = 1.0_dp
  !  
  !    Count G-vectors such that (\hbar^2/2m)G^2 < 4*Ecut
  !    nm1,nm2,nm3 are estimate of max index n1,n2,n3 used in the generation
  !    of G-vectors: G(n1,n2,n3) =  n1*h1 + n2*h2 + n3*h3. Since n1=G*a1 etc
  !    then abs(n1) <= |G||a1| etc
  !
  nm1 = nint ( sqrt(4.0_dp*ecut)/tpi*a*sqrt(a1(1)**2+a1(2)**2+a1(3)**2) + 0.5 )
  nm2 = nint ( sqrt(4.0_dp*ecut)/tpi*a*sqrt(a2(1)**2+a2(2)**2+a2(3)**2) + 0.5 )
  nm3 = nint ( sqrt(4.0_dp*ecut)/tpi*a*sqrt(a3(1)**2+a3(2)**2+a3(3)**2) + 0.5 )
  !
  ngm = 0
  do n1 = -nm1, nm1
     do n2 = -nm2, nm2
        do n3 = -nm3, nm3
           g0(:) = n1*h1(:) + n2*h2(:) + n3*h3(:)
           gg = (tpi/a)**2 * ( g0(1)**2 + g0(2)**2 + g0(3)**2 )
           if ( gg <= 4.0_dp*ecut ) ngm = ngm + 1
        end do
     end do
  end do
  write(*, '("Number of G-vectors=",i6)') ngm
  !
  !     Allocate and fill arrays for G-vectors: G, G^2, Miller indices
  !
  allocate ( g(3,ngm), g2(ngm), mill(3,ngm) ) 
  !
  ng = 0
  do n1 = -nm1, nm1
     do n2 = -nm2, nm2
        do n3 = -nm3, nm3
           g0(:) = n1*h1(:) + n2*h2(:) + n3*h3(:)
           gg = (tpi/a)**2 * ( g0(1)**2 + g0(2)**2 + g0(3)**2 )
           if ( gg <= 4.0_dp*ecut ) then
              ng = ng + 1
              g (:,ng) = g0(:)
              g2(ng) = gg
              mill(1,ng) = n1; mill(2,ng) = n2; mill(3,ng) = n3
           end if
        end do
     end do
  end do
  if ( ng /= ngm ) then
     print *, 'Mismatch in number of G-vectors!',ng,ngm
     stop
  end if
  !
  !   indg gives the index of a G vector from its Miller indices
  !
  nm1 = maxval ( mill(1,:) )
  nm2 = maxval ( mill(2,:) )
  nm3 = maxval ( mill(3,:) )
  allocate ( indg (-nm1:nm1, -nm2:nm2, -nm3:nm3) )
  indg (:,:,:) = 0
  do ng = 1, ngm
     indg ( mill(1, ng), mill(2,ng), mill(3,ng) ) = ng
  end do
  !
  !    Geometrical structure factors S(G)
  !
  allocate (sg(ngm))
  sg(:) = 0.0_dp
  do ng = 1, ngm
     sg(ng) = 2.0_dp*cos(tpi*(g(1,ng)*tau1 + g(2,ng)*tau2 + g(3,ng)*tau3))
  end do
  sg(:) = sg(:)/omega
  !
  !    For real-space quantities: dimensions of real-space grids
  !
  nr1 = 2*nm1 + 2
  nr2 = 2*nm2 + 2
  nr3 = 2*nm3 + 2
  write(*,'("Real-space grid nr1,nr2,nr3=",3i4)') nr1,nr2,nr3
  !
  do n=1, nk
     !
     !    Count plane waves such that (\hbar^2/2m)(k+G)^2 < Ecut
     !
     npw(n) = 0
     do ng = 1, ngm
        ! this is k+G in 2pi/a units
        kg0(:) = k(:,n) + g(:,ng)
        kg2 = (tpi/a)**2 * ( kg0(1)**2 + kg0(2)**2 + kg0(3)**2 )
        if ( kg2 <= ecut ) npw(n) = npw(n) + 1
     end do
     write(*, '("Number of plane waves for k-vector ",i3,"=",i4)') n,npw(n)
  end do
  !
  !     Max number of plane waves over all k-points
  !
  npwx = MAXVAL ( npw(:) )
  !
  !     Allocate and fill array of indices for k+G
  !
  allocate ( igk(npwx,nk) )
  do n=1, nk
     nn = 0
     do ng = 1, ngm
        ! this is k+G in 2pi/a units
        kg0(:) = k(:,n) + g(:,ng)
        kg2 = (tpi/a)**2 * ( kg0(1)**2 + kg0(2)**2 + kg0(3)**2 )
        if ( kg2 <= ecut ) then
           nn = nn + 1
           igk(nn,n) = ng  
        end if
     end do
     if ( nn /= npw(n)) then
        print *, 'Mismatch in number of plane waves!',n,nn,npw(n)
        stop
     end if
     !
  end do
  !
  allocate ( e(nbands,nk), h(npwx,npwx), psi(npwx,nbands) )
  allocate ( rhoin(nr1,nr2,nr3),  rhoout(nr1,nr2,nr3), vg(ngm) )
  !
  !     Self-consistency loop - set initial rhoin to homogeneous charge
  !     (not a smart choice, integrates to the correct total charge)
  !
  rhoin(:,:,:) = nelec/omega
  vg(:) = 0.0_dp
  !
  do iter=1, max_iter
     !
     !     Loop on k-vectors
     !
     rhoout(:,:,:) = 0.0_dp
     do n=1, nk
        !
        call fill_h ( npwx,  npw(n), igk(1,n), k(1,n), ngm, g, mill, &
                      nm1, nm2, nm3, indg, g2, vg, sg, tpi, a, h )
        !
        !       Solution [expansion coefficients are stored into psi(j,i)
        !                 j=basis function index, i= eigenvalue index]
        !
        call diag_h ( npwx, npw(n), h, nbands, e(1,n), psi )
        !
        call sum_charge ( npwx, psi, npw(n), igk(1,n), nbands, ngm, mill, &
                          nr1,nr2,nr3, wk(n), omega, rhoout )
        !
     end do
     !
     drho2 = mix_charge ( alpha, omega, nr1,nr2,nr3, rhoin, rhoout )
     !
     if ( drho2 < threshold ) then
        write(*,'("Convergence threshold (",ES10.3,") reached")') drho2
        do n=1, nk
           write (*,'("k = ",3f10.4,"  2pi/a")') k(:,n)
           write (*,'("e = ",4f10.4,"  eV")') e(1:nbands,n)*13.6058
        end do
        exit 
     else
        write(*,'("Iteration # ",i2,"  Delta rho:",ES10.3)') iter,drho2
     end if
     !
     ! new charge is now in rhoin - calculate new potential     
     !
     call v_of_rho ( omega, nelec, nr1, nr2, nr3, rhoin, ngm, mill, g2, vg )
     !
  end do
  !
  deallocate ( rhoin,  rhoout, vg )
  deallocate ( psi, h, e )
  deallocate ( igk )
  deallocate ( sg )
  deallocate ( indg )
  deallocate ( g, g2, mill ) 
  deallocate ( k, npw, wk )
  !
  stop
end program ah

!---------------------------------------------------------------
subroutine fill_h ( npwx, npw, igk, k, ngm, g, mill, nm1, nm2,nm3, indg, &
                    g2, vg, sg, tpi, a, h )
  !---------------------------------------------------------------
  !
  !        Fill the Hamiltonian matrix
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent(in) :: npwx, npw, ngm, nm1,nm2,nm3
  integer, intent(in) :: mill(3,ngm), igk(npwx)
  integer, intent(in) :: indg(-nm1:nm1,-nm2:nm2,-nm3:nm3)
  real (dp), intent(in) :: tpi, a, k(3), vg(ngm), sg(ngm), g(3,ngm), g2(ngm)
  real (dp), intent(out):: h(npwx, npwx)
  integer :: i, j, ik, jk, ng, n1,n2,n3
  real(dp) :: kg0(3), vsg
  real(dp), external :: form_factor
  !
  h(:,:) = 0.0_dp
  !
  !       Assign values of the matrix elements of the hamiltonian 
  !       on the plane wave basis
  !
  do i=1,npw
     ik = igk(i)
     kg0(:) = k(:) + g(:,ik)
     !  j >= i because only upper triangle needed
     do j=i,npw
        jk = igk(j)
        ! find G=G_i-G_j in the list of G-vectors via Miller indices
        n1 = mill(1,ik) - mill(1,jk)
        n2 = mill(2,ik) - mill(2,jk)
        n3 = mill(3,ik) - mill(3,jk)
        ng = indg(n1,n2,n3)
        if ( ng == 0 ) then
           print *, ' something seriously wrong here'
           stop
        end if
        vsg = form_factor ( g2(ng) )
        if ( i ==j ) then
           h(i,j) =  (tpi/a)**2 * (kg0(1)**2 + kg0(2)**2 + kg0(3)**2) + &
                      vsg * sg(ng) + vg(ng)
        else
           h(i,j) = vsg * sg(ng) + vg(ng) 
        end if
        !print  '(2i4,f12.6)', i,j, h(i,j)
     end do
  end do
end subroutine fill_h
!
!---------------------------------------------------------------
function form_factor (g2)
  !---------------------------------------------------------------
  !
  !  form factor for Appelbaum-Hamann pseudopotential
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi=3.14159265358979_dp, tpi=2.0_dp*pi, fpi=4.0_dp*pi
  ! a small quantity
  real(dp), parameter :: eps=1.0d-8
  ! conversion factor Hartree-Rydberg
  real(dp), parameter :: e2 = 2.0_dp
  ! parameters of the pseudopotential: PRB 8, 1777 (1973) - v1,v2 in Ha
  real(dp), parameter :: v1=3.042, v2=-1.372, alpha=0.6102, zv=4.0_dp
  real(dp), intent(in):: g2
  real(dp) :: form_factor

  if ( g2 < eps ) then
     !
     ! G=0: divergent Hartree and pseudopotential terms cancel
     !      what is left can be analytically calculated
     !
     form_factor = e2 * pi * zv / alpha &
            + e2 * (pi/alpha)**(3.0_dp/2.0_dp) * (v1 + 3.0_dp/2.0_dp*v2/alpha)
  else
     form_factor = e2 * exp(-g2/4.0_dp/alpha) * &
               ( - fpi*zv/g2 + (pi/alpha)**(3.0_dp/2.0_dp)* &
               ( v1 + v2/alpha * (3.0_dp/2.0_dp - g2/4.0_dp/alpha) ) ) 
  end if
end function form_factor

!---------------------------------------------------------------
subroutine sum_charge ( npwx, psi, npw, igk, nbands, ngm, mill, &
                        nr1,nr2,nr3, wk, omega, rho )
  !---------------------------------------------------------------
  !
  !        Calculate charge density in real space, sum over k-points
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent(in) :: npwx, npw, nbands, ngm, nr1, nr2, nr3
  integer, intent(in) :: mill(3,ngm), igk(npwx)
  real (dp), intent(in) :: psi(npwx,nbands), wk, omega
  real (dp), intent(inout) :: rho(nr1,nr2,nr3)
  !
  !    aux = complex array storing a valence band, or the charge density,
  !          or the potential, in real or G space, used to perform FFT's
  !
  complex (dp), allocatable :: aux(:,:,:)
  integer :: nb, n1,n2,n3, i, ik

  allocate (aux(nr1,nr2,nr3))
  do nb=1,nbands
     aux(:,:,:) = (0.0_dp, 0.0_dp)
     do i=1,npw
        ik = igk(i)
        n1 = mill(1,ik)
        if ( n1 < 0 ) then
           n1 = n1 + nr1 + 1
        else
           n1 = n1 + 1
        end if
        n2 = mill(2,ik)
        if ( n2 < 0 ) then
           n2 = n2 + nr2 + 1
        else
           n2 = n2 + 1
        end if
        n3 = mill(3,ik)
        if ( n3 < 0 ) then
           n3 = n3 + nr3 + 1
        else
           n3 = n3 + 1
        end if
        !
        ! Miller indices run from -n/2 to + n/2
        ! Indices n1,n2,n3 run from 1 to nr1,nr2,nr3
        ! with negative values refolded so that they
        ! fall in the "other side of the cell" in G space
        !
        aux(n1,n2,n3) = cmplx ( psi(i,nb), 0.0_dp, kind=dp )
     end do
     !
     call fft3d ( 1, nr1, nr2, nr3, aux )
     !
     ! Factor 2 is for spin degeneracy: two electrons per band
     ! Factor omega comes from the definition of plane waves
     !
     rho(:,:,:) = rho(:,:,:) + 2.0_dp*wk*abs(aux(:,:,:))**2/omega
     !
  end do
  !
  deallocate (aux)
  !
end subroutine sum_charge
!
!---------------------------------------------------------------
function mix_charge ( alpha, omega, nr1,nr2,nr3, rhoin, rhoout ) result (drho2)
  !---------------------------------------------------------------
  !
  ! simple mixing to get a new estimate of the charge density
  ! mixed charge density is stored in rhoin, rhoout unchanged
  ! drho2 = sqrt( \int |rhoout(r)-rhoin(r)|^2 dr )
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent(in) :: nr1, nr2, nr3
  real (dp), intent(in) :: alpha, omega, rhoout(nr1,nr2,nr3)
  real (dp), intent(inout) :: rhoin(nr1,nr2,nr3)
  !
  real (dp) :: drho2
  integer :: n1,n2,n3
  !
  drho2 = 0.0_dp
  do n3=1,nr3
     do n2=1,nr2
        do n1=1,nr1
           drho2 = drho2 + abs(rhoout(n1,n2,n3)-rhoin(n1,n2,n3))**2
        end do
     end do
  end do
  drho2 = sqrt(drho2 * omega/(nr1*nr2*nr3)) 
  !
  rhoin(:,:,:) = alpha*rhoin(:,:,:) + (1.0_dp - alpha)*rhoout(:,:,:)
  !
  return
end function mix_charge

!---------------------------------------------------------------
subroutine v_of_rho ( omega, nelec, nr1, nr2, nr3, rho, ngm, mill, g2, vg )
  !---------------------------------------------------------------
  !
  ! Calculate Hartree and exchange-correlation potential
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  ! conversion factor Hartree-Rydberg
  real(dp), parameter :: e2 = 2.0_dp
  real(dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi
  ! a small quantity
  real(dp), parameter :: eps=1.0d-8
  ! Slater alpha
  real (dp), parameter :: alphax = 2.0_dp/3.0_dp
  !
  integer, intent(in) :: nr1, nr2, nr3, ngm, nelec
  integer, intent(in) :: mill(3,ngm)
  real (dp), intent(in) :: omega
  real (dp), intent(in) :: rho(nr1,nr2,nr3), g2(ngm)
  real (dp), intent(out) :: vg(ngm)
  !
  !    rhog= charge density  in reciprocal space
  real (dp), allocatable :: rhog(:)
  !    aux = complex array used to perform FFT's
  complex (dp), allocatable :: aux(:,:,:)
  !    vr  = exchange and coulomb potential, in real space
  real (dp), allocatable :: vr(:,:,:)
  real (dp) :: charge
  integer :: n1,n2,n3,ng
  !
  ! Exchange (Slater) in real space
  !
  allocate ( vr(nr1,nr2,nr3) )
  do n3=1,nr3
     do n2=1,nr2
        do n1=1,nr1
           vr(n1,n2,n3) = - alphax * 3.0_dp / 2.0_dp * e2 * &
                ( 3.0_dp*rho(n1,n2,n3)/pi )**(1.0_dp/3.0_dp)
        end do
     end do
  end do
  !
  ! bring now v(r) to reciprocal space 
  !
  allocate (aux(nr1,nr2,nr3) )
  aux(:,:,:) = cmplx ( vr(:,:,:), 0.0_dp, kind=dp )
  call fft3d ( -1, nr1, nr2, nr3, aux )
  deallocate (vr)
  !
  vg(:) = 0.0_dp
  do ng=1,ngm
     n1 = mill(1,ng)
     if ( n1 < 0 ) then
        n1 = n1 + nr1 + 1
     else
        n1 = n1 + 1
     end if
     n2 = mill(2,ng)
     if ( n2 < 0 ) then
        n2 = n2 + nr2 + 1
     else
        n2 = n2 + 1
     end if
     n3 = mill(3,ng)
     if ( n3 < 0 ) then
        n3 = n3 + nr3 + 1
     else
        n3 = n3 + 1
     end if
     !
     vg(ng) = real( aux(n1,n2,n3), kind=dp )
     !
  end do
  !
  ! Now compute Hartree potential - we need rho(G) in reciprocal space
  !
  charge = sum( abs(rho(:,:,:)) ) * omega/(nr1*nr2*nr3)
  if ( abs(charge-nelec) > eps ) &
       write(*,'("Check: charge (real space)=",f12.8)') charge
  !
  aux(:,:,:) = cmplx ( rho(:,:,:), 0.0_dp, kind=dp )
  call fft3d ( -1, nr1, nr2, nr3, aux )
  ! now fill rho(G)
  allocate (rhog(ngm))
  rhog(:) = 0.0_dp
  do ng=1,ngm
     n1 = mill(1,ng)
     if ( n1 < 0 ) then
        n1 = n1 + nr1 + 1
     else
        n1 = n1 + 1
     end if
     n2 = mill(2,ng)
     if ( n2 < 0 ) then
        n2 = n2 + nr2 + 1
     else
        n2 = n2 + 1
     end if
     n3 = mill(3,ng)
     if ( n3 < 0 ) then
        n3 = n3 + nr3 + 1
     else
        n3 = n3 + 1
     end if
     !
     ! Miller indices run from -n/2 to + n/2
     ! Indices n1,n2,n3 run from 1 to nr1,nr2,nr3
     ! with negative values refolded so that they
     ! fall in "the far side of the cell" in G space
     !
     rhog(ng) = real( aux(n1,n2,n3), kind=dp )
     !
  end do

  ! Check: is rho(G=0)=Nelec/Omega ?
  charge = omega * real( aux(1,1,1), kind=dp )
  if ( abs(charge-nelec) > eps ) &
       write(*,'("Check: charge (G-space)   =",f12.8)') charge
  
  ! Hartree term in G space - beware the G=0 divergence!

  do ng =1, ngm
     if ( g2(ng) > eps ) vg(ng) = vg(ng) + fpi*e2*rhog(ng)/g2(ng)
  end do
  !
  deallocate (aux)
  deallocate (rhog)
  !
end subroutine v_of_rho
!
!---------------------------------------------------------------
subroutine diag_h ( npwx, npw, h, nbands, e, psi )
  !---------------------------------------------------------------
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: npwx, npw, nbands
  real(dp), intent(inout) :: h(npwx, npw)
  real(dp), intent(out) :: e(nbands), psi(npwx, nbands) 
  !
  integer :: lwork, info
  real(dp), allocatable :: work(:), ework(:)

  lwork = 3*npw
  allocate ( work(lwork), ework (npwx) )

  call dsyev ( 'V', 'U', npwx, h, npw, ework, work, lwork, info )
  if (info /= 0) stop 'H-matrix diagonalization failed '

  ! copy desired results to output variables

  e(:) = ework (1:nbands)
  psi(:,1:nbands) = h(:,1:nbands) 

  deallocate (ework, work)

end subroutine diag_h
!
subroutine fft3d ( dir, n1, n2, n3, cin )
  !
  ! interface for 3D FFTW v.3
  ! dir is the sign of the exponential:
  !   dir  = -1  forward  (R->G) fft,
  !   dir  = +1  backward (G->R) fft
  ! n1, n2, n3   FFT dimension
  ! cin (1:n1*n2*n3) = complex array containing input data 
  !                    On output, cin contains transformed data
  !
  ! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
  use, intrinsic :: iso_c_binding
  implicit none
  ! fftw3.f03 provides definitions and variables used by FFTW v.3
  include 'fftw3.f03'
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent(in) :: dir
  integer, intent(in) :: n1, n2, n3
  complex(dp), intent(inout) :: cin (n1,n2,n3)
  !
  !   Keep track of fft dimension for plan initialization
  integer, save :: n1s = -1, n2s = -1, n3s = -1
  !   Pointers to the C structures (plan) containing FFT factors
  type(C_PTR), save :: fw_plan = C_NULL_PTR
  type(C_PTR), save :: bw_plan = C_NULL_PTR
  !
  !
  if( n1 < 1 .or. n2 < 1 .or. n3 < 1 ) then
     print *, "Error: fft dimensions out of range ", n1,n2,n3
     stop
  end if
  if( dir /= -1 .and. dir /= 1 ) then
     print *, "Error: incorrect value for dir ", dir
     stop
  end if
  !
  if( n1 /= n1s .or. n2 /= n2s .or. n3 /= n3s ) then
     !
     !   initialize plans if not previously initialized,
     !   or if a new initialization is needed
     !   Note the transposition of the three dimensions!
     !   It is needed for C-Fortran interoperability
     !
     IF( C_ASSOCIATED(fw_plan) ) CALL fftw_destroy_plan( fw_plan )
     IF( C_ASSOCIATED(bw_plan) ) CALL fftw_destroy_plan( bw_plan )
     ! The reversal of n1, n2, n3 takes care of the different 
     ! fortran-C ordering of arrays (fftw assumes the C logic)
     fw_plan = fftw_plan_dft_3d ( n3, n2, n1, cin, cin, FFTW_FORWARD, &
          FFTW_ESTIMATE)
     bw_plan = fftw_plan_dft_3d ( n3, n2, n1, cin, cin,FFTW_BACKWARD, &
          FFTW_ESTIMATE)
     !
     n1s = n1; n2s = n2; n3s=n3
     !
  END IF
  !
  IF ( dir < 0 ) THEN
     CALL fftw_execute_dft( fw_plan, cin, cin )
     cin = cin / dble(n1*n2*n3)
  ELSE IF ( dir > 0 ) THEN
     CALL fftw_execute_dft( bw_plan, cin, cin )
  END IF
  !
END SUBROUTINE fft3d
