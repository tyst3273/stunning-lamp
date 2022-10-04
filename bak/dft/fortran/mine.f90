! --------------------------------------------------------------------------------------------------
program mine
! --------------------------------------------------------------------------------------------------
    ! 
    !   SCF program for Si based on Paolo Gianozzi's example from one of his courses
    !   listed on his website. 
    !
    ! ---------------------------------------------------------------------------------------------- 
    !
    implicit none
    !
    integer,  parameter :: dp  = selected_real_kind(14,200)
    real(dp), parameter :: pi  = 3.14159265_dp, twopi = 2.0_dp*pi
    real(dp), parameter :: eps = 1.0d-8
    ! 
    !   a       = the lattice parameter in Bohr
    !   alpha   = scf mixing parameter
    !   v_tol   = drho2 threshold for scf convergence
    !
    real(dp) :: a = 10.26, alpha = 0.5, v_tol = 1.0d-6
    !
    !   atom positions are +- tau
    !
    real(dp) :: tau1 = 0.125, tau2 = 0.125, tau3 = 0.125
    !
    !   nelec       = number of electrons in unitcell
    !   nk          = number of k points
    !   nbnd        = number of occupied bands
    !   max_iter    = maximum number of allowed scf iterations
    !   
    integer :: nelec = 8, nk = 1, nbnd = 4, max_iter = 50
    !
    !   npwx    = max number of plane waves across k-points
    !   ngm     = number of G vectors with G^2 < 4*Ecut
    !   nm*     = Max value of Miller indices
    !   nr*     = real-space grid dimensions
    !
    integer :: npwx, ngm, nm1, nm2, nm3, nr1, nr2, nr3
    !
    !   npw     = number of planewaves per kpoint
    !   mill    = miller indicies of G vectors
    !   indg    = ind of G vectors for given miller index
    !   igk     = ind of G vector in list of k+G such that (k+G)**2 < Ecut
    !
    integer, allocatable :: npw(:), mill(:,:), indg(:,:,:), igk(:,:)
    !
    !   ecut        = pw energy cutoff in ry
    !   a1, a2, a3  = direct lattice vectors
    !   b1, b2, b3  = reciprocal lattice vectors
    !   omega       = volume of unit cell
    !
    real(dp) :: ecut, a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), omega
    !
    !   k   = k-points, in units of 2*pi/a
    !   wk  = normalized kpoint weights
    !   e   = energy eigenvals
    !   h   = eigenvectors (overwrites hamiltonian)
    !   g   = G vectors in units of 2*pi/a
    !   g2  = G^2 in atomic units
    !   sg  = geometrical structure factor
    !
    real(dp), allocatable :: k(:,:), wk(:), e(:,:), h(:,:), g(:,:), g2(:), sg(:)
    !
    !   rhoi    = input charge density in a given scf step, real space
    !   rhoo    = output charge density in a given scf step, real space
    !   vg      = exhange/hartree potential in G space
    !
    real(dp), allocatable :: rhoi(:,:,:), rhoo(:,:,:), vg(:)
    !
    !   psi     = valence bands (in planwave basis). taken to be real due to inversion symmetry
    !   
    real(dp), allocatable :: psi(:,:)
    !
    !   work variables
    !
    real(dp) :: vsg, gg, g0(3), kg2, kg0(3), r(3), drho2
    real(dp), external  :: form_factor, mix_charge
    integer :: ii, jj, ik, jk, n, ng, nn, n1, n2, n3, iter
    !
    !   set up the plane wave basis set
    !   G_n = n*2*pi/L. expand up to hbar**2*G**2/2/m <= Ecut
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !   get input args and set up G vectors, kpoints, etc
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !   get ecut from user (or use hardcoded version below)
    !
!    write(*,"('Cutoff for plane waves: ecut (Ry) > ')",advance='no')
!    read(*,*) ecut
!    if( ecut <= 0.0_dp ) stop ' wrong cutoff'
    !
    !   hardcoded !!!
    !
    ecut = 1
    !
    write(*,*)
    write(*,'("energy cutoff (Ry): ",f6.2)') ecut
    write(*,*)
    !
    !   allocate (and set) the kpoint arrays
    !
    allocate( k(3,nk), wk(nk), npw(nk) )
    k(1,1) = 0.6223_dp
    k(2,1) = 0.2953_dp
    k(3,1) = 0.0_dp
    wk(1)  = 1.0_dp
    !
    !   set up the fcc lattice
    !
    a1(1) =  0.5_dp; a1(2) =  0.5_dp; a1(3) =  0.0_dp
    a2(1) =  0.5_dp; a2(2) =  0.0_dp; a2(3) =  0.5_dp
    a3(1) =  0.0_dp; a3(2) =  0.5_dp; a3(3) =  0.5_dp
    omega = a**3/4.0_dp
    !
    !   get fcc reciprocal lattice
    !
    b1(1) =  1.0_dp; b1(2) =  1.0_dp; b1(3) = -1.0_dp
    b2(1) =  1.0_dp; b2(2) = -1.0_dp; b2(3) =  1.0_dp
    b3(1) = -1.0_dp; b3(2) =  1.0_dp; b3(3) =  1.0_dp
    !
    !   count the G vectors such that hbar**2/2/m*G**2 < 4*Ecut
    !   this calculates |b|=sqrt(b*b) and finds n_max from |G_max| = n_max*|b|
    !
    nm1 = nint(sqrt(4.0_dp*ecut)/twopi*a*sqrt(a1(1)**2+a1(2)**2+a1(3)**2)+0.5)
    nm2 = nint(sqrt(4.0_dp*ecut)/twopi*a*sqrt(a2(1)**2+a2(2)**2+a2(3)**2)+0.5)
    nm3 = nint(sqrt(4.0_dp*ecut)/twopi*a*sqrt(a3(1)**2+a3(2)**2+a3(3)**2)+0.5)
    !
    write(*,*)
    write(*,'("estimates for n_max: ",3I6)') nm1, nm2, nm3
    write(*,*)
    !
    !   find largest G vector
    !
    ngm = 0
    do n1 =  -nm1, nm1
        do n2 = -nm2, nm2 
            do n3 = -nm3, nm3
                g0(:) = n1*b1(:)+n2*b2(:)+n3*b3(:)
                gg = (twopi/a)**2*(g0(1)**2+g0(2)**2+g0(3)**2)
                if( gg <= 4.0_dp*ecut ) ngm = ngm+1
            end do
        end do
    end do
    write(*,*)
    write(*,'("number of G-vectors: ",i6)') ngm
    write(*,*)
    !
    !   now that we know the size, allocate the arrays and then fill them
    !
    allocate( g(3,ngm), g2(ngm), mill(3,ngm) ) 
    !
    ng = 0
    do n1 = -nm1, nm1
        do n2 = -nm2, nm2
            do n3 = -nm3, nm3
                g0(:) = n1*b1(:)+n2*b2(:)+n3*b3(:)
                gg = (twopi/a)**2*(g0(1)**2+g0(2)**2+g0(3)**2)
                if( gg <= 4.0_dp*ecut ) then
                    ng = ng+1
                    g(:,ng) = g0(:)
                    g2(ng) = gg
                    mill(1,ng) = n1 ; mill(2,ng) = n2 ; mill(3,ng) = n3
                end if
            end do
        end do
    end do
    if( .not. ng == ngm ) then
        write(*,*)
        write(*,*) '** Error **'
        write(*,'("Mismatch in number of G-vectors! ",2i6)') ng, ngm
        write(*,*) 
        stop
    end if
    !
    !   indg givs index of a G from its miller indiceis
    !
    nm1 = maxval(mill(1,:))
    nm2 = maxval(mill(2,:))
    nm3 = maxval(mill(3,:))
    allocate( indg(-nm1:nm1, -nm2:nm2, -nm3:nm3) )
    indg (:,:,:) = 0
        do ng = 1, ngm
            indg ( mill(1, ng), mill(2,ng), mill(3,ng) ) = ng
    end do
    !
    !   Geometrical structure factors S(G) for the atoms.
    !   the pseudopotential is set at the origin, so in G space, it comes with a phase
    !   S(G) = cos(G*T)
    !
    allocate( sg(ngm) )
    sg(:) = 0.0_dp
    do ng = 1, ngm
        sg(ng) = 2.0_dp*cos(twopi*((g(1,ng)*tau1+g(2,ng)*tau2+g(3,ng)*tau3)))
    end do
    sg(:) = sg(:)/omega
    !
    !   get dimensions of real space grid to store real space quantities   
    !
    nr1 = 2*nm1+2
    nr2 = 2*nm2+2
    nr3 = 2*nm3+2
    write(*,*)
    write(*,'("Real space grid: ",3i6)') nr1, nr2, nr3
    write(*,*)
    !
    !   count the plane waves with \hbar**2/2/m*(k+G)**2 < Ecut
    !
    write(*,*)
    do n = 1, nk
        !
        npw(n) = 0
        do ng = 1, ngm
            ! 
            !   get k+G in units of 2*pi/a
            !
            kg0(:) = k(:,n)+g(:,ng)
            kg2 = (twopi/a)**2*(kg0(1)**2+kg0(2)**2+kg0(3)**2)
            if( kg2 <= ecut ) npw(n) = npw(n)+1
        end do
        write(*,'("num pw for k vector: ",i3," = ",i6)') n, npw(n)
    end do
    write(*,*)
    !
    !   get max num pw over all k points
    !
    npwx = maxval(npw(:))
    !
    !   allocate the k+G indicies arrays
    !
    allocate( igk(npwx,nk) )
    do n = 1, nk
        nn = 0
        do ng = 1, ngm
            kg0(:) = k(:,n)+g(:,ng)
            kg2 = (twopi/a)**2*(kg0(1)**2+kg0(2)**2+kg0(3)**2)
            if( kg2 <= ecut ) then
                nn = nn+1
                igk(nn,n) = ng
            end if
        end do
        if( .not. nn == npw(n) ) then
            write(*,*)
            write(*,*) '** Error **'
            write(*,'("Mismatch in number of plane waves! ",3i6)') n, nn, npw(n)
            write(*,*)
            stop
        end if
    end do
    !
    !   allocate the basis set and hamiltonian, 
    !   realspace density grids
    !
    allocate( e(nbnd,nk), h(npwx,npwx), psi(npwx,nbnd) )
    allocate( rhoi(nr1,nr2,nr3), rhoo(nr1,nr2,nr2), vg(ngm) )
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !   self-consistent-loop  -- set initial rhoi to homogeneous charge normalized to
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !   equal the total charge. ( its a bad choice, should be atomic orbitals )
    !   also, since grad rho == 0, vg = 0
    !
    rhoi(:,:,:) = nelec/omega
    vg(:) = 0.0_dp
    !
    do iter = 1, max_iter
        !
        !   loop over the k-vectors
        !   zero rhoo before filling the required matrix elements
        !
        rhoo(:,:,:) = 0.0_dp
        do n = 1, nk
            !
            !   fill the hamiltonian for this k-vector
            !
            call fill_h( npwx, npw(n), igk(1,n), k(1,n), ngm, g, mill, & 
                         nm1, nm2, nm3, indg, g2, vg, sg, twopi, a, h )
            !
            !   diagonalize the hamiltonian. expansion coeffs. go into psi(ii,jj)
            !   ii = basis funcion index, jj = eigenvalue index
            !
            call diag_h( npwx, npw(n), h, nbnd, e(1,n), psi )
            !
            !   sum k-points over 
            !
            call sum_charge( npwx, psi, npw(n), igk(1,n), nbnd, ngm, mill, &
                             nr1, nr2, nr3, wk(n), omega, rhoo )
            !
        end do
        write(*,*) e
        !
        !   now mix the charge for the next scf step
        !
        drho2 = mix_charge( alpha, omega, nr1, nr2, nr3, rhoi, rhoo )
        !
        !   check if convergence threshold is achieved
        !
        if( drho2 < v_tol ) then
            !
            !   if so, print the results and exit!
            !
            write(*,*)
            write(*,'("Convergence threshold: ",es10.3," reached ")') drho2
            write(*,*)
            do n = 1, nk
                write (*,'("k = ",3f10.4,"  2pi/a")') k(:,n)
                write (*,'("e = ",4f10.4,"  eV")')    e(1:nbnd,n)*13.6058 ! Ry to eV
            end do
            write(*,*)
            write(*,*) 'Goodbye!'
            write(*,*)
            exit
        else
            write(*,'("Iter # : ",i3," drho2: ",es10.3)') iter, drho2
        end if
        !
        !   new charge was put in rhoi. now calculate new potential potential
        !
        call v_from_rho( omega, nelec, nr1, nr2, nr3, rhoi, ngm, mill, g2, vg )
        !
    end do
    !
    deallocate( rhoi, rhoo, vg )
    deallocate( psi, h, e )
    deallocate( sg )
    deallocate( indg )
    deallocate( g, g2, mill )
    deallocate( k, npw, wk )
    !
    stop
    !
end program mine
!
! --------------------------------------------------------------------------------------------------
!
!   subroutines that are used above
!
! --------------------------------------------------------------------------------------------------
subroutine fill_h( npwx, npw, igk, k, ngm, g, mill, nm1, nm2, nm3, indg, g2, vg, sg, twopi, a, h )
! --------------------------------------------------------------------------------------------------
    !
    !   subroutine that is used to calculate matrix elements of the hamiltonian
    !
    ! ----------------------------------------------------------------------------------------------
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    !
    !   input variables
    !
    integer, intent(in) :: npwx, npw, ngm, nm1, nm2, nm3
    integer, intent(in) :: mill(3,ngm), igk(npwx)
    integer, intent(in) :: indg( -nm1:nm1, -nm2:nm2, -nm3:nm3 )
    real(dp), intent(in) :: twopi, a, k(3), vg(ngm), sg(ngm), g(3,ngm), g2(ngm)
    !
    !   output variables
    !
    real(dp), intent(out) :: h(npwx, npwx)
    !
    !   local work variables
    !
    integer :: ii, jj, ik, jk, ng, n1, n2, n3
    real(dp) :: kg0(3), vsg 
    !
    !   function handles
    !
    real(dp), external :: form_factor
    !
    !   zero the hamiltonian
    !
    h(:,:) = 0.0_dp
    !
    !   calculate the matrix elements
    !
    do ii = 1, npw
        ik = igk(ii) ! get the index of this k point in the gk array
        kg0(:) = k(:)+g(:,ik)
        !
        !   only do upper triangle, save time
        !
        do jj = ii, npw
            jk = igk(jj)
            !
            !   get G=G-i-G_j in the list of G vectors. (units of G are miller indicies)
            !
            n1 = mill(1,ik)-mill(1,jk)
            n2 = mill(2,ik)-mill(2,jk)
            n3 = mill(3,ik)-mill(3,jk)
            ng = indg(n1,n2,n3)
            !
            !   check that something bad didnt happen
            !
            if( ng == 0 ) then
                write(*,*)
                write(*,*) ' ** ERROR ** '
                write(*,*) 'something is seriously wrong ... '
                write(*,*)
                stop
            end if
            !
            !   get the fourier transform of the pseudopotential (i.e its matrix element)
            !   note, its calculated with the atom at the origin, so we need the structure 
            !   factor sg too
            !
            vsg = form_factor( g2(ng) )
            !
            !   for diagonal elements, we need kinetic energy term
            !
            if( ii == jj ) then
                h(ii,jj) = (twopi/a)**2*(kg0(1)**2+kg0(2)**2+kg0(3)**2)+vsg*sg(ng)+vg(ng)
            !
            ! otherwise, non diagonal elements
            !
            else
                h(ii,jj) = vsg*sg(ng)+vg(ng)
            end if
            !
        end do
    end do
    !
end subroutine fill_h
!
! --------------------------------------------------------------------------------------------------
function form_factor(g2)
! --------------------------------------------------------------------------------------------------
    !
    !   form factor for the Appelbaum Hamann pseudopotential
    !
    ! ----------------------------------------------------------------------------------------------
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    real(dp), parameter :: pi  = 3.14159265_dp, twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
    real(dp), parameter :: eps = 1.0d-8
    real(dp), parameter :: h2ry = 2.0_dp
    !
    !   parameters from the pseudopotential
    !
    real(dp), parameter :: v1 = 3.042, v2 = -1.372, alpha = 0.6102, zv = 4.0_dp
    !
    !   input variables
    !
    real(dp), intent(in) :: g2
    !
    !   output variables
    !
    real(dp) :: form_factor
    !
    !
    !
    if( g2 < eps ) then
        !
        !   G = 0 is divergent. Hartree and pseudocancel eachother and the result at G==0
        !   can be solved analytically
        !
        form_factor = h2ry*(pi*zv/alpha+(pi/alpha)**(3.0_dp/2.0_dp)*(v1+3.0_dp/2.0_dp*v2/alpha))
        !
    else
        form_factor = h2ry*exp(-g2/4.0_dp/alpha)* &
               (-fourpi*zv/g2+(pi/alpha)**(3.0_dp/2.0_dp)*(v1+v2/alpha*(3.0_dp/2.0_dp-g2/4.0_dp/alpha)))
    end if
    !
end function form_factor
!
! --------------------------------------------------------------------------------------------------
subroutine diag_h( npwx, npw, h, nbnd, e, psi )
! --------------------------------------------------------------------------------------------------
    !
    !   diagonalize the hamiltonian using lapack routines
    !
    ! ----------------------------------------------------------------------------------------------
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    !
    !   input vars
    !
    integer, intent (in) :: npwx, npw, nbnd
    !
    !   modified in place vars
    !
    real(dp), intent(inout) :: h(npwx, npw)
    !
    !   output vars
    !
    real(dp), intent(out) :: e(nbnd), psi(npwx, nbnd)
    !
    !   local variables
    !
    integer :: lwork, info
    real(dp), allocatable :: work(:), ework(:)
    !
    !   i dont know why we need these
    !
    lwork = 3*npw
    allocate( work(lwork), ework(npwx) )
    !
    !   diagonalize the hamiltonian using lapack call
    !
    call dsyev( 'V', 'U', npwx, h, npw, ework, work, lwork, info )
    if (info /= 0) then
        write(*,*)
        write(*,*) ' ** ERROR **'
        write(*,*) ' H-matrix diagonalization failed! fuck! '
        write(*,*)
        stop
    end if
    !
    !   copy desired results to output variable
    !
    e(:) = ework(1:nbnd)
    psi(:,1:nbnd) = h(:,1:nbnd)
    !
    !   free up the memory
    !
    deallocate( ework, work)
    !
end subroutine diag_h
!
! --------------------------------------------------------------------------------------------------
subroutine sum_charge( npwx, psi, npw, igk, nbnd, ngm, mill, nr1, nr2, nr3, wk, omega, rho )
! --------------------------------------------------------------------------------------------------
    !
    !   calculate charge density in real space by summing over k points
    !
    ! ----------------------------------------------------------------------------------------------
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    !
    !   input vars
    !
    integer, intent(in) :: npwx, npw, nbnd, ngm, nr1, nr2, nr3
    integer, intent(in) :: mill(3,ngm), igk(npwx)
    real(dp), intent(in) :: psi(npwx,nbnd), wk, omega
    !
    !   in place args
    !
    real(dp), intent(inout) :: rho(nr1,nr2,nr3)
    !
    !   local/work args
    !
    !   aux = complex array storing a valence band, or the charge density,
    !          or the potential, in real or G space, used to perform FFT's
    !
    complex (dp), allocatable :: aux(:,:,:)
    integer :: nb, n1, n2, n3, ii, ik
    !
    allocate( aux( nr1, nr2, nr3 ) )
    !
    !   loop over occupied bands and for each band, over the orbitals, to get thier coeffs.
    !   the coeffs are put into a G space grid, FFT'd to R space, then the density is computed
    !   for those orbtials and added to the total density
    !
    do nb = 1, nbnd
        aux(:,:,:) = (0.0_dp, 0.0_dp)
        !
        !   loop over the basis functions in this orbital
        !
        do ii = 1, npw
            ! 
            !   which k point we are doing
            !
            ik = igk(ii)
            !
            !   assemble to wfns into a 'grid' in G space
            !   
            n1 = mill(1,ik)
            if( n1 < 0 ) then
                n1 = n1+nr1+1
            else 
                n1 = n1+1
            end if
            !
            n2 = mill(2,ik)
            if( n2 < 0 ) then
                n2 = n2+nr2+1
            else
                n2 = n2+1
            end if
            !
            n3 = mill(3,ik)
            if( n3 < 0 ) then
                n3 = n3+nr3+1
            else
                n3 = n3+1
            end if
            !
            !   put the wavefunctions into 'grid' in G space
            !
            aux(n1,n2,n3) = cmplx( psi(ii,nb), 0.0_dp, kind=dp )
            !
        end do
        !
        !   now get them on 'grid' in R space
        !
        call fft3d( 1, nr1, nr2, nr3, aux )
        !
        !   now calculate rho as sum |psi|**2
        !
        rho(:,:,:) = rho(:,:,:)+2.0_dp*wk*abs(aux(:,:,:))**2/omega
        !
    end do
    !
end subroutine sum_charge
!
! --------------------------------------------------------------------------------------------------
function mix_charge( alpha, omega, nr1, nr2, nr3, rhoi, rhoo ) result( drho2 )
! --------------------------------------------------------------------------------------------------
    !
    !   simple mixing of the input/output charges to get new estimate of input
    !   charge for next iter step. mixhed charge goes into rhoi. 
    !   drho2 = sqrt( \integral |rhoo(r)-rhoi(r)|**2 dr )
    !   
    ! ----------------------------------------------------------------------------------------------
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    !
    !   input vars
    !
    integer, intent(in) :: nr1, nr2, nr3
    real (dp), intent(in) :: alpha, omega, rhoo(nr1,nr2,nr3)
    !
    !   in-place vars
    !
    real (dp), intent(inout) :: rhoi(nr1,nr2,nr3)
    !
    !   local work variables
    !
    real (dp) :: drho2
    integer :: n1,n2,n3
    ! 
    drho2 = 0.0_dp
    do n3 = 1, nr3
        do n2 = 1, nr2
            do n1 = 1, nr1
                drho2 = drho2+abs(rhoo(n1,n2,n3)-rhoi(n1,n2,n3))**2
            end do
        end do
    end do
    drho2 = sqrt(drho2*omega/(nr1*nr2*nr3)) ! dr = omega/(nr1*nr2*nr3)
    !
    !   mix the new and old charge density
    !
    rhoi(:,:,:) = alpha*rhoi(:,:,:)+(1.0_dp-alpha)*rhoo(:,:,:)
    !
    return
    !
end function mix_charge
!
! --------------------------------------------------------------------------------------------------
subroutine v_from_rho( omega, nelec, nr1, nr2, nr3, rho, ngm, mill, g2, vg )
! --------------------------------------------------------------------------------------------------
    !
    !   calculate hartree and xc potentials from the density
    !
    ! ----------------------------------------------------------------------------------------------
    !
    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    real(dp), parameter :: pi  = 3.14159265_dp, twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
    real(dp), parameter :: eps = 1.0d-8
    real(dp), parameter :: h2ry = 2.0_dp
    real (dp), parameter :: alphax = 2.0_dp/3.0_dp 
    !
    !   input params
    !
    integer, intent(in) :: nr1, nr2, nr3, ngm, nelec
    integer, intent(in) :: mill(3,ngm)
    real (dp), intent(in) :: omega
    real (dp), intent(in) :: rho(nr1,nr2,nr3), g2(ngm)
    !
    !   output params
    !
    real (dp), intent(out) :: vg(ngm)
    !
    !   rhog = charge density in reciprocal space
    !
    real (dp), allocatable :: rhog(:)
    !
    !   local/work args
    !
    !   aux = complex array storing a valence band, or the charge density,
    !          or the potential, in real or G space, used to perform FFT's
    !   vr  = exchange and coulomb potential, in real space
    !
    complex (dp), allocatable :: aux(:,:,:)
    real (dp), allocatable :: vr(:,:,:)
    real (dp) :: charge
    integer :: n1,n2,n3,ng
    !
    !   Slater exhange calculated in real space
    !
    allocate( vr( nr1, nr2, nr3 ) )
    do n3 = 1, nr3
        do n2 = 1, nr2
            do n1 = 1, nr1
                vr(n1,n2,n3) = -alphax*3.0_dp/2.0_dp*h2ry*(3.0_dp*rho(n1,n2,n3))**(1.0_dp/3.0_dp)
            end do
        end do
    end do
    !
    !   now get v(r) in G space
    !
    allocate( aux( nr1, nr2, nr3 ) )
    aux(:,:,:) = cmplx ( vr(:,:,:), 0.0_dp, kind=dp )
    call fft3d ( -1, nr1, nr2, nr3, aux )
    deallocate( vr )
    !
    vg(:) = 0.0_dp
    do ng = 1, ngm
        !
        n1 = mill(1, ng)
        if( n1 < 0 ) then
            n1 = n1+nr1+1
        else
            n1 = n1 + 1
        end if
        !
        n2 = mill(2,ng)
        if ( n2 < 0 ) then
            n2 = n2 + nr2 + 1
        else
            n2 = n2 + 1
        end if
        !
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
    !   now we need to calculate hartree potential -- use rho(G) in reciprocal space
    !
    charge = sum( abs( rho(:,:,:) )) * omega/(nr1*nr2*nr3)
    if ( abs(charge-nelec) > eps ) then
        write(*,*)
        write(*,*) ' ** WARNING ** '
        write(*,'(" check charge (real space): ",f12.8)') charge
        write(*,*)
    end if
    aux(:,:,:) = cmplx ( rho(:,:,:), 0.0_dp, kind=dp )
    call fft3d ( -1, nr1, nr2, nr3, aux )
    !
    !   now fill rho(G)
    !
    allocate( rhog(ngm) )
    rhog(:) = 0.0_dp
    do ng=1,ngm
        !
        n1 = mill(1,ng)
        if( n1 < 0 ) then
            n1 = n1 + nr1 + 1
        else
            n1 = n1 + 1
        end if
        !
        n2 = mill(2,ng)
        if( n2 < 0 ) then
            n2 = n2 + nr2 + 1
        else
            n2 = n2 + 1
        end if
        !
        n3 = mill(3,ng)
        if( n3 < 0 ) then
            n3 = n3 + nr3 + 1
        else
            n3 = n3 + 1
        end if
        !
    end do
    !
    !   check that rho(G=0) is correct
    !
    charge = omega*real(aux(1,1,1), kind=dp)
    if ( abs(charge-nelec) > eps ) then
        write(*,*) ' ** WARNING ** '
        write(*,'(" check charge (recip. space): ",f12.8)') charge
        write(*,*)
    end if
    !
    !   compute hartree potential in G space (solve poisson's equation)
    !
    do ng = 1, ngm
        if( g2(ng) > eps ) vg(ng) = vg(ng) + fourpi*h2ry*rhog(ng)/g2(ng)
    end do
    !
end subroutine v_from_rho
!
! --------------------------------------------------------------------------------------------------
subroutine fft3d ( dir, n1, n2, n3, cin )
! --------------------------------------------------------------------------------------------------
    !
    !   interface for 3D FFTW v.3
    !   dir is the sign of the exponential:
    !       dir  = -1  forward  (R->G) fft,
    !       dir  = +1  backward (G->R) fft
    !   n1, n2, n3   FFT dimension
    !   cin (1:n1*n2*n3) = complex array containing input data 
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
        !
        !   The reversal of n1, n2, n3 takes care of the different
        !   fortran-C ordering of arrays (fftw assumes the C logic)
        !
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
end subroutine fft3d
!
!
