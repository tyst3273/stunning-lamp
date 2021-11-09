!---------------------------------------------------------------
program heisenberg_exact
!---------------------------------------------------------------
  !
  !  Exact solution of isotropic 1-d Heisenberg model for s=1/2
  !  Periodic Boundary Conditions. energy in units of J.
  !  Lanczos algorithm vs conventional Hamiltonian diagonalization
  !  (Hamiltonian stored as a matrix - limited to N<15 spins or so)
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14, 200)
  integer, parameter :: Jsign =-1 ! sign of J: -1 AFM, +1 FM
  integer :: iseed = 1234567 ! not sure whether this is a good value
  integer :: N, nup, nhil, nl, ncount, nonzero
  integer :: ii, jj, lwork, info, i, j, k, ibip
  integer, allocatable :: states(:), conf(:), neigh(:)
  real(dp), allocatable :: H(:,:), W(:), work(:), fat(:)
  real(dp), allocatable :: e(:), d(:), eaux(:), daux(:), v0(:), v1(:), v2(:)
  real(dp), external :: dnrm2, ddot
  integer, external :: search, bin
  !
  !
10 continue
  write(6,'("Number of sites, number of up spins > ")', advance='no')
  read(5,*) N, nup
  if ( N <= 0 ) stop
  if ( N > 31 ) then 
     write(6,'("N too big, not enough bits in an integer")')
     stop
  end if
  !
  ! nhil = N!/nup!/ndw! is the dimension of the Hilbert space
  ! Directly compute binomial coefficient by recursion
  !
  nhil = bin(N,nup)
  !
  ! deallocate ( fat )
  !
  allocate ( conf(N), neigh(N) )
  allocate ( H(nhil,nhil), states(nhil) )

  do j=1,N
     i=0
     ! periodic boundary conditions
     if ( j == N ) i=N
     ! neigh(j) is the index of the neighbor of the j-th spin
     !          (only the one at the right to prevent double counting)
     neigh(j)=j+1-i
   enddo

   ncount=0
   ! the do loop runs on all possible 2^N states
   ! the bit representation is used to represent up and down spins
   do k=0, 2**N-1
      ! count all up spin 
      ibip=0
      do j=0, N-1
         ! btest(k,j) = .true. if bit j in k is 1 (i.e. spin j is up)
         ! following line is bogus - intel compiler bug
         if ( j  == N ) print *, k,j,btest(k,j)
         if ( btest(k,j) ) ibip=ibip+1
      enddo
      if ( ibip == nup ) then
         ! k labels a basis set state with Nup spin up
         ncount=ncount+1
         states(ncount)=k
         ! write(6,*) ncount,states(ncount)
      endif
      if ( ncount > nhil ) then
         write(6,'("Dimension of Hilbert space:",i6," is insufficient")') nhil
         stop
      endif
   enddo
   if ( ncount /= nhil ) then
      write(6,'("Dimension of Hilbert space mismatch: ",2i6)') nhil, ncount
      stop
   else
      write(*,'("Dimension of Hilbert space:",i6)') nhil 
   endif

   H=0.0_dp
   
   ! Fill the Hamiltonian: (J/2)\sum_{ij} S_-(i)S_+(j) + h.c. term
   nonzero = 0
   do ii=1,nhil
      ! do on all states (ii) of the Hilbert space
      do j=0,N-1
         ! conf(j+1) = -1 if j-th bit of state ii is 0 (i.e. spin down)
         ! conf(j+1) = +1 if j-th bit of state ii is 1 (i.e. spin up)
         conf(j+1)=-1
         if ( btest(states(ii),j) ) conf(j+1) = 1
      enddo
      ! conf(j) = sequence of + and - spins for state ii
      do j=1,N
         i=neigh(j)
         ! i is the index of the neighbor of spin j-th
         if ( conf(j) == -1 .and. conf(i) == 1 ) then
            ! state |...sigma_i...sigma_j...> with sigma_j=-1, sigma_i=+1
            ! nonzero contribution from term (1/2)\sum_{ij} S_-(i)S_+(j)
            k = states(ii)+2**(j-1)-2**(i-1)
            ! k is the label of state <...sigma_i...sigma_j...|
            !                    with sigma_j=+1, sigma_i=-1
            jj = search ( states, nhil, k )
            ! jj is the index of state k: states(jj) = k
            H(ii,jj)=H(ii,jj) - 0.5_dp*Jsign
            nonzero = nonzero + 1
         else if ( conf(j) == 1 .and. conf(i) == -1 ) then
            ! as above for (1/2)\sum_{ij} S_-(i)S_+(j)
            k = states(ii)-2**(j-1)+2**(i-1)
            jj = search ( states, nhil, k )
            H(ii,jj)=H(ii,jj) - 0.5_dp*Jsign
            nonzero = nonzero + 1
         endif
      enddo
   enddo
   ! fill the Hamiltonian: J\sum_{ij} S_z(i)S_z(j) term
   ! (contributes only to the diagonal of H)
   do ii=1,nhil
      do j=0,N-1
         conf(j+1)=-1
         if ( btest(states(ii),j) ) conf(j+1)=1
      enddo
      do j=1,N
         i=neigh(j)
         H(ii,ii)=H(ii,ii) - Jsign*0.25_dp*conf(j)*conf(i)
      enddo
      nonzero = nonzero + 1
   enddo
   !
   write(6,'("Number of nonzero elements:",i6," (",f5.2,"% of total)")') &
        nonzero, float(nonzero)/nhil**2*100
   !
   write(6,'("Number of Lanczos steps > ")', advance='no')
   read (5,*) nl
   nl = min ( nl, nhil )
   ! d, e vectors containing the Hamiltonian in tridiagonal form
   allocate ( d(nl), e(0:nl) )
   ! vectors used in the Lanczos algorithm
   allocate ( v0(nhil), v1(nhil), v2(nhil) )
   ! fill starting vector with random numbers
   call random_seed(size=iseed)
   call random_number ( harvest=v1 )
   e(0) = dnrm2 ( nhil, v1, 1 )
   v1(:) = v1(:) / e(0)
   v0(:) = 0.0_dp
   e(0) = 0.0_dp
   ! Lanczos procedure starts here 
   do j=1, nl
      ! v_{j-1} == v0
      ! v_j     == v1
      ! v_{j+1} = H*v_j - beta*v_{j-1}  (overwritten to v0)
      call dgemv ( 'N', nhil, nhil, 1.0_dp, H, nhil, v1, 1,-e(j-1), v0, 1 )
      ! alpha = <v_{j+1}|v_j>
      d(j) = ddot ( nhil, v0, 1, v1, 1 )
      ! v_{j+1} = v_{j+1} - alpha*v_j   (is orthogonal to v_j)
      v2(:) = v0(:) - d(j)*v1(:)
      ! beta = |v_{j+1}|
      e(j) = dnrm2 ( nhil, v2, 1 ) 
      ! v_{j+1} = v_j, v_j = v_{j+1}/beta
      v0(:) = v1(:)
      v1(:) = v2(:) / e(j) 
   end do
   deallocate ( v2, v1, v0 )

   ! Hamiltonian is now in tridiagonal form: solve

   allocate ( daux(nl), eaux(nl) )
   do i=1,nl
      ! arrays containing diagonal and subdiagonal are destroyed by dsterf
      daux(:) = d(:)
      ! note that e(0) is zero and is not copied
      eaux(:) = e(1:)
      ! dsterf calculates eigenvalues only of a tridiagonal matrix
      call dsterf( i, daux, eaux, info )
      if (info /= 0) write(6,'("Error in dsterf, info=",i4)') info
      !write(6,'("nl=",i3," E=",f15.10,",   E/N=",f15.10)') i,daux(1),daux(1)/N
      write(6,'("nl=",i3," E=",4f14.8)') i,daux(1:4)
   end do
   deallocate ( eaux, daux, e, d )
   
   ! check: conventional diagonalization
   
   lwork = 3*nhil
   allocate ( W(nhil), work(lwork) )
   call dsyev ('N','U', nhil, H, nhil, W, work, lwork, info )
   deallocate ( work )
   if (info /= 0) write(6,'("Error in dsyev, info=",i4)') info
   write(6,'("exact: E=",f15.10,",   E/N=",f15.10)') W(1),W(1)/N
   !
   deallocate ( W, states, H, neigh, conf )
   go to 10
   !
end program heisenberg_exact

!---------------------------------------------------------------
integer function search ( p, nhil, k )
!---------------------------------------------------------------
    ! On input:
    !      p(1:nhil) = labels for all states (ordered: p(i+1) > p(i) )
    !      k         = a state label
    ! On output:
    !      index i such that p(i) = k
    implicit none

    integer, intent(in) :: k, nhil
    integer, intent(in) :: p(nhil)

    integer, parameter :: dp = selected_real_kind(14, 200)
    integer :: imin, imax, lim, l, ii

    lim = log(float(nhil))/log(2.0_dp) + 1
    ! 2^lim > Nhil: lim = max number of steps needed to locate p(i)=k

    imin = 1
    imax = nhil
bisec: do l=1,lim
       ii=(imin+imax)/2
       if (p(ii) == k) then
          imin=ii
          exit bisec
       else
          if ( p(ii) > k ) then
             imax=ii-1
          else
             imin=ii+1
          endif
       endif
    enddo bisec 

    if ( ii /= imin ) then
       write(6,'("Something wrong: search not converged")')
       stop
    else
       search = ii
       return
    end if
end function search

recursive function bin(n,k) result (bun)
  ! binomial coefficient n!/k!/(n-k)! using recursion
  ! (and not using factorials that overflow for n>12)
  implicit none
  integer :: bun, n, k
  !
  if ( n < 1 .or. k < 0 .or. k > n ) then
     bun = 0
  else if ( n == 1 ) then
     bun = 1
  else
     bun = bin(n-1, k-1) + bin(n-1,k)
  end if
  !
end function bin
