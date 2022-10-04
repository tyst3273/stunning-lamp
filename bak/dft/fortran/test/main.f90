
program main

    use m_template

    implicit none
    integer :: ii, jj, kk

!    write(*,*)
!    write(*,'("give 3 integers > ")',advance='no')
!    read(*,*) ii, jj, kk

    ii = 2
    jj = 2
    kk = 2

    write(*,*)
    write(*,"('you entered: ii = ',i3,', jj = ',i3,', kk = ',i3)") ii, jj, kk
    write(*,*)

    call alloc_mat(ii,jj,kk)
    
    write(*,*) mat3d

    call fill_mat3d()
    call print_mat3d()

    call clean_up()

end program main
