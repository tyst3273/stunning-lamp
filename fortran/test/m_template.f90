module m_template

    implicit none

    private
    public  :: alloc_mat, clean_up, print_mat3d, fill_mat3d
    public  :: mat3d

    real, allocatable :: mat3d(:,:,:)
    integer :: n1, n2, n3

    contains


    subroutine alloc_mat( ii, jj, kk )

        integer, intent(in) :: ii, jj, kk
        n1 = ii
        n2 = jj
        n3 = kk
        allocate( mat3d(ii,jj,kk) )
        mat3d(:,:,:) = 0

    end subroutine alloc_mat


    subroutine clean_up

        call destroy_mat()

    end subroutine clean_up


    subroutine destroy_mat

        write(*,*)
        write(*,*) ' destroy! '
        write(*,*)
        deallocate( mat3d )

    end subroutine destroy_mat


    subroutine fill_mat3d

        integer :: ii, jj, kk

        do ii = 1, n1
            do jj = 1, n2
                do kk = 1, n2
                    mat3d(ii,jj,kk) = ii+jj+kk
                end do
            end do
        end do

        write(*,*)
        write(*,*) 'done filling matrix'
        write(*,*)

    end subroutine fill_mat3d


    subroutine print_mat3d

        integer :: ii, jj, kk

        do ii = 1, n1
            do jj = 1, n2
                do kk = 1, n3
                    write(*,*) ii, jj, kk, mat3d(ii,jj,kk)
                end do
            end do
        end do

    end subroutine print_mat3d


end module m_template
