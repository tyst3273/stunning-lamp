! --------------------------------------------------------------------------------------------------
module m_template
! --------------------------------------------------------------------------------------------------
    !
    !   name: m_template
    !
    !   whats it do:
    !       this is just an example placeholder. it shows all of the things that should go
    !       into defining a new module file and subroutine.
    !
    !   subroutines:
    !       sub_1
    !
    !   authors:
    !       Ty Sterling
    !
    ! ----------------------------------------------------------------------------------------------
    !
    !   import external modules here
    !
    use m_external1                             ! this brings the whole thing in
    use m_external2, only   ::  extern_subrt    ! this only brings in extern_subrt
    !
    !   parameters, etc
    !
    implicit none
    !
    !   default everything to private, set the ones we want to be public to public
    !
    private
    !
    public  :: sub_1    ! an example 
    public  :: sub_2    ! another examples
    !
    !   all subroutines etc. go below here
    !
    contains
    !
    ! ----------------------------------------------------------------------------------------------
    subroutine sub_1( arg1, arg2, arg3 )
    ! ----------------------------------------------------------------------------------------------
        !
        !   name: sub_1
        !
        !   whats it do:
        !       this is just an example placeholder. it shows all of the things that should go 
        !       into defining a new module file and subroutine.
        !
        !   input:
        !       arg1: integer, example input arg
        !
        !   in/out (modified in place):
        !       arg2: real, example modified in place arg   
        !
        !   output:
        !       arg3: real, example arg to be 'returned'
        !
        !   authors:
        !       Ty Sterling
        !
        ! ----------------------------------------------------------------------------------------------
        !
        !   input args
        !
        integer, intent(in) :: arg1
        !
        !   in place args
        !
        real, intent(in,out) :: arg2
        !
        !   output args
        !
        real, intent(out) :: arg3
        !
        !   local/work variables
        !
        integer :: ii, jj, kk
        real, allocatable :: vec3d(:,:,:)
        !
        !   put what ever else it needs to do below here
        !
    ! ----------------------------------------------------------------------------------------------
    end sub_1
    ! ----------------------------------------------------------------------------------------------
    !
! --------------------------------------------------------------------------------------------------
end module m_template
! --------------------------------------------------------------------------------------------------
