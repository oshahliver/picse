module mesh

implicit none
type :: geometry
    real(8) :: xmin, xmax, dx            ! coordinates of origin and grid size
    integer :: nx                        ! number of grid points
    real(8), dimension(:), pointer :: x  ! coordinates of points
end type geometry

contains

subroutine create(geom, xmin, xmax, nx)

    !f2py integer(8), intent(out) :: geom
    type(geometry), pointer :: geom
    real(8), intent(in) :: xmin, xmax
    integer, intent(in) :: nx
            
    real(8) :: dx
            
    integer :: i
            
    allocate(geom)
    geom%xmin = xmin
    geom%xmax = xmax
    geom%dx = ( xmax - xmin ) / (nx-1) 
    geom%nx = nx
    allocate(geom%x(nx))
    do i=1,nx
        geom%x(i)=geom%xmin+(i-1)*geom%dx
    end do

end subroutine create

subroutine view(geom)
    !f2py integer(8), intent(in) :: geom
    type(geometry), pointer :: geom
    print*, 'nx = ', geom%nx
    print*, geom%xmin, geom%xmax
    print*, geom%x(:)
end subroutine view

subroutine get_size(geom, nx)

    !f2py integer(8), intent(in) :: geom
    type(geometry), pointer :: geom
    integer, intent(out) :: nx
    
    nx = geom%nx
    
end subroutine get_size


subroutine create_type(geom)
    !f2py integer(8), intent(in) :: geom
    type(geometry), pointer :: geom
 

end subroutine create_type

end module mesh