program test_reshape

   implicit none

   real(8), dimension(2, 2) :: arr = reshape((/1d0, 3d0, 2d0, 4d0/), (/2, 2/))
   print *, 'arr =', arr(1, :)
   print *, 'arr =', arr(2, :)

end program test_reshape
