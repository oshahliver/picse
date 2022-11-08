MODULE test

use class_table

implicit none

!Why do I need to do it like this? Why can I not just
!Define type(table_list) :: tables ??
private
public :: initiate, compute
type(table_list), public :: tables

contains 

SUBROUTINE initiate(n_tables)

implicit none

integer, intent(in) :: n_tables
integer :: ord=1
integer :: i,j

allocate(tables%objects(n_tables))
allocate(tables%n_in_list(n_tables))

do i=1, n_tables
  call init(self=tables%objects(i), label=i)
enddo

print *, 'We have the following objects:'
do i=1, n_tables
  write(*,*) tables%objects(i)%label
enddo


do i=1, n_tables
  call load(tables%objects(i), file_name='table.tab')
  print *, 'The axes_lengths are:'
  print *, tables%objects(i)%axes_lengths
enddo


END SUBROUTINE initiate

END MODULE test
