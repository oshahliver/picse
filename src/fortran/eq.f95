PROGRAM eq

implicit none

logical :: trans = .false.
integer :: i, ph1, ph2

ph1 = 0
ph2 = 0

read (*,*) ph1
read (*,*) ph2

if(ph1.eq.ph2)then
	print*, 'Equal'
else
	print *, 'Not equal'
endif

END PROGRAM eq
