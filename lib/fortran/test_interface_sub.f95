MODULE my_module

    
    contains 

    SUBROUTINE simulation(a,b,c,d,e,f)

        implicit none
        
        real(8), intent(in) :: a, b, c
        real(8), intent(out) :: d, e, f

        ! Do lot's of science stuff with the input
        ! ...

        ! Create some outputs

        d = a
        e = b
        f = c

    END SUBROUTINE simulation

END MODULE my_module