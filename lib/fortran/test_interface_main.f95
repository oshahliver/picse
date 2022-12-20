MODULE interface

    use my_module
    use class_dict

    contains

    SUBROUTINE do_some_science_stuff(a, b, c, d, e, f, x, y, z)

        implicit none

        real(8), intent(in) :: a, b, c
        real(8), intent(out) :: d, e, f
        real(8), intent(out) :: x, y, z
        
        call simulation(a=a, b=b, c=c, d=d, e=e, f=f)

        ! Compute some more paramters from the outputs
        x = a + d
        y = b + e
        z = c + f

    END SUBROUTINE do_some_science_stuff

END module interface