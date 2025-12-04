#include "macros.h"
module stats
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use vfield
    use rhs
    use timestep

    real(dp) :: ekin, powerin, enstrophy, dissip, norm_rhs, &
                power_unit, ekin_perturb

    integer(i4) :: stats_stat_ch,
    logical :: stats_stat_written = .false.
    
    character(255) :: stats_stat_file = 'stat.gp'

    contains 

!==============================================================================

    subroutine stats_compute(vfieldk, fvfieldk)
        complex(dpc), intent(in), dimension(:, :, :, :)  :: &
             vfieldk, fvfieldk

        real(dp) :: norm2_hor

        ! Kinetic energy
        call vfield_norm2(vfieldk, ekin, .false.)

        ! Power input
        ! power_unit is the inner product with the laminar state...
        call vfield_inprod(vfieldk, laminar_vfieldk, power_unit, .false.)
        power_unit = 2.0_dp * power_unit ! get rid of the 1/2 factor
        ! ...which is proportional to the inner product with the forcing
        powerin = (amp / (4.0_dp * Re)) * power_unit

        ! Perturbation kinetic energy, (1/2)|u - u_lam|^2
        ekin_perturb = ekin + ekin_lam - power_unit
                
        ! Viscous dissipation
        call vfield_enstrophy(vfieldk, enstrophy, .false.)
        dissip = 2.0_dp * enstrophy / Re

        ! norm of rhs
        call vfield_norm(fvfieldk,norm_rhs,.false.)

    end subroutine stats_compute

!==============================================================================

    subroutine stats_write
        
        ! outputting statistics
        
        if (my_id==0) then
            
           ! outputting all this in the stat file

            inquire(file=TRIM(stats_stat_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_stat_ch,file=TRIM(stats_stat_file),form='formatted')
                write(stats_stat_ch,"(A2,"//i4_len//","//"6"//sp_len//")") &
                    "# ", "itime", "time", "ekin", "powerin", "dissip", "norm_rhs", "power_unit"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_stat_ch,file=TRIM(stats_stat_file),position='append')
            end if
            write(stats_stat_ch,"(A2,"//i4_f//","//"6"//sp_f//")")&
                "  ", itime, time, ekin, powerin, dissip, norm_rhs, power_unit

           stats_stat_written = .true.

        end if

    end subroutine stats_write

!==============================================================================   

end module stats