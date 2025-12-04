#include "macros.h"
module run
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use diffops
    use symmops
    use vfield
    use timestep
    use stats
    use projector
    use lyap
    use symred

    real(sp) :: cput_now, cput_start, cput_stop

    complex(dpc), allocatable, dimension(:, :, :, :) :: &
        vel_vfieldk_now, fvel_vfieldk_now! u and F(u) now

    real(dp), allocatable :: vel_vfieldxx_now(:, :, :, :)

    logical     :: kill_switch = .false.

    contains
    
    subroutine run_init
     
        call openmpi_init
        call io_init
        call parameters_init
        call fftw_init
        call vfield_init

        allocate(vel_vfieldk_now(nx_perproc, ny_half, nz, 3))
        allocate(fvel_vfieldk_now(nx_perproc, ny_half, nz, 3))
        allocate(vel_vfieldxx_now(nyy, nzz_perproc, nxx, 3))

        ! Initial time
        itime = i_start
        time  = t_start

        write(file_ext, "(i6.6)") IC
        fname = 'state.'//file_ext
        call fieldio_read(vel_vfieldk_now)

        call fftw_vk2x(vel_vfieldk_now, vel_vfieldxx_now)

        call rhs_all(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)

        call cpu_time(cput_start)

    end subroutine run_init

!==============================================================================

    subroutine run_exit
        real(dp)    :: run_time_wall, dt_sim
        integer(i4) :: num_time_step_sim

        call run_close_channels
        call cpu_time(cput_stop)

        run_time_wall =  cput_stop - cput_start
        dt_sim = time - t_start
        num_time_step_sim = itime - i_start
        write(out, *) "sec:", run_time_wall
        if (dt_sim > 0) write(out, *) "sec/t:", run_time_wall / dt_sim
        if (num_time_step_sim > 0) write(out, *) "sec/step:", run_time_wall / num_time_step_sim
        write(out, *) "cpusec:", run_time_wall * num_procs
        if (dt_sim > 0) write(out, *) "cpusec/t:", run_time_wall * num_procs / dt_sim
        if (num_time_step_sim > 0) write(out, *) "cpusec/step:", run_time_wall * num_procs / num_time_step_sim
        if (run_time_wall > 0) then
            write(out, *) "fft/sec:", time_fftw / run_time_wall
            write(out, *) "local_transpose/sec:", time_local_transpose / run_time_wall
            write(out, *) "global_transpose/sec:", time_global_transpose / run_time_wall
        end if

        call io_exit    
        call openmpi_exit
        stop

    end subroutine run_exit

!==============================================================================

    subroutine run_close_channels
        if (stats_stat_written) close(stats_stat_ch)
        if (stats_frac_written) close(stats_frac_ch)
        if (steps_written) close(steps_ch)
    end subroutine run_close_channels

!==============================================================================

end module run