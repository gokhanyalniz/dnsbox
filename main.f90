program main
    ! modules:
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use fieldio
    use rhs
    use vfield
    use timestep
    use stats
    use run
    use solver
    use projector
    use lyap
    use symred
    
    ! initialization:
    call run_init

    write(out, *) "Starting time stepping."

    do
        if (i_finish > 0 .and. itime > i_finish) exit

        if (i_save_fields > 0 .and. itime > i_start .and. mod(itime, i_save_fields) == 0) then
            write(file_ext, "(i6.6)") itime/i_save_fields
            fname = 'state.'//file_ext
            call fieldio_write(vel_vfieldk_now)
        end if

        if (i_print_stats > 0 .and. mod(itime, i_print_stats) == 0) then
            call stats_compute(vel_vfieldk_now, fvel_vfieldk_now)
        end if

        ! Write stats
        if (i_print_stats > 0 .and. mod(itime, i_print_stats) == 0) then
            call stats_write
        end if

        if (i_print_steps > 0 .and. mod(itime, i_print_steps) == 0) call timestep_write

        ! wall clock limit
        call cpu_time(cput_now)
        if (wall_clock_limit > 0 .and. cput_now - cput_start > wall_clock_limit) then
            write(out, *) "Runtime limit reached, stopping."
            kill_switch = .true.
        end if

        ! Broadcast the kill_switch status
        call MPI_BCAST(kill_switch, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)

        if (kill_switch) then
            call run_exit
        end if

        call timestep_precorr(vel_vfieldxx_now, vel_vfieldk_now, fvel_vfieldk_now)

        time = time + dt
        itime = itime + 1

    end do

    write(out, *) "Finished time stepping."
    call run_exit
    
end program main

