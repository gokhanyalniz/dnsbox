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

    integer :: i_benchmark_iter = 0
    real(sp) :: bench_0, bench_t
    
    ! initialization:
    call run_init

    write(out, *) "Starting time stepping."

    do
        if (i_finish > 0 .and. itime > i_finish) exit

        if (i_save_fields > 0 .and. itime > i_start .and. mod(itime, i_save_fields) == 0) then
            write(file_ext, "(i6.6)") itime/i_save_fields
            fname = 'state.'//file_ext
            call fieldio_write(vel_vfieldk_now)
            call run_flush_channels
        end if

        if (i_save_phys > 0 .and. itime > i_start .and. mod(itime, i_save_phys) == 0) then
            write(file_ext, "(i6.6)") itime/i_save_phys
            fname = 'phys.'//file_ext
            ! call fftw_vk2x_0(vel_vfieldk_now, vel_vfieldx_now)
            call fieldio_write_phys(vel_vfieldxx_now)
            call run_flush_channels
        end if

        if (adaptive_dt .or. (i_print_steps > 0 .and. mod(itime, i_print_steps) == 0)) then
            call timestep_courant(vel_vfieldxx_now)
        end if

        if (adaptive_dt) call timestep_set_dt

        if (i_print_stats > 0 .and. mod(itime, i_print_stats) == 0) then
            
            call stats_compute(vel_vfieldk_now, fvel_vfieldk_now)

        end if

        ! Write stats
        if (i_print_stats > 0 .and. mod(itime, i_print_stats) == 0) then
            call stats_write

            ! Stop if laminarized
            if (my_id==0 .and. terminate_laminar) then

                e_diff = abs(ekin - ekin_lam)/ekin_lam
                input_diff = abs(powerin - powerin_lam)/powerin_lam
                diss_diff = abs(dissip - dissip_lam )/dissip_lam 
                
                if (e_diff < relerr_lam .and. &
                    input_diff < relerr_lam .and. diss_diff < relerr_lam) then
                    kill_switch = .true.
                    write(out, *) "Laminarized, stopping."
                    open(newunit=laminarized_ch,file='LAMINARIZED',position='append')
                    write(laminarized_ch,*) time
                    close(laminarized_ch)
                end if
            end if
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

        ! Flush manually if desired
        if (i_flush > 0 .and. mod(itime, i_flush) == 0) then
            call run_flush_channels
        end if

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

