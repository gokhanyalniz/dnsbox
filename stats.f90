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
                dissip_mhd, input_mhd, dissip_ray, input_ray

    integer(i4) :: stats_stat_ch, stats_specx_ch, stats_specy_ch, &
                   stats_specz_ch, stats_stat_mhd_ch, stats_stat_ray_ch
    logical :: stats_stat_written = .false., stats_specs_written = .false., &
               stats_stat_mhd_written = .false., stats_stat_ray_written = .false.
    
    character(255) :: stats_stat_file = 'stat.gp', &
                      stats_specx_file = 'specs_x.gp', &
                      stats_specy_file = 'specs_y.gp', &
                      stats_specz_file = 'specs_z.gp', &
                      stats_stat_mhd_file = 'stat_mhd.gp', &
                      stats_stat_ray_file = 'stat_ray.gp'

    contains 

!==============================================================================

    subroutine stats_compute_powerin_unit(vfieldk, res, allreduce)
        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        real(dp), intent(out) :: res
        logical, intent(in) :: allreduce
        real(dp) :: my_powerin

        my_powerin = 0
        if (ix_zero /= -1) then
            if (tilting) then
                if (forcing == 1) then ! sine
                    my_powerin = -cos(tilt_angle * PI / 180.0_dp) &
                                    * vfieldk(ix_zero,iy_force,1,1)%im &
                                 -sin(tilt_angle * PI / 180.0_dp) &
                                    * vfieldk(ix_zero,iy_force,1,3)%im
                elseif (forcing == 2) then ! cosine
                    ! This may require thinking in the presence of drag
                    my_powerin = cos(tilt_angle * PI / 180.0_dp) &
                                                        * vfieldk(ix_zero,iy_force,1,1)%re  &
                                 + sin(tilt_angle * PI / 180.0_dp) &
                                                        * vfieldk(ix_zero,iy_force,1,3)%re
                end if
            else
                if (forcing == 1) then ! sine
                    my_powerin = -vfieldk(ix_zero,iy_force,1,1)%im
                elseif (forcing == 2) then ! cosine
                    ! This may require thinking in the presence of drag
                    my_powerin = vfieldk(ix_zero,iy_force,1,1)%re 
                end if
            end if
        end if

        if (.not. allreduce) then
            call MPI_REDUCE(my_powerin, res, 1, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
        else
            call MPI_ALLREDUCE(my_powerin, res, 1, MPI_REAL8, MPI_SUM, &
            MPI_COMM_WORLD, mpi_err)
        end if

    end subroutine stats_compute_powerin_unit

!==============================================================================

    subroutine stats_compute(vfieldk, fvfieldk, cur_vfieldk)
        complex(dpc), intent(in), dimension(:, :, :, :)  :: &
             vfieldk, fvfieldk

        complex(dpc), optional, intent(in), dimension(:, :, :, :) :: &
            cur_vfieldk

        real(dp) :: norm2_hor, power_unit

        ! Kinetic energy
        call vfield_norm2(vfieldk, ekin, .false.)

        ! Power input
        call stats_compute_powerin_unit(vfieldk, power_unit, .false.)
        powerin = (amp / (4.0_dp * Re)) * power_unit
                
        ! Viscous dissipation
        call vfield_enstrophy(vfieldk, enstrophy, .false.)
        dissip = 2.0_dp * enstrophy / Re

        if (rayleigh_friction .or. MHD) then
            call vfield_norm2_horizontal(vfieldk, norm2_hor, .false.)
        end if

        ! Input and dissipation due to Rayleigh friction
        if (rayleigh_friction) then
            input_ray = sigma_R * power_unit
            powerin = powerin + input_ray

            dissip_ray = 2 * sigma_R * norm2_hor
            dissip = dissip + dissip_ray
        end if

        if (MHD) then

            ! MHD dissipation
            dissip_mhd = 2*(Ha**2/Re) * norm2_hor
            dissip = dissip + dissip_mhd

            ! MHD total power
            call vfield_power_mhd(vfieldk, cur_vfieldk, input_mhd, .false.)
            
            ! Input due to MHD (above includes the -dissipation)
            input_mhd = input_mhd + dissip_mhd
            powerin = powerin + input_mhd

        end if

        ! norm of rhs
        call vfield_norm(fvfieldk,norm_rhs,.false.)

    end subroutine stats_compute

!==============================================================================

    subroutine stats_spectra(vfieldk)
        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        real(dp) :: specx(nx_half), my_specx(nx_half), &
                    specy(ny_half), my_specy(ny_half), &
                    specz(nz_half), my_specz(nz_half)
        complex(dpc) :: spec_
        character(255) :: formatStr

        _indices

        my_specx(:) = 0
        my_specy(:) = 0
        my_specz(:) = 0

        _loop_spec_begin
            spec_ = sum(conjg(vfieldk(ix,iy,iz,1:3))*vfieldk(ix,iy,iz,1:3))
            if (iy==1) spec_ = spec_ / 2
            my_specx(abs(qx(ix))+1) = my_specx(abs(qx(ix))+1) + spec_%re
            my_specy(abs(qy(iy))+1) = my_specy(abs(qy(iy))+1) + spec_%re
            my_specz(abs(qz(iz))+1) = my_specz(abs(qz(iz))+1) + spec_%re
        _loop_spec_end

        call MPI_REDUCE(my_specx, specx, nx_half, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(my_specy, specy, ny_half, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(my_specz, specz, nz_half, MPI_REAL8, MPI_SUM, 0, &
        MPI_COMM_WORLD, mpi_err)

        if (my_id==0) then
            
            ! outputting all this in the stat file
 
             inquire(file=TRIM(stats_specx_file), exist=there, opened=there2)
             if (.not.there) then
                open(newunit=stats_specx_ch,file=TRIM(stats_specx_file),form='formatted')
                 write(stats_specx_ch,"(A2,"//i4_len//","//"2"//sp_len//")") &
                     "# ", "itime", "time", "specs_x"
             end if
             if(there.and..not.there2) then
                open(newunit=stats_specx_ch,file=TRIM(stats_specx_file),position='append')
             end if
             write(formatStr,*) nx_half + 1
             write(stats_specx_ch,"(A2,"//i4_f//","//TRIM(formatStr)//sp_f//")")&
                 "  ", itime, time, specx

            inquire(file=TRIM(stats_specy_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_specy_ch,file=TRIM(stats_specy_file),form='formatted')
                write(stats_specy_ch,"(A2,"//i4_len//","//"2"//sp_len//")") &
                    "# ", "itime", "time", "specs_y"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_specy_ch,file=TRIM(stats_specy_file),position='append')
            end if
            write(formatStr,*) 1 + ny_half
            write(stats_specy_ch,"(A2,"//i4_f//","//TRIM(formatStr)//sp_f//")")&
                "  ", itime, time, specy

            inquire(file=TRIM(stats_specz_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_specz_ch,file=TRIM(stats_specz_file),form='formatted')
                write(stats_specz_ch,"(A2,"//i4_len//","//"2"//sp_len//")") &
                    "# ", "itime", "time", "specs_z"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_specz_ch,file=TRIM(stats_specz_file),position='append')
            end if
            write(formatStr,*) 1 + nz_half
            write(stats_specz_ch,"(A2,"//i4_f//","//TRIM(formatStr)//sp_f//")")&
                "  ", itime, time, specz
 
            stats_specs_written = .true.
 
         end if
    end subroutine stats_spectra

!==============================================================================

    subroutine stats_write
        
        ! outputting statistics
        
        if (my_id==0) then
            
           ! outputting all this in the stat file

            inquire(file=TRIM(stats_stat_file), exist=there, opened=there2)
            if (.not.there) then
            open(newunit=stats_stat_ch,file=TRIM(stats_stat_file),form='formatted')
                write(stats_stat_ch,"(A2,"//i4_len//","//"5"//sp_len//")") &
                    "# ", "itime", "time", "ekin", "powerin", "dissip", "norm_rhs"
            end if
            if(there.and..not.there2) then
            open(newunit=stats_stat_ch,file=TRIM(stats_stat_file),position='append')
            end if
            write(stats_stat_ch,"(A2,"//i4_f//","//"5"//sp_f//")")&
                "  ", itime, time, ekin, powerin, dissip, norm_rhs

           stats_stat_written = .true.

            if (rayleigh_friction) then
                inquire(file=TRIM(stats_stat_ray_file), exist=there, opened=there2)
                if (.not.there) then
                open(newunit=stats_stat_ray_ch,file=TRIM(stats_stat_ray_file),form='formatted')
                    write(stats_stat_ray_ch,"(A2,"//i4_len//","//"3"//sp_len//")") &
                        "# ", "itime", "time", "input_ray", "dissip_ray"
                end if
                if(there.and..not.there2) then
                open(newunit=stats_stat_ray_ch,file=TRIM(stats_stat_ray_file),position='append')
                end if
                write(stats_stat_ray_ch,"(A2,"//i4_f//","//"3"//sp_f//")")&
                    "  ", itime, time, input_ray, dissip_ray

                stats_stat_ray_written = .true.
            end if

           if (MHD) then
                inquire(file=TRIM(stats_stat_mhd_file), exist=there, opened=there2)
                if (.not.there) then
                open(newunit=stats_stat_mhd_ch,file=TRIM(stats_stat_mhd_file),form='formatted')
                    write(stats_stat_mhd_ch,"(A2,"//i4_len//","//"3"//sp_len//")") &
                        "# ", "itime", "time", "input_mhd", "dissip_mhd"
                end if
                if(there.and..not.there2) then
                open(newunit=stats_stat_mhd_ch,file=TRIM(stats_stat_mhd_file),position='append')
                end if
                write(stats_stat_mhd_ch,"(A2,"//i4_f//","//"3"//sp_f//")")&
                    "  ", itime, time, input_mhd, dissip_mhd

                stats_stat_mhd_written = .true.
           end if

        end if

    end subroutine stats_write

!==============================================================================

    subroutine stats_worst_divergence(vfieldk)
        complex(dpc), intent(in)  :: vfieldk(:, :, :, :)
        complex(dpc) :: div_vfieldk(nx_perproc, ny_half, nz)
        real(dp)     :: div_sfieldxx(nyy, nzz_perproc, nxx)
        real(dp)     :: div, my_div
    
        call diffops_div(vfieldk, div_vfieldk)
        call fftw_sk2x(div_vfieldk, div_sfieldxx)
    
        my_div = maxval(abs(div_sfieldxx(:, :, :)))
        call MPI_REDUCE(my_div, div, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)
    
        if (my_id == 0 .and. div >  divergence_th) then
            write(out,*) 'stats: Time = ', time
            write(out,*) 'stats: Worst divergence: ', div
        end if
            
    end subroutine stats_worst_divergence    

end module stats