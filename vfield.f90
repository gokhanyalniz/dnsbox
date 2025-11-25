#include "macros.h"
#include "mt19937_64.f90"
module vfield
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio
    use fftw
    use diffops
    use symmops
    
    logical     :: seed_init = .false.
    integer(i8) :: seed
    complex(dpc), allocatable, dimension(:, :, :, :) :: laminar_vfieldk ! Laminar solution
    
    contains 

!==============================================================================

    subroutine vfield_init

        allocate(laminar_vfieldk(nx_perproc, ny_half, nz, 3))
        call vfield_laminar(laminar_vfieldk)

    end subroutine vfield_init

!==============================================================================

    subroutine vfield_inprod(vfieldk1, vfieldk2, res, allreduce)
        ! res = <vfieldk1, vfieldk2>, returns only the real part
        ! of the result
        complex(dpc), intent(in) :: vfieldk1(:, :, :, :)
        complex(dpc), intent(in) :: vfieldk2(:, :, :, :)
        real(dp), intent(out) :: res
        logical, intent(in) :: allreduce
        real(dp)          :: res1
        complex(dpc)      :: dummy
        _indices
        integer(i4) :: n

        res1 = 0
        do n=1,3
            _loop_spec_begin
                dummy = conjg(vfieldk1(ix,iy,iz,n)) * vfieldk2(ix,iy,iz,n)
                ! correct for double counting
                if (iy == 1) dummy = dummy * 0.5_dp
                res1 = res1 + dummy%re
            _loop_spec_end
        end do

        if (.not. allreduce) then
            call MPI_REDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
        else
            call MPI_ALLREDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, &
            MPI_COMM_WORLD, mpi_err)
        end if
        
    end subroutine vfield_inprod

!==============================================================================

    subroutine vfield_enstrophy(vfieldk, res, allreduce)

        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        real(dp), intent(out) :: res
        logical, intent(in) :: allreduce
        real(dp)          :: res1
        complex(dpc)      :: dummy
        _indices
        integer(i4) :: n

        res1 = 0
        do n=1,3
            _loop_spec_begin
                dummy = laplacian(ix,iy,iz) * &
                            conjg(vfieldk(ix,iy,iz,n)) * vfieldk(ix,iy,iz,n)
                ! correct for double counting
                if (iy == 1) dummy = dummy * 0.5_dp
                res1 = res1 + dummy%re
            _loop_spec_end
        end do

        if (.not. allreduce) then
            call MPI_REDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
        else
            call MPI_ALLREDUCE(res1, res, 1, MPI_REAL8, MPI_SUM, &
            MPI_COMM_WORLD, mpi_err)
        end if
        
    end subroutine vfield_enstrophy

!==============================================================================

    subroutine vfield_norm2(vfieldk, res, allreduce)
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        logical, intent(in) :: allreduce
        real(dp), intent(out) :: res

        call vfield_inprod(vfieldk, vfieldk, res, allreduce)
        
    end subroutine vfield_norm2

!==============================================================================

    subroutine vfield_norm(vfieldk, res, allreduce)
        complex(dpc), intent(in) :: vfieldk(:, :, :, :)
        logical, intent(in) :: allreduce
        real(dp), intent(out) :: res
        
        call vfield_norm2(vfieldk, res, allreduce)
        res = sqrt(res)
        
    end subroutine vfield_norm

!==============================================================================

    subroutine vfield_laminar(vfieldk)
        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        
        vfieldk(:,:,:,1:3) = 0

        if (ix_zero /= -1 .and. iy_force /= -1) then
            if (forcing == 1) then
                vfieldk(ix_zero,iy_force,1,1)%im = -0.5_dp
            elseif (forcing == 2) then
                vfieldk(ix_zero,iy_force,1,1)%re = 0.5_dp
            endif
        endif
    end subroutine vfield_laminar

!==============================================================================

    subroutine vfield_pressure(vfieldk)
        ! Apply pressure solver on vfieldk(:,:,:,1:3) 
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)
        _indices
        integer(i4)  :: d
        complex(dpc) :: dpressure
        
        _loop_spec_begin
            dpressure = 0
            do d=1,3
                dpressure = dpressure + vfieldk(ix,iy,iz,d) * &
                                        vfield_coordinatek(ix,iy,iz,d) * &
                                        inverse_laplacian(ix,iy,iz)
            end do
            
            do d=1,3
                vfieldk(ix,iy,iz,d) = vfieldk(ix,iy,iz,d) - dpressure * &
                                                vfield_coordinatek(ix,iy,iz,d)
            end do
        _loop_spec_end
                
    end subroutine vfield_pressure

!==============================================================================

    subroutine vfield_galinv(vfieldk)
        ! Set (0,0,0)-mode to 0 (Galilean invariance)
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)

        ! index ix_zero has the zero mode in x
        if(ix_zero /= -1) vfieldk(ix_zero,1,1,1:3) = 0
    end subroutine vfield_galinv

!==============================================================================

    subroutine vfield_solvediv(vfieldk)
        complex(dpc), intent(inout) :: vfieldk(:, :, :, :)

        _indices

        ! uses
        !   u k_x + v k_y + w k_z = 0
        ! to solve for one velocity component in terms of the other two

        if (nx - 1 >= ny_half .and. nx - 1 >= nz - 1) then
            _loop_spec_begin
                if (ix_zero /= -1 .and. ix==ix_zero) cycle
                vfieldk(ix,iy,iz,1) = -(vfieldk(ix,iy,iz,2)*ky(iy) &
                                            + vfieldk(ix,iy,iz,3)*kz(iz)) / kx(ix)
            _loop_spec_end
        elseif (ny_half >= nx - 1 .and. ny_half >= nz - 1) then
            _loop_spec_begin
                if (iy==1) cycle
                vfieldk(ix,iy,iz,2) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                            + vfieldk(ix,iy,iz,3)*kz(iz)) / ky(iy)
            _loop_spec_end
        else
            _loop_spec_begin
                if (iz==1) cycle
                vfieldk(ix,iy,iz,3) = -(vfieldk(ix,iy,iz,1)*kx(ix) &
                                                + vfieldk(ix,iy,iz,2)*ky(iy)) / kz(iz)
            _loop_spec_end
        end if

    end subroutine vfield_solvediv

!==============================================================================

    real(dp) function vfield_coordinatek(ix,iy,iz,d)
        integer(i4), intent(in) :: ix, iy, iz, d
        select case (d)
            case(1)
                vfield_coordinatek = kx(ix)
            case(2)
                vfield_coordinatek = ky(iy)
            case(3)
                vfield_coordinatek = kz(iz)
        end select
    end function vfield_coordinatek

!==============================================================================

end module vfield
