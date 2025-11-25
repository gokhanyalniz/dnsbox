#include "macros.h"
module rhs
    use numbers
    use openmpi
    use io
    use parameters
    use fieldio 
    use fftw
    use vfield

    contains

!==============================================================================

    subroutine rhs_nonlin_term(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)

        real(dp), intent(in)  :: vel_vfieldxx(:, :, :, :)
        complex(dpc), intent(in) :: vel_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: fvel_vfieldk(:, :, :, :)

        complex(dpc) :: rhs_vfieldk(nx_perproc, ny_half, nz, 6), advect(3), advect_, div
        real(dp) :: rhs_vfieldxx(nyy, nzz_perproc, nxx, 5), u_out_u(6), trace

        integer(i4) :: n, i, j

        _indices
        _indicess
    
        ! get 6 velocity products:
        _loop_phys_begin

            do n = 1, 6
                u_out_u(n) = vel_vfieldxx(iyy,izz,ixx,isym(n)) &
                                * vel_vfieldxx(iyy,izz,ixx,jsym(n))
            end do

            ! Basdevant 1983, subtract the trace
            trace = u_out_u(nsym(1,1)) + u_out_u(nsym(2,2)) &
                                                + u_out_u(nsym(3,3))
            u_out_u(nsym(1,1)) = u_out_u(nsym(1,1)) - trace / 3.0_dp
            u_out_u(nsym(2,2)) = u_out_u(nsym(2,2)) - trace / 3.0_dp
            ! No need to do
            !   u_out_u(nsym(3,3)) = u_out_u(nsym(3,3)) - trace / 3.0_dp
            ! as it doesn't get used.
            
            ! No need to write to rhs_vfieldxx(:,:,:,6) as it is not used.
            do n = 1, 5
                rhs_vfieldxx(iyy,izz,ixx,n) = u_out_u(n)
            end do
        _loop_phys_end

        do n = 1,5
            call fftw_sx2k(rhs_vfieldxx(:,:,:,n), rhs_vfieldk(:,:,:,n))
        end do

        ! Get one element on the diagonal from the tracelessness
        rhs_vfieldk(:,:,:,nsym(3,3)) = &
            -(rhs_vfieldk(:,:,:,nsym(1,1)) + rhs_vfieldk(:,:,:,nsym(2,2)))
        
        ! Now rhs_vfieldk(:,:,:,:, 1:6) holds the products uu, uv, uw, ...
        ! in Fourier space.  

        _loop_spec_begin                  
            ! Nonlinear term:

            ! advect(j) = - F{u_i d_i u_j}
            do j = 1, 3
                advect(j) = 0

                do i = 1 ,3
                    advect(j) = advect(j) - nabla(ix, iy, iz, i) &
                                * rhs_vfieldk(ix, iy, iz, nsym(i,j))
                end do

            end do

            advect_ = 0

            ! Pressure terms
            div = 0
            do n = 1, 3
                div = div + vfield_coordinatek(ix,iy,iz,n) * advect(n)
            end do
            div = div + advect_

            do n = 1, 3
                fvel_vfieldk(ix,iy,iz,n) = advect(n) &
                                    - vfield_coordinatek(ix,iy,iz,n) * div &
                                        * inverse_laplacian(ix,iy,iz)
            end do

        _loop_spec_end
        
        ! Add forcing
        if (forcing /= 0 .and. ix_zero /= -1) then
            if (forcing == 1) then ! sine
                fvel_vfieldk(ix_zero,iy_force,1,1) = fvel_vfieldk(ix_zero,iy_force,1,1) &
                        - imag_1 * 0.5_dp * amp / (4.0_dp*Re)
            elseif (forcing == 2) then ! cosine
                fvel_vfieldk(ix_zero,iy_force,1,1) = fvel_vfieldk(ix_zero,iy_force,1,1) &
                        + 0.5_dp * amp / (4.0_dp*Re)
            end if      
        end if

    end subroutine rhs_nonlin_term

!==============================================================================

    subroutine rhs_all(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)
        real(dp), intent(in)  :: vel_vfieldxx(:, :, :, :)
        complex(dpc), intent(in)  :: vel_vfieldk(:, :, :, :)
        complex(dpc), intent(out) :: fvel_vfieldk(:, :, :, :)

        _indices

        call rhs_nonlin_term(vel_vfieldxx, vel_vfieldk, fvel_vfieldk)

        _loop_spec_begin
            fvel_vfieldk(ix,iy,iz,1:3) = fvel_vfieldk(ix,iy,iz,1:3) &
                                    + (- laplacian(ix,iy,iz) / Re) & 
                                    * vel_vfieldk(ix,iy,iz,1:3)
        _loop_spec_end

    end subroutine rhs_all

!==============================================================================

end module rhs
