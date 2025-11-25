#include "macros.h"

module diffops
    use numbers
    use openmpi
    use io
    use parameters
    use fftw

    contains 

!==============================================================================

    complex(dpc) function nabla(ix, iy, iz, d)
        integer(i4), intent(in) :: ix, iy, iz, d

        select case (d)
            case(1)
                nabla = imag_1 * kx(ix)
            case(2)
                nabla = imag_1 * ky(iy)
            case(3)
                nabla = imag_1 * kz(iz)
        end select
    end function nabla

!==============================================================================

    complex(dpc) function laplacian(ix, iy, iz)
        integer(i4), intent(in) :: ix, iy, iz

        laplacian = kx(ix)**2 + ky(iy)**2 + kz(iz)**2
    end function laplacian

!==============================================================================

    complex(dpc) function inverse_laplacian(ix, iy, iz)
        integer(i4), intent(in) :: ix, iy, iz
        
        if (qx(ix) == 0 .and. qy(iy) == 0 .and. qz(iz) == 0) then
            inverse_laplacian = 0
        else
            inverse_laplacian = 1.0_dp / laplacian(ix,iy,iz)
        end if
    end function inverse_laplacian

!==============================================================================  

end module diffops