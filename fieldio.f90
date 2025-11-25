#include "macros.h"

module fieldio
    use numbers
    use openmpi
    use io
    use parameters
    use fftw

    contains 

!==============================================================================

    subroutine fieldio_read(vfieldk)

        complex(dpc), intent(out) :: vfieldk(:, :, :, :)
        _indices
        integer(i4) :: forcing1, nx1, ny1, nz1, n, un, count
        TYPE(MPI_File) :: fh
        integer(MPI_OFFSET_KIND)  :: offset
        complex(dpc), allocatable, dimension(:, :, :) :: &
            buf_sfieldk, buf2_sfieldk

        inquire(file=TRIM(fname),exist=there)
        if(.not.there) then
            write(out,*) 'fieldio_read: Stopping, cannot find file : '//TRIM(fname)
            flush(out)
            error stop
        end if

        if (ix_max /= -1) then
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc - 1))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc - 1))
        else
            allocate(buf_sfieldk(ny_half, nz - 1, nx_perproc))
            allocate(buf2_sfieldk(ny_half, 1, nx_perproc))
        end if

        if (my_id==0) then
            open(newunit=un,file=TRIM(fname),form='unformatted',access='stream')
            read(un) forcing1, nx1, ny1, nz1
            close(un)
        end if

        call MPI_BCAST(forcing1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(nx1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ny1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(nz1,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

        ! if (forcing /= forcing1) then
        !     write(out,*) 'fieldio_read: Error, forcing is different.'
        !     write(out,*) 'fieldio_read:     .in file: ', forcing
        !     write(out,*) 'fieldio_read: Restart file: ', forcing1
        !     flush(out)
        !     error stop
        ! end if 

        if (nx/=nx1 .or. ny/=ny1 .or. nz/=nz1) then
            write(out,*) 'fieldio_read: Eror, grid is different.'
            write(out,*) 'fieldio_read:     .in file: ',nx,ny,nz
            write(out,*) 'fieldio_read: Restart file: ',nx1,ny1,nz1
            flush(out)
            error stop
        end if

        ! opening the file
        call MPI_INFO_CREATE(mpi_info_var, mpi_err)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(fname), MPI_MODE_RDONLY, mpi_info_var, &
                           fh, mpi_err)

        ! In all cases x was made to be the slowest index while saving to disk
        ! And zeroed modes in x and z were skipped

        do n = 1, 3

            offset = 68 + (n-1)*(nx-1)*ny_half*(nz-1)*16 &
                                                + my_id*nx_perproc*ny_half*(nz-1)*16
            if (my_id > my_id_ix_max) offset = offset - ny_half*(nz-1)*16
            if (ix_max /= -1) then
                count = (nx_perproc-1) * ny_half * (nz-1)
            else
                count = nx_perproc * ny_half * (nz-1)
            end if
            call MPI_FILE_READ_AT_ALL(fh, offset, buf_sfieldk, count, &
                                        MPI_COMPLEX16, mpi_status_var, mpi_err)

            _loop_spec_begin
                if (ix_max /= -1) then
                    if (ix < ix_max) then
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                        end if
                    else
                        if (iz < iz_max) then
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix-1)
                        else
                            vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix-1)
                        end if
                    end if
                else 
                    if (iz < iz_max) then
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz,ix)
                    else
                        vfieldk(ix,iy,iz,n) = buf_sfieldk(iy,iz-1,ix)
                    end if
                end if
            _loop_spec_end

        end do

        deallocate(buf_sfieldk)
        deallocate(buf2_sfieldk)

        if (ix_max /= -1) vfieldk(ix_max,:,:,1:3) = 0
        vfieldk(:,:,iz_max,1:3) = 0

        call MPI_FILE_CLOSE(fh, mpi_err)
        call MPI_INFO_FREE(mpi_info_var, mpi_err)


    end subroutine fieldio_read

!==============================================================================

end module fieldio
