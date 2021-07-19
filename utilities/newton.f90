#include "../macros.h"
module mnewton
    use numbers
    use openmpi
    use io
    use run
    use solver

    contains

    include "NewtonHook.f90"
    include "GMRESm.f90"

    subroutine getrhs(x, y)

        ! function to be minimised

        real(dp), intent(in)  :: x(nnewt)
        real(dp), intent(out) :: y(nnewt)
        real(dp) :: y_(nnewt)
        integer(i4) :: ims, ims_

        do ims = 0, ms -1
            ims_ = modulo(ims+1,ms)
            call solver_steporbit(ndtss(ims+1), ims, x)
            if (ims == ms -1) call solver_relative_symmetries_apply(vel_vfieldxx_now, vel_vfieldk_now)
            call solver_vectorize(vel_vfieldk_now, ims_, y_)
        end do
        y = y_ - x            ! diff
        
        if (my_id == 0 .and. nscalars > 0) then
            do ims=0, ms-1
                y(ims*nnewt_pershot+1:ims*nnewt_pershot+nscalars) = 0 ! constraints, rhs=0
            end do
        end if
        
    end subroutine getrhs

!==============================================================================

    subroutine multJ(x, y)
        
        !  Jacobian of function + lhs of constraints on update
        
        real(dp), intent(in)     :: x(nnewt)
        real(dp), intent(out)    :: y(nnewt)   
        real(dp) :: eps, s(nnewt), d, dt_new, dt_old
        integer(i4)  :: ims
                        ! (F(x0+eps.x)-F(x0))/eps
        
        write(out, *) 'Jacobian call '
        
        eps = sqrt(solver_dotprod_ms(x,x))
        
        write(out, *) 'eps = ', eps
        if(abs(eps) < small)  then
            write(out,*) 'multJ: eps=0 (1)'
            flush(out)
            stop
        end if
        eps = epsJ * sqrt(solver_dotprod_ms(new_x,new_x)) / eps
        write(out, *) 'eps_ = ', eps
        if(abs(eps) < small)  then 
            write(out,*) 'multJ: eps=0 (2)'
            flush(out)
            stop
        end if
        y = new_x + eps*x
        ! This is for debugging
        d = solver_dotprod_ms(y, y)
        write(out, *) 'Jacobian: normy2 = ', d
        
        call getrhs(y, s)
        call getrhs(new_x, new_fx)
        y = (s - new_fx) / eps

        do ims = 0, ms -1
            if (my_id == 0 .and. find_period) then 
                dt_new = new_x(ims*nnewt_pershot + 1) * scaleT / ndtss(ims + 1)
            end if
            
            if (find_period) then
                call MPI_BCAST(dt_new, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
                write(out, *) 'MultJ: dt_new = ', dt_new

                dt_old = dt
                dt = dt_new
            end if

            ! contstraint, 
            ! no update in trajectory direction

            call solver_steporbit(1,ims,new_x)

            call solver_vectorize(vel_vfieldk_now, ims, s)
            s(ims*nnewt_pershot+1:(ims+1)*nnewt_pershot) = &
                (s(ims*nnewt_pershot+1:(ims+1)*nnewt_pershot) - new_x(ims*nnewt_pershot+1:(ims+1)*nnewt_pershot)) / dt
            d = solver_dotprod(ims,s,x)

            if (find_period) dt = dt_old
            
            if (my_id == 0 .and. nscalars > 0) then
                if (find_period) y(ims*nnewt_pershot + 1) = d
            end if
        end do
        
    end subroutine multJ

!==============================================================================

    subroutine saveorbit
        
        real(dp) :: norm_x
        integer(i4) :: KILL_SWITCH_PERIOD = 0,un, ims
        character*1 :: ims_str
        complex(dpc) :: vfieldk(nx_perproc, ny_half, nz, 3)

        ndts = sum(ndtss)
        norm_x = sqrt(solver_dotprod_ms(new_x,new_x))
        
        if (my_id == 0) then
            open(newunit=un,status='unknown',access='append',file='newton.dat')
            if(new_nits==0) then
                write(un,"(A2,"//"4"//i4_f//")") "# ", ndts, mgmres, nnewt, ms
                write(un, "(A2,"//"2"//i4_len//","//"5"//sp_len//")") "# ", "nits", "gits", &
                                            "rel_err", "tol_ratio", "del", "tol", "norm_x"
            end if
            write(un, "(A2,"//"2"//i4_f//","//"5"//sp_f//")") "  ", new_nits, new_gits, &
                                new_tol / norm_x, new_tol / tol, new_del, new_tol, norm_x
            close(un)
        end if
        
        if (my_id == 0 .and. nscalars > 0) then
            if (find_period) then
                period = 0
                do ims = 0, ms -1
                    period = period + new_x(ims*nnewt_pershot+1) * scaleT
                    ! Kill switch for negative period guesses
                    if (new_x(ims*nnewt_pershot+1) * scaleT < 0) then
                        KILL_SWITCH_PERIOD = 1
                    end if
                end do
            end if

            open(newunit=un,status='unknown',access='append',file='guesses.dat')
            if (ms > 1) then
                write(un,"(A2,"//dp_f//","//i4_f//")") "# ", period, ndts
                do ims = 0, ms -1
                    write(un, "("//i4_f//","//dp_f//","//i4_f//")") &
                            new_nits, new_x(ims*nnewt_pershot+1) * scaleT, ndtss(ims+1)
                end do
            else
                if(new_nits==0) write(un,"(A2,"//i4_f//","//i4_f//","//i4_f//")") "# ", ndts
                write(un, "("//i4_f//","//dp_f//","//dp_f//")") new_nits, period
            end if
            close(un)

        end if

        ! Broadcast the KILL_SWITCH status
        call MPI_BCAST(KILL_SWITCH_PERIOD, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

        if (KILL_SWITCH_PERIOD == 1) then
            write(out, *) "Period guess is negative, stopping."
            call newton_signal_not_converged
            call run_exit
        end if
        
        ! quit run if already converged and do not start again
        if (new_nits == 0 .and. new_tol / tol < 1) then
            call newton_done
            call run_exit
        end if

        ! Save the state file:
        write(file_ext, "(i6.6)") new_nits
        do ims = 0, ms - 1
            write(ims_str, "(i1.1)") ims
            if (ms > 1) then
                fname = 'newton.'//file_ext//'-'//ims_str
            else
                fname = 'newton.'//file_ext
            end if
            call solver_tensorize(vfieldk, ims, new_x)
            call fieldio_write(vfieldk)
        end do

    end subroutine saveorbit

!==============================================================================

    subroutine newton_done
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='NEWTON_DONE',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine newton_done

!==============================================================================

    subroutine newton_signal_converged
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='NEWTON_CONVERGED',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine newton_signal_converged

!==============================================================================

    subroutine newton_signal_not_converged
        
        integer(i4) :: un
        if (my_id == 0) then
            open(newunit=un,file='NEWTON_NOT_CONVERGED',position='append')
            write(un,*)
            close(un)
        end if
    end subroutine newton_signal_not_converged

end module mnewton

!==============================================================================

program newton
    !*************************************************************************
    !  Example with one extra parameter T.
    !  Can put parameter and constraint = 0 when not required.
    !- - - - - - - - - - - - - - -
    !  Newton vector:
    !    x(1)   = T / scaleT
    !    x(2:) = state vector x
    !
    !  Extra constraint:
    !    (F(x)-x). dx/dt = 0 .
    !       no update along direction of trajectory.
    !
    !  Jacobian approximation:
    !    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
    !
    !*************************************************************************
    
    ! modules:
    use numbers
    use openmpi
    use io
    use parameters
    use fftw
    use symmops
    use fieldio
    use vfield
    use timestep
    use run
    use solver
    use mnewton
    
    real(dp)           :: d
    logical            :: forbidNewton
    
    integer(i4) :: ims, info
    logical     :: fexist
    character*1 :: ims_str

    call run_init
    IC = 0
    adaptive_dt = .false.
    
    inquire(file='NEWTON_NOT_CONVERGED', exist=forbidNewton)
    if (forbidNewton) call run_exit

    inquire(file='NEWTON_DONE', exist=forbidNewton)
    if (forbidNewton) call run_exit

    call solver_data_read
    call solver_set_problem_size
    
    ! allocate vectors:
    allocate(new_x(nnewt))
    allocate(new_fx(nnewt))
    new_x = 0
    
    do ims = 1, ms
        if (ndtss(ims) < 0) ndtss(ims) = nint(periods(ims)/dt)
    end do
    
    if (my_id == 0 .and. nscalars > 0) then
        if (find_period) then
            do ims = 0, ms - 1
                new_x(ims*nnewt_pershot + 1) = periods(ims + 1)
            end do
        end if
    end if
    
    ! Load the state    
    ! Initial time
    if (ms > 1) then
        write(file_ext, "(i6.6)") 0
        fname = 'state.'//file_ext//'-0'
        inquire(file=fname, exist=fexist)
        if (fexist) then
            ! read from disk
            do ims = 0, ms - 1
                write(ims_str, "(i1.1)") ims
                fname = 'state.'//file_ext//'-'//ims_str
                call fieldio_read(vel_vfieldk_now)
                call solver_vectorize(vel_vfieldk_now, ims, new_x)
            end do
        else
            ! construct from time zero
            write(file_ext, "(i6.6)") 0
            fname = 'state.'//file_ext
            call fieldio_read(vel_vfieldk_now)
            call solver_vectorize(vel_vfieldk_now, 0, new_x)
            scaleT = 1.0_dp
            do ims = 1, ms - 1
                call solver_steporbit(ndtss(ims), ims-1, new_x)
                call solver_vectorize(vel_vfieldk_now, ims, new_x)
            end do
        end if
    else
        write(file_ext, "(i6.6)") 0
        fname = 'state.'//file_ext
        call fieldio_read(vel_vfieldk_now)
        call solver_vectorize(vel_vfieldk_now, 0, new_x)
    end if
    time = 0
    itime = 0

    ! Set the scales
    d = solver_dotprod(0,new_x,new_x)
    if(abs(d) < small) d=1.0_dp
    scaleT = periods(1) / sqrt(d)
    if (my_id == 0 .and. nscalars > 0) then
        if (find_period) then
            do ims  = 0, ms - 1
                new_x(ims*nnewt_pershot+1) = new_x(ims*nnewt_pershot+1) / scaleT
            end do
        end if
    end if

    tol  = rel_err * sqrt(d)
    del  = del     * sqrt(d)
    mndl = mndl    * sqrt(d)
    mxdl = mxdl    * sqrt(d)
    
    info = 1
    call newtonhook(getrhs, multJ, saveorbit, solver_dotprod_ms, &
                    mgmres, nnewt, gtol, tol, del, mndl, mxdl, nits, info, out)

    if (info == 0) then
        call newton_signal_converged
    else
        call newton_signal_not_converged
    end if
     
    call run_exit   
    
end program newton
