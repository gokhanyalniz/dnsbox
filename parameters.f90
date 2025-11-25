module parameters
    use numbers
    use openmpi
    use io

    !# Geometry & discretization
    integer(i4) :: &
        nx = 48, ny = 48, nz = 48, & ! number of grid points
        nxx, nyy, nzz, & ! number of super-sampled grid points
        nx_half, ny_half, nz_half, & ! shorthands
        nyy_half_pad1, & ! nyy/2+1 are frequently used in r2c 
                                 ! transforms see "2.3 One-Dimensional DFTs of
                                 ! Real Data" / fftw3 manual
        nx_perproc, nzz_perproc, &
        nz_0, nz_perproc, ny_half_pad1

    integer(i4), parameter :: supsamp_fac = 3 ! supsamp_fac / 2 dealiasing

    real(dp)            :: Lx = 4.0_dp, Lz = 4.0_dp ! grid dimensions
    real(dp), parameter :: Ly = 4.0_dp ! fixed by non-dimensionalization

    !# Physics
    integer(i4) :: forcing = 1 ! none = 0, sine = 1, cosine = 2
    integer(i4), parameter :: qF = 1 ! Forcing wave number
    real(dp) :: Re = 630.0_dp

    !# Initiation
    integer(i4) :: IC = -1, &     ! Initial condition 
                               ! (-3 shapiro, -2 laminar, -1 random, 
                               !   0 state.{istart / i_save_fields)
                   i_start = 0 ! Starting time step  
    integer(i8) :: random_seed = -1 ! -1 reads from /dev/urandom
    real(dp) :: random_energy = 0.1_dp, & ! (ekin_laminar) energy of the random IC
                random_smooth = 0.9_dp, & ! smoothness factor for random IC
                t_start = 0 ! Starting time
    
    !# Results -- note i_* convention for variables in time step units
    integer(i4) :: i_print_stats = 20, &   ! Print energy, dissipation etc.
                   i_print_steps = 20, &   ! Print Courant, step error etc.
                   i_save_fields = 2000, & ! Save fields (restart files)
                   i_flush = -1, &
                   i_save_phys = -1

    !# Time stepping 
    real(dp)    :: dt = 0.025_dp, &
                   implicitness = 0.5_dp, &
                   steptol = 1.0e-9_dp, &
                   dtmax = 0.1_dp, &
                   courant_target = 0.25_dp
    integer(i4) :: ncorr = 10 
    logical :: adaptive_dt = .true.
    
    !# Termination
    logical  :: terminate_laminar = .true. ! terminate on laminarization
    real(dp) :: relerr_lam = 0.1_dp ! conclude laminarization if stats within
                                       ! relerr_lam of those of laminar
    real(sp) :: wall_clock_limit = -1.0_sp
    integer(i4) :: i_finish = -1 ! time step limit

    integer(i4) :: itime ! simulation timestep count
    real(dp)    :: time, &  ! simulation time
                   dx, dy, dz, &! grid spacings in the expanded physical space
                   norm_fft, &! fft forward+backward needs to be normalized
                   kF ! forcing wave number

    ! laminar values
    real(dp) :: ekin_lam, powerin_lam, dissip_lam

    ! forcing coefficient
    real(dp) :: amp

    ! Given 3x3 symmetric matrix M, entries M_{ij} will be used
    integer(i4), parameter, dimension(6) :: isym = (/1, 1, 1, 2, 2, 3/), &
                                            jsym = (/1, 2, 3, 2, 3, 3/)
    integer(i4) :: nsym(3,3)

    namelist /grid/ nx, ny, nz, Lx, Lz
    namelist /physics/ forcing, Re
    namelist /initiation/ IC, random_seed, random_energy, random_smooth, &
                          t_start, i_start
    namelist /output/ i_print_stats, i_print_steps, i_save_fields, &
                      i_flush, i_save_phys
    namelist /time_stepping/ dt, implicitness, steptol, ncorr, adaptive_dt, &
                             dtmax, courant_target
    namelist /termination/ terminate_laminar, relerr_lam, wall_clock_limit, i_finish

    contains 

!==============================================================================

    subroutine parameters_init

        integer(i4) :: in, i, j, n

        do n = 1,6
            i = isym(n)
            j = jsym(n)
            nsym(i,j) = n
            nsym(j,i) = n
        end do

        ! Is parameters.in here?
        inquire(file = 'parameters.in', exist=there)
        if(.not.there) then
            write(out,*) 'cannot find parameters.in'
            flush(out)
            error stop
        end if
        
        ! Read parameters from the input file
        open(newunit=in, file='parameters.in', status="old")
        read(in, nml=grid)
        read(in, nml=physics)
        read(in, nml=symmetries)
        read(in, nml=initiation)
        read(in, nml=output)
        read(in, nml=time_stepping)
        read(in, nml=termination)
        close(in)

        write(out, *) 'dnsbox revision: ', revision
        write(out, *) 'num_procs: ', num_procs
        write(out, *) 'forcing = ', forcing
        write(out, *) 'epsilon = ', epsilon
        write(out, *) 'small = ', small

        write(out, '(79(''=''))')

        ! ---------------------------------------------------------------------

        write(out, *) 'Re = ', Re

        ! ---------------------------------------------------------------------

        write(out, "('nx, ny, nz', 3i4)") nx, ny, nz
        write(out, *) 'Lx = ', Lx
        write(out, *) 'Ly = ', Ly
        write(out, *) 'Lz = ', Lz

        if (mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0 .or. mod(nz, 2) /= 0) then
            write(out, *) 'nx, ny and nz must be even'
            flush(out)
            error stop
        end if
        
        nx_perproc = nx / num_procs
        nx_half = nx / 2 
        ny_half = ny / 2 
        nz_half = nz / 2 

        if (nx_perproc * num_procs /= nx) then
            write(out, *) '*** wrong nx:', nx, & 
                          '*** should be divisible by num_procs:', num_procs
            flush(out)
            error stop
        end if

        ! compute the size of the supersampled grid
        nxx = supsamp_fac * nx_half
        nyy = supsamp_fac * ny_half
        nzz = supsamp_fac * nz_half

        nyy_half_pad1 = nyy / 2 + 1
        ny_half_pad1 = ny / 2 + 1

        nzz_perproc  = nzz / num_procs
        nz_0 = nz
        nz_perproc  = nz_0 / num_procs

        ! Pad z in the physical dimension if necessary to make it divisible
        ! by num_procs
        if (nzz_perproc * num_procs /= nzz) then
            nzz_perproc = nzz_perproc + 1
            nzz = nzz_perproc * num_procs

        end if

        if (nz_perproc * num_procs /= nz_0) then
            nz_perproc = nz_perproc + 1
            nz_0 = nz_perproc * num_procs

        end if

        ! Not sure what happens to the transpose algorithm if in either
        ! wavenumber or physical space we end up with just one mode in the
        ! distributed dimension
        if (nx_perproc <= 1 .or. nzz_perproc <= 1) then
            write(out, *) '*** num_procs too large:', num_procs, & 
                          '*** nx_perproc:', nx_perproc, &
                          '*** nzz_perproc:', nzz_perproc
            flush(out)
            error stop
        end if

        if (i_save_phys > 0 .and. nz_perproc <= 1) then
            write(out, *) '*** num_procs too large:', num_procs, & 
                          '*** nz_perproc:', nz_perproc
            flush(out)
            error stop
        end if

        write(out, '(79(''=''))')

        ! ---------------------------------------------------------------------

        write(out, "('nxx, nyy, nzz', 3i4)") nxx, nyy, nzz

        norm_fft = 1.0_dp / (nxx * nyy * nzz)
        
        dx = Lx / nxx
        dy = Ly / nyy
        dz = Lz / nzz
        
        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        amp = PI**2

        ! Set the forcing wavenumber
        kF = 2.0_dp * PI * qF / Ly

        ! compute laminar values
        ekin_lam    = 1.0_dp / 4.0_dp
        powerin_lam = amp / (8 * Re)
        dissip_lam  = powerin_lam
        
        write(out, '(79(''=''))')

        ! ---------------------------------------------------------------------  

        write(out, *) 'IC = ', IC
        write(out, *) 'random_seed = ', random_seed
        write(out, *) 'random_energy = ', random_energy
        write(out, *) 'random_smooth = ', random_smooth
        
        write(out, '(79(''=''))')
                       
        ! ---------------------------------------------------------------------

        write(out, *) 't_start = ', t_start
        write(out, *) 'i_start = ', i_start
        write(out, *) 'i_print_stats = ', i_print_stats
        write(out, *) 'i_print_steps = ', i_print_steps
        write(out, *) 'i_save_fields = ', i_save_fields
        write(out, *) 'i_flush = ', i_flush
        write(out, *) 'i_save_phys = ', i_save_phys

        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        write(out, *) 'dt = ', dt
        write(out, *) 'implicitness = ', implicitness
        write(out, *) 'steptol = ', steptol
        write(out, *) 'ncorr = ', ncorr
        write(out, *) 'adaptive_dt = ', adaptive_dt
        write(out, *) 'dtmax = ', dtmax
        write(out, *) 'courant_target = ', courant_target

        write(out, '(79(''=''))')
        
        ! ---------------------------------------------------------------------

        write(out, *) 'terminate_laminar = ', terminate_laminar
        write(out, *) 'relerr_lam = ', relerr_lam
        write(out, *) 'wall_clock_limit = ', wall_clock_limit
        write(out, *) 'i_finish = ', i_finish

        write(out, '(79(''=''))')

    end subroutine parameters_init

!==============================================================================
    
end module parameters
