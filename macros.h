#ifndef macros_h
#define macros_h
#define _indices integer(i4) :: ix, iy, iz
#define _indicess integer(i4) :: ixx, iyy, izz
#define _indices0 integer(i4) :: ix0, iy0, iz0

#define _loop_spec_begin /*
*/ do iz = 1, nz; if(iz == iz_max) cycle; do iy = 1, ny_half; do ix = 1, nx_perproc; if(ix_max /= -1 .and. ix == ix_max) cycle;
         
#define _loop_spec_end end do; end do; end do

#define _loop_phys_begin /*
*/ do ixx = 1, nxx; do izz = 1, nzz_perproc; do iyy = 1, nyy;
         
#define _loop_phys_end end do; end do; end do

#define _loop_phys_0_begin /*
*/ do ix0 = 1, nx; do iz0 = 1, nz_perproc; do iy0 = 1, ny;
         
#define _loop_phys_0_end end do; end do; end do
#endif
