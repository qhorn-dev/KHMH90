module mod_grid
  use mod_kinds, only: i4, rk
  implicit none
  private

  type, public :: grid_info_t
    integer(i4) :: nx = 0
    integer(i4) :: ny = 0
    integer(i4) :: nt = 0
    real(rk) :: dx = 0.0_rk
    real(rk) :: dy = 0.0_rk
    real(rk) :: dz = 0.0_rk
    real(rk), allocatable :: x(:), y(:), z(:)
  end type grid_info_t

  public :: build_r_values, valid_i_bounds

contains

  subroutine build_r_values(r_min_pts, r_max_pts, r_step_pts, r_pts)
    integer(i4), intent(in) :: r_min_pts, r_max_pts, r_step_pts
    integer(i4), allocatable, intent(out) :: r_pts(:)
    integer(i4) :: nr, n

    if (r_step_pts <= 0) stop 'r_step_pts must be > 0'
    if (r_max_pts < r_min_pts) stop 'r_max_pts must be >= r_min_pts'

    n = (r_max_pts - r_min_pts) / r_step_pts + 1
    allocate(r_pts(n))
    do nr = 1, n
      r_pts(nr) = r_min_pts + (nr - 1) * r_step_pts
    end do
  end subroutine build_r_values

  subroutine valid_i_bounds(n, r, i_start, i_end)
    integer(i4), intent(in) :: n, r
    integer(i4), intent(out) :: i_start, i_end

    i_start = 1
    i_end = n - r
    if (i_end < i_start) then
      i_start = 1
      i_end = 0
    end if
  end subroutine valid_i_bounds

end module mod_grid
