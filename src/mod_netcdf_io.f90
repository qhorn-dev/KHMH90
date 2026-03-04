module mod_netcdf_io
  use mod_kinds, only: i4, rk
  use mod_grid, only: grid_info_t
  use netcdf
  implicit none
  private

  type, public :: run_reader_t
    character(len=256) :: path = ''
    integer(i4) :: nt = 0
    integer(i4) :: ny = 0
    integer(i4) :: nx = 0
    integer :: ncid = -1
    integer :: varid_u = -1
    integer :: varid_v = -1
    integer :: varid_w = -1
  end type run_reader_t

  public :: load_stats_and_grid
  public :: open_run_reader, close_run_reader
  public :: read_run_block

contains

  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where
    if (status /= nf90_noerr) then
      write(*,'(a)') 'NetCDF error in '//trim(where)//': '//trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

  subroutine load_stats_and_grid(stats_file, grid, mean_u, mean_v, mean_w)
    character(len=*), intent(in) :: stats_file
    type(grid_info_t), intent(out) :: grid
    real(rk), allocatable, intent(out) :: mean_u(:, :), mean_v(:, :), mean_w(:, :)

    integer :: ncid, varid, dimid_y, dimid_x
    integer :: ny, nx
    real(rk), allocatable :: gx(:, :), gy(:, :), mu(:, :), mv(:, :), mw(:, :)
    real(rk), allocatable :: gx_c(:, :), gy_c(:, :), mu_c(:, :), mv_c(:, :), mw_c(:, :)
    integer :: start2(2), count2(2)

    call nc_check(nf90_open(trim(stats_file), nf90_nowrite, ncid), 'open stats_file')
    call nc_check(nf90_inq_dimid(ncid, 'dim_y', dimid_y), 'inq dim_y')
    call nc_check(nf90_inquire_dimension(ncid, dimid_y, len=ny), 'len dim_y')
    call nc_check(nf90_inq_dimid(ncid, 'dim_x', dimid_x), 'inq dim_x')
    call nc_check(nf90_inquire_dimension(ncid, dimid_x, len=nx), 'len dim_x')

    allocate(gx(ny, nx), gy(ny, nx), mu(ny, nx), mv(ny, nx), mw(ny, nx))
    allocate(gx_c(nx, ny), gy_c(nx, ny), mu_c(nx, ny), mv_c(nx, ny), mw_c(nx, ny))
    allocate(mean_u(ny, nx), mean_v(ny, nx), mean_w(ny, nx))
    allocate(grid%x(nx), grid%y(ny))
    start2 = (/ 1, 1 /)
    count2 = (/ nx, ny /)

    call nc_check(nf90_inq_varid(ncid, 'grid_x', varid), 'varid grid_x')
    call nc_check(nf90_get_var(ncid, varid, gx_c, start=start2, count=count2), 'read grid_x')
    call nc_check(nf90_inq_varid(ncid, 'grid_y', varid), 'varid grid_y')
    call nc_check(nf90_get_var(ncid, varid, gy_c, start=start2, count=count2), 'read grid_y')
    call nc_check(nf90_inq_varid(ncid, 'vel_x_mean', varid), 'varid vel_x_mean')
    call nc_check(nf90_get_var(ncid, varid, mu_c, start=start2, count=count2), 'read vel_x_mean')
    call nc_check(nf90_inq_varid(ncid, 'vel_y_mean', varid), 'varid vel_y_mean')
    call nc_check(nf90_get_var(ncid, varid, mv_c, start=start2, count=count2), 'read vel_y_mean')
    call nc_check(nf90_inq_varid(ncid, 'vel_z_mean', varid), 'varid vel_z_mean')
    call nc_check(nf90_get_var(ncid, varid, mw_c, start=start2, count=count2), 'read vel_z_mean')
    call nc_check(nf90_close(ncid), 'close stats_file')

    gx = transpose(gx_c)
    gy = transpose(gy_c)
    mu = transpose(mu_c)
    mv = transpose(mv_c)
    mw = transpose(mw_c)

    grid%ny = ny
    grid%nx = nx
    grid%nt = 0
    grid%x = gx(1, :)
    grid%y = gy(:, 1)
    if (nx > 1) grid%dx = grid%x(2) - grid%x(1)
    if (ny > 1) grid%dy = grid%y(2) - grid%y(1)
    grid%dz = 0.0_rk

    mean_u = mu
    mean_v = mv
    mean_w = mw

    deallocate(gx, gy, mu, mv, mw)
    deallocate(gx_c, gy_c, mu_c, mv_c, mw_c)
  end subroutine load_stats_and_grid

  subroutine open_run_reader(path, rr)
    character(len=*), intent(in) :: path
    type(run_reader_t), intent(out) :: rr
    integer :: dimid_t, dimid_y, dimid_x

    rr%path = path
    call nc_check(nf90_open(trim(path), nf90_nowrite, rr%ncid), 'open run file')
    call nc_check(nf90_inq_dimid(rr%ncid, 'dim_t', dimid_t), 'inq dim_t')
    call nc_check(nf90_inquire_dimension(rr%ncid, dimid_t, len=rr%nt), 'len dim_t')
    call nc_check(nf90_inq_dimid(rr%ncid, 'dim_y', dimid_y), 'inq dim_y')
    call nc_check(nf90_inquire_dimension(rr%ncid, dimid_y, len=rr%ny), 'len dim_y')
    call nc_check(nf90_inq_dimid(rr%ncid, 'dim_x', dimid_x), 'inq dim_x')
    call nc_check(nf90_inquire_dimension(rr%ncid, dimid_x, len=rr%nx), 'len dim_x')
    call nc_check(nf90_inq_varid(rr%ncid, 'vel_x', rr%varid_u), 'varid vel_x')
    call nc_check(nf90_inq_varid(rr%ncid, 'vel_y', rr%varid_v), 'varid vel_y')
    call nc_check(nf90_inq_varid(rr%ncid, 'vel_z', rr%varid_w), 'varid vel_z')
  end subroutine open_run_reader

  subroutine close_run_reader(rr)
    type(run_reader_t), intent(inout) :: rr
    if (rr%ncid >= 0) call nc_check(nf90_close(rr%ncid), 'close run file')
    rr%path = ''
    rr%nt = 0
    rr%ny = 0
    rr%nx = 0
    rr%ncid = -1
    rr%varid_u = -1
    rr%varid_v = -1
    rr%varid_w = -1
  end subroutine close_run_reader

  subroutine read_run_block(rr, t_start, t_count, u, v, w)
    type(run_reader_t), intent(in) :: rr
    integer(i4), intent(in) :: t_start, t_count
    real(rk), allocatable, intent(out) :: u(:, :, :), v(:, :, :), w(:, :, :)

    integer :: start(3), count(3), t
    real(rk), allocatable :: ux(:, :, :), vx(:, :, :), wx(:, :, :)

    if (t_count <= 0) stop 'read_run_block: invalid t_count'
    if (t_start < 1 .or. t_start + t_count - 1 > rr%nt) stop 'read_run_block: time range out of bounds'

    allocate(u(t_count, rr%ny, rr%nx), v(t_count, rr%ny, rr%nx), w(t_count, rr%ny, rr%nx))
    allocate(ux(rr%nx, rr%ny, t_count), vx(rr%nx, rr%ny, t_count), wx(rr%nx, rr%ny, t_count))

    start = (/ 1, 1, t_start /)
    count = (/ rr%nx, rr%ny, t_count /)
    call nc_check(nf90_get_var(rr%ncid, rr%varid_u, ux, start=start, count=count), 'read vel_x block')
    call nc_check(nf90_get_var(rr%ncid, rr%varid_v, vx, start=start, count=count), 'read vel_y block')
    call nc_check(nf90_get_var(rr%ncid, rr%varid_w, wx, start=start, count=count), 'read vel_z block')

    do t = 1, t_count
      u(t, :, :) = transpose(ux(:, :, t))
      v(t, :, :) = transpose(vx(:, :, t))
      w(t, :, :) = transpose(wx(:, :, t))
    end do

    deallocate(ux, vx, wx)
  end subroutine read_run_block

end module mod_netcdf_io
