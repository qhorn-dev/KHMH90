module mod_interscale_transfer_terms
  use mod_kinds, only: i4, i8, rk
  use mod_grid, only: grid_info_t
  implicit none
  private

  type, public :: interscale_transfer_accumulator_t
    real(rk), allocatable :: divrx_rx(:), divry_rx(:), divr_rx(:)
    real(rk), allocatable :: divrx_ry(:), divry_ry(:), divr_ry(:)
    integer(i8), allocatable :: count_rx(:), count_ry(:)
  end type interscale_transfer_accumulator_t

  type, public :: interscale_transfer_outputs_t
    real(rk), allocatable :: divrx_rx(:), divry_rx(:), divr_rx(:)
    real(rk), allocatable :: divrx_ry(:), divry_ry(:), divr_ry(:)
  end type interscale_transfer_outputs_t

  public :: init_interscale_transfer_accumulator
  public :: accumulate_interscale_transfer_block
  public :: finalize_interscale_transfer

contains

  subroutine central_diff_3d_x(f, dx, df)
    real(rk), intent(in) :: f(:, :, :), dx
    real(rk), allocatable, intent(out) :: df(:, :, :)
    integer :: nt, ny, nx
    nt = size(f, 1)
    ny = size(f, 2)
    nx = size(f, 3)
    allocate(df(nt, ny - 2, nx - 2))
    df = (f(:, 2:ny-1, 3:nx) - f(:, 2:ny-1, 1:nx-2)) / (2.0_rk * dx)
  end subroutine central_diff_3d_x

  subroutine central_diff_3d_y(f, dy, df)
    real(rk), intent(in) :: f(:, :, :), dy
    real(rk), allocatable, intent(out) :: df(:, :, :)
    integer :: nt, ny, nx
    nt = size(f, 1)
    ny = size(f, 2)
    nx = size(f, 3)
    allocate(df(nt, ny - 2, nx - 2))
    df = (f(:, 3:ny, 2:nx-1) - f(:, 1:ny-2, 2:nx-1)) / (2.0_rk * dy)
  end subroutine central_diff_3d_y

  subroutine init_interscale_transfer_accumulator(acc, nrx, nry)
    type(interscale_transfer_accumulator_t), intent(inout) :: acc
    integer(i4), intent(in) :: nrx, nry

    allocate(acc%divrx_rx(nrx), acc%divry_rx(nrx), acc%divr_rx(nrx))
    allocate(acc%divrx_ry(nry), acc%divry_ry(nry), acc%divr_ry(nry))
    allocate(acc%count_rx(nrx), acc%count_ry(nry))

    acc%divrx_rx = 0.0_rk
    acc%divry_rx = 0.0_rk
    acc%divr_rx = 0.0_rk
    acc%divrx_ry = 0.0_rk
    acc%divry_ry = 0.0_rk
    acc%divr_ry = 0.0_rk
    acc%count_rx = 0_i8
    acc%count_ry = 0_i8
  end subroutine init_interscale_transfer_accumulator

  subroutine accumulate_interscale_transfer_block(u, v, w, grid, rx_pts, ry_pts, acc)
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    type(grid_info_t), intent(in) :: grid
    integer(i4), intent(in) :: rx_pts(:), ry_pts(:)
    type(interscale_transfer_accumulator_t), intent(inout) :: acc

    real(rk), allocatable :: uc(:, :, :), vc(:, :, :), wc(:, :, :)
    real(rk), allocatable :: du_dx(:, :, :), du_dy(:, :, :), dv_dx(:, :, :), dv_dy(:, :, :), dw_dx(:, :, :), dw_dy(:, :, :)
    integer(i4) :: ir, r
    integer(i8) :: c

    allocate(uc(size(u,1), grid%ny-2, grid%nx-2))
    allocate(vc(size(v,1), grid%ny-2, grid%nx-2))
    allocate(wc(size(w,1), grid%ny-2, grid%nx-2))
    uc = u(:, 2:grid%ny-1, 2:grid%nx-1)
    vc = v(:, 2:grid%ny-1, 2:grid%nx-1)
    wc = w(:, 2:grid%ny-1, 2:grid%nx-1)

    call central_diff_3d_x(u, grid%dx, du_dx)
    call central_diff_3d_y(u, grid%dy, du_dy)
    call central_diff_3d_x(v, grid%dx, dv_dx)
    call central_diff_3d_y(v, grid%dy, dv_dy)
    call central_diff_3d_x(w, grid%dx, dw_dx)
    call central_diff_3d_y(w, grid%dy, dw_dy)

    !$omp parallel do default(none) shared(rx_pts,uc,vc,wc,du_dx,du_dy,dv_dx,dv_dy,dw_dx,dw_dy,acc) &
    !$omp private(ir,r,c) schedule(static)
    do ir = 1, size(rx_pts)
      r = rx_pts(ir)
      call accumulate_one_rx(ir, r, uc, vc, wc, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, acc, c)
      acc%count_rx(ir) = acc%count_rx(ir) + c
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(ry_pts,uc,vc,wc,du_dx,du_dy,dv_dx,dv_dy,dw_dx,dw_dy,acc) &
    !$omp private(ir,r,c) schedule(static)
    do ir = 1, size(ry_pts)
      r = ry_pts(ir)
      call accumulate_one_ry(ir, r, uc, vc, wc, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, acc, c)
      acc%count_ry(ir) = acc%count_ry(ir) + c
    end do
    !$omp end parallel do
  end subroutine accumulate_interscale_transfer_block

  subroutine accumulate_one_rx(ir, r, u, v, w, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(rk), intent(in) :: du_dx(:, :, :), du_dy(:, :, :), dv_dx(:, :, :), dv_dy(:, :, :), dw_dx(:, :, :), dw_dy(:, :, :)
    type(interscale_transfer_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: ddx(:, :, :), ddy(:, :, :), vdx(:, :, :), vdy(:, :, :), wdx(:, :, :), wdy(:, :, :)
    real(rk) :: s1, s2

    allocate(du(size(u,1), size(u,2), size(u,3)-r))
    allocate(dv(size(v,1), size(v,2), size(v,3)-r))
    allocate(dw(size(w,1), size(w,2), size(w,3)-r))
    allocate(ddx(size(du_dx,1), size(du_dx,2), size(du_dx,3)-r))
    allocate(ddy(size(du_dy,1), size(du_dy,2), size(du_dy,3)-r))
    allocate(vdx(size(dv_dx,1), size(dv_dx,2), size(dv_dx,3)-r))
    allocate(vdy(size(dv_dy,1), size(dv_dy,2), size(dv_dy,3)-r))
    allocate(wdx(size(dw_dx,1), size(dw_dx,2), size(dw_dx,3)-r))
    allocate(wdy(size(dw_dy,1), size(dw_dy,2), size(dw_dy,3)-r))

    du = u(:, :, 1+r:size(u,3)) - u(:, :, 1:size(u,3)-r)
    dv = v(:, :, 1+r:size(v,3)) - v(:, :, 1:size(v,3)-r)
    dw = w(:, :, 1+r:size(w,3)) - w(:, :, 1:size(w,3)-r)

    ddx = 0.5_rk * (du_dx(:, :, 1+r:size(du_dx,3)) + du_dx(:, :, 1:size(du_dx,3)-r))
    ddy = 0.5_rk * (du_dy(:, :, 1+r:size(du_dy,3)) + du_dy(:, :, 1:size(du_dy,3)-r))
    vdx = 0.5_rk * (dv_dx(:, :, 1+r:size(dv_dx,3)) + dv_dx(:, :, 1:size(dv_dx,3)-r))
    vdy = 0.5_rk * (dv_dy(:, :, 1+r:size(dv_dy,3)) + dv_dy(:, :, 1:size(dv_dy,3)-r))
    wdx = 0.5_rk * (dw_dx(:, :, 1+r:size(dw_dx,3)) + dw_dx(:, :, 1:size(dw_dx,3)-r))
    wdy = 0.5_rk * (dw_dy(:, :, 1+r:size(dw_dy,3)) + dw_dy(:, :, 1:size(dw_dy,3)-r))

    s1 = 2.0_rk * sum(du * du * ddx + du * dv * vdx + du * dw * wdx)
    s2 = 2.0_rk * sum(dv * du * ddy + dv * dv * vdy + dv * dw * wdy)

    acc%divrx_rx(ir) = acc%divrx_rx(ir) + s1
    acc%divry_rx(ir) = acc%divry_rx(ir) + s2
    acc%divr_rx(ir) = acc%divr_rx(ir) + s1 + s2

    c = int(size(du,1), i8) * int(size(du,2), i8) * int(size(du,3), i8)
  end subroutine accumulate_one_rx

  subroutine accumulate_one_ry(ir, r, u, v, w, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(rk), intent(in) :: du_dx(:, :, :), du_dy(:, :, :), dv_dx(:, :, :), dv_dy(:, :, :), dw_dx(:, :, :), dw_dy(:, :, :)
    type(interscale_transfer_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: ddx(:, :, :), ddy(:, :, :), vdx(:, :, :), vdy(:, :, :), wdx(:, :, :), wdy(:, :, :)
    real(rk) :: s1, s2

    allocate(du(size(u,1), size(u,2)-r, size(u,3)))
    allocate(dv(size(v,1), size(v,2)-r, size(v,3)))
    allocate(dw(size(w,1), size(w,2)-r, size(w,3)))
    allocate(ddx(size(du_dx,1), size(du_dx,2)-r, size(du_dx,3)))
    allocate(ddy(size(du_dy,1), size(du_dy,2)-r, size(du_dy,3)))
    allocate(vdx(size(dv_dx,1), size(dv_dx,2)-r, size(dv_dx,3)))
    allocate(vdy(size(dv_dy,1), size(dv_dy,2)-r, size(dv_dy,3)))
    allocate(wdx(size(dw_dx,1), size(dw_dx,2)-r, size(dw_dx,3)))
    allocate(wdy(size(dw_dy,1), size(dw_dy,2)-r, size(dw_dy,3)))

    du = u(:, 1+r:size(u,2), :) - u(:, 1:size(u,2)-r, :)
    dv = v(:, 1+r:size(v,2), :) - v(:, 1:size(v,2)-r, :)
    dw = w(:, 1+r:size(w,2), :) - w(:, 1:size(w,2)-r, :)

    ddx = 0.5_rk * (du_dx(:, 1+r:size(du_dx,2), :) + du_dx(:, 1:size(du_dx,2)-r, :))
    ddy = 0.5_rk * (du_dy(:, 1+r:size(du_dy,2), :) + du_dy(:, 1:size(du_dy,2)-r, :))
    vdx = 0.5_rk * (dv_dx(:, 1+r:size(dv_dx,2), :) + dv_dx(:, 1:size(dv_dx,2)-r, :))
    vdy = 0.5_rk * (dv_dy(:, 1+r:size(dv_dy,2), :) + dv_dy(:, 1:size(dv_dy,2)-r, :))
    wdx = 0.5_rk * (dw_dx(:, 1+r:size(dw_dx,2), :) + dw_dx(:, 1:size(dw_dx,2)-r, :))
    wdy = 0.5_rk * (dw_dy(:, 1+r:size(dw_dy,2), :) + dw_dy(:, 1:size(dw_dy,2)-r, :))

    s1 = 2.0_rk * sum(du * du * ddx + du * dv * vdx + du * dw * wdx)
    s2 = 2.0_rk * sum(dv * du * ddy + dv * dv * vdy + dv * dw * wdy)

    acc%divrx_ry(ir) = acc%divrx_ry(ir) + s1
    acc%divry_ry(ir) = acc%divry_ry(ir) + s2
    acc%divr_ry(ir) = acc%divr_ry(ir) + s1 + s2

    c = int(size(du,1), i8) * int(size(du,2), i8) * int(size(du,3), i8)
  end subroutine accumulate_one_ry

  subroutine finalize_interscale_transfer(acc, out)
    type(interscale_transfer_accumulator_t), intent(in) :: acc
    type(interscale_transfer_outputs_t), intent(out) :: out

    allocate(out%divrx_rx(size(acc%divrx_rx)))
    allocate(out%divry_rx(size(acc%divry_rx)))
    allocate(out%divr_rx(size(acc%divr_rx)))
    allocate(out%divrx_ry(size(acc%divrx_ry)))
    allocate(out%divry_ry(size(acc%divry_ry)))
    allocate(out%divr_ry(size(acc%divr_ry)))

    out%divrx_rx = acc%divrx_rx / real(acc%count_rx, rk)
    out%divry_rx = acc%divry_rx / real(acc%count_rx, rk)
    out%divr_rx = acc%divr_rx / real(acc%count_rx, rk)
    out%divrx_ry = acc%divrx_ry / real(acc%count_ry, rk)
    out%divry_ry = acc%divry_ry / real(acc%count_ry, rk)
    out%divr_ry = acc%divr_ry / real(acc%count_ry, rk)
  end subroutine finalize_interscale_transfer

end module mod_interscale_transfer_terms
