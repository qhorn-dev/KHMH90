module mod_linear_transport_terms
  use mod_kinds, only: i4, i8, rk
  use mod_grid, only: grid_info_t
  implicit none
  private

  type, public :: linear_transport_context_t
    real(rk), allocatable :: mean_u(:, :), mean_v(:, :), mean_w(:, :)
    real(rk), allocatable :: dmean_u_dx(:, :), dmean_u_dy(:, :)
    real(rk), allocatable :: dmean_v_dx(:, :), dmean_v_dy(:, :)
    real(rk), allocatable :: dmean_w_dx(:, :), dmean_w_dy(:, :)
  end type linear_transport_context_t

  type, public :: linear_transport_accumulator_t
    real(rk), allocatable :: ux_dS2_dXx_rx(:), uy_dS2_dXy_rx(:)
    real(rk), allocatable :: dux_dS2_drx_rx(:), duy_dS2_dry_rx(:)
    real(rk), allocatable :: ux_dS2_dXx_ry(:), uy_dS2_dXy_ry(:)
    real(rk), allocatable :: dux_dS2_drx_ry(:), duy_dS2_dry_ry(:)
    integer(i8), allocatable :: count_rx(:), count_ry(:)
  end type linear_transport_accumulator_t

  type, public :: linear_transport_outputs_t
    real(rk), allocatable :: ux_dS2_dXx_rx(:), uy_dS2_dXy_rx(:)
    real(rk), allocatable :: dux_dS2_drx_rx(:), duy_dS2_dry_rx(:)
    real(rk), allocatable :: ux_dS2_dXx_ry(:), uy_dS2_dXy_ry(:)
    real(rk), allocatable :: dux_dS2_drx_ry(:), duy_dS2_dry_ry(:)
  end type linear_transport_outputs_t

  public :: init_linear_transport_context
  public :: init_linear_transport_accumulator
  public :: accumulate_linear_transport_block
  public :: finalize_linear_transport

contains

  subroutine central_diff_2d_x(f, dx, df)
    real(rk), intent(in) :: f(:, :), dx
    real(rk), allocatable, intent(out) :: df(:, :)
    integer :: ny, nx
    ny = size(f, 1)
    nx = size(f, 2)
    allocate(df(ny - 2, nx - 2))
    df = (f(2:ny-1, 3:nx) - f(2:ny-1, 1:nx-2)) / (2.0_rk * dx)
  end subroutine central_diff_2d_x

  subroutine central_diff_2d_y(f, dy, df)
    real(rk), intent(in) :: f(:, :), dy
    real(rk), allocatable, intent(out) :: df(:, :)
    integer :: ny, nx
    ny = size(f, 1)
    nx = size(f, 2)
    allocate(df(ny - 2, nx - 2))
    df = (f(3:ny, 2:nx-1) - f(1:ny-2, 2:nx-1)) / (2.0_rk * dy)
  end subroutine central_diff_2d_y

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

  subroutine init_linear_transport_context(ctx, grid, mean_u_full, mean_v_full, mean_w_full)
    type(linear_transport_context_t), intent(inout) :: ctx
    type(grid_info_t), intent(in) :: grid
    real(rk), intent(in) :: mean_u_full(:, :), mean_v_full(:, :), mean_w_full(:, :)

    allocate(ctx%mean_u(grid%ny - 2, grid%nx - 2))
    allocate(ctx%mean_v(grid%ny - 2, grid%nx - 2))
    allocate(ctx%mean_w(grid%ny - 2, grid%nx - 2))

    ctx%mean_u = mean_u_full(2:grid%ny-1, 2:grid%nx-1)
    ctx%mean_v = mean_v_full(2:grid%ny-1, 2:grid%nx-1)
    ctx%mean_w = mean_w_full(2:grid%ny-1, 2:grid%nx-1)

    call central_diff_2d_x(mean_u_full, grid%dx, ctx%dmean_u_dx)
    call central_diff_2d_y(mean_u_full, grid%dy, ctx%dmean_u_dy)
    call central_diff_2d_x(mean_v_full, grid%dx, ctx%dmean_v_dx)
    call central_diff_2d_y(mean_v_full, grid%dy, ctx%dmean_v_dy)
    call central_diff_2d_x(mean_w_full, grid%dx, ctx%dmean_w_dx)
    call central_diff_2d_y(mean_w_full, grid%dy, ctx%dmean_w_dy)
  end subroutine init_linear_transport_context

  subroutine init_linear_transport_accumulator(acc, nrx, nry)
    type(linear_transport_accumulator_t), intent(inout) :: acc
    integer(i4), intent(in) :: nrx, nry

    allocate(acc%ux_dS2_dXx_rx(nrx), acc%uy_dS2_dXy_rx(nrx))
    allocate(acc%dux_dS2_drx_rx(nrx), acc%duy_dS2_dry_rx(nrx))
    allocate(acc%ux_dS2_dXx_ry(nry), acc%uy_dS2_dXy_ry(nry))
    allocate(acc%dux_dS2_drx_ry(nry), acc%duy_dS2_dry_ry(nry))
    allocate(acc%count_rx(nrx), acc%count_ry(nry))

    acc%ux_dS2_dXx_rx = 0.0_rk
    acc%uy_dS2_dXy_rx = 0.0_rk
    acc%dux_dS2_drx_rx = 0.0_rk
    acc%duy_dS2_dry_rx = 0.0_rk
    acc%ux_dS2_dXx_ry = 0.0_rk
    acc%uy_dS2_dXy_ry = 0.0_rk
    acc%dux_dS2_drx_ry = 0.0_rk
    acc%duy_dS2_dry_ry = 0.0_rk
    acc%count_rx = 0_i8
    acc%count_ry = 0_i8
  end subroutine init_linear_transport_accumulator

  subroutine accumulate_linear_transport_block(u, v, w, grid, rx_pts, ry_pts, ctx, acc)
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    type(grid_info_t), intent(in) :: grid
    integer(i4), intent(in) :: rx_pts(:), ry_pts(:)
    type(linear_transport_context_t), intent(in) :: ctx
    type(linear_transport_accumulator_t), intent(inout) :: acc

    real(rk), allocatable :: uc(:, :, :), vc(:, :, :), wc(:, :, :)
    real(rk), allocatable :: du_dx(:, :, :), du_dy(:, :, :), dv_dx(:, :, :), dv_dy(:, :, :), dw_dx(:, :, :), dw_dy(:, :, :)
    integer(i4) :: ir, r
    integer(i8) :: c

    uc = u(:, 2:grid%ny-1, 2:grid%nx-1)
    vc = v(:, 2:grid%ny-1, 2:grid%nx-1)
    wc = w(:, 2:grid%ny-1, 2:grid%nx-1)

    call central_diff_3d_x(u, grid%dx, du_dx)
    call central_diff_3d_y(u, grid%dy, du_dy)
    call central_diff_3d_x(v, grid%dx, dv_dx)
    call central_diff_3d_y(v, grid%dy, dv_dy)
    call central_diff_3d_x(w, grid%dx, dw_dx)
    call central_diff_3d_y(w, grid%dy, dw_dy)

    !$omp parallel do default(none) shared(rx_pts,uc,vc,wc,du_dx,du_dy,dv_dx,dv_dy,dw_dx,dw_dy,ctx,acc) &
    !$omp private(ir,r,c) schedule(static)
    do ir = 1, size(rx_pts)
      r = rx_pts(ir)
      call accumulate_one_rx(ir, r, uc, vc, wc, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, ctx, acc, c)
      acc%count_rx(ir) = acc%count_rx(ir) + c
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(ry_pts,uc,vc,wc,du_dx,du_dy,dv_dx,dv_dy,dw_dx,dw_dy,ctx,acc) &
    !$omp private(ir,r,c) schedule(static)
    do ir = 1, size(ry_pts)
      r = ry_pts(ir)
      call accumulate_one_ry(ir, r, uc, vc, wc, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, ctx, acc, c)
      acc%count_ry(ir) = acc%count_ry(ir) + c
    end do
    !$omp end parallel do
  end subroutine accumulate_linear_transport_block

  subroutine accumulate_one_rx(ir, r, u, v, w, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, ctx, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(rk), intent(in) :: du_dx(:, :, :), du_dy(:, :, :), dv_dx(:, :, :), dv_dy(:, :, :), dw_dx(:, :, :), dw_dy(:, :, :)
    type(linear_transport_context_t), intent(in) :: ctx
    type(linear_transport_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: mean_du(:, :), mean_dv(:, :), ux_mean(:, :), uy_mean(:, :)
    real(rk), allocatable :: ddx(:, :, :), ddy(:, :, :), vdx(:, :, :), vdy(:, :, :), wdx(:, :, :), wdy(:, :, :)
    real(rk) :: s1, s2, s3, s4

    allocate(du(size(u,1), size(u,2), size(u,3)-r))
    allocate(dv(size(v,1), size(v,2), size(v,3)-r))
    allocate(dw(size(w,1), size(w,2), size(w,3)-r))
    allocate(mean_du(size(ctx%mean_u,1), size(ctx%mean_u,2)-r))
    allocate(mean_dv(size(ctx%mean_v,1), size(ctx%mean_v,2)-r))
    allocate(ux_mean(size(ctx%mean_u,1), size(ctx%mean_u,2)-r))
    allocate(uy_mean(size(ctx%mean_v,1), size(ctx%mean_v,2)-r))
    allocate(ddx(size(du_dx,1), size(du_dx,2), size(du_dx,3)-r))
    allocate(ddy(size(du_dy,1), size(du_dy,2), size(du_dy,3)-r))
    allocate(vdx(size(dv_dx,1), size(dv_dx,2), size(dv_dx,3)-r))
    allocate(vdy(size(dv_dy,1), size(dv_dy,2), size(dv_dy,3)-r))
    allocate(wdx(size(dw_dx,1), size(dw_dx,2), size(dw_dx,3)-r))
    allocate(wdy(size(dw_dy,1), size(dw_dy,2), size(dw_dy,3)-r))

    du = u(:, :, 1+r:size(u,3)) - u(:, :, 1:size(u,3)-r)
    dv = v(:, :, 1+r:size(v,3)) - v(:, :, 1:size(v,3)-r)
    dw = w(:, :, 1+r:size(w,3)) - w(:, :, 1:size(w,3)-r)

    mean_du = ctx%mean_u(:, 1+r:size(ctx%mean_u,2)) - ctx%mean_u(:, 1:size(ctx%mean_u,2)-r)
    mean_dv = ctx%mean_v(:, 1+r:size(ctx%mean_v,2)) - ctx%mean_v(:, 1:size(ctx%mean_v,2)-r)
    ux_mean = 0.5_rk * (ctx%mean_u(:, 1+r:size(ctx%mean_u,2)) + ctx%mean_u(:, 1:size(ctx%mean_u,2)-r))
    uy_mean = 0.5_rk * (ctx%mean_v(:, 1+r:size(ctx%mean_v,2)) + ctx%mean_v(:, 1:size(ctx%mean_v,2)-r))

    ddx = du_dx(:, :, 1+r:size(du_dx,3)) + du_dx(:, :, 1:size(du_dx,3)-r)
    ddy = du_dy(:, :, 1+r:size(du_dy,3)) + du_dy(:, :, 1:size(du_dy,3)-r)
    vdx = dv_dx(:, :, 1+r:size(dv_dx,3)) + dv_dx(:, :, 1:size(dv_dx,3)-r)
    vdy = dv_dy(:, :, 1+r:size(dv_dy,3)) + dv_dy(:, :, 1:size(dv_dy,3)-r)
    wdx = dw_dx(:, :, 1+r:size(dw_dx,3)) + dw_dx(:, :, 1:size(dw_dx,3)-r)
    wdy = dw_dy(:, :, 1+r:size(dw_dy,3)) + dw_dy(:, :, 1:size(dw_dy,3)-r)

    s1 = sum(ux_mean * sum(du * (4.0_rk * 0.25_rk * (du_dx(:, :, 1+r:size(du_dx,3)) - du_dx(:, :, 1:size(du_dx,3)-r))), dim=1) &
         + ux_mean * sum(dv * (4.0_rk * 0.25_rk * (dv_dx(:, :, 1+r:size(dv_dx,3)) - dv_dx(:, :, 1:size(dv_dx,3)-r))), dim=1) &
         + ux_mean * sum(dw * (4.0_rk * 0.25_rk * (dw_dx(:, :, 1+r:size(dw_dx,3)) - dw_dx(:, :, 1:size(dw_dx,3)-r))), dim=1))
    s2 = sum(uy_mean * sum(du * (4.0_rk * 0.25_rk * (du_dy(:, :, 1+r:size(du_dy,3)) - du_dy(:, :, 1:size(du_dy,3)-r))), dim=1) &
         + uy_mean * sum(dv * (4.0_rk * 0.25_rk * (dv_dy(:, :, 1+r:size(dv_dy,3)) - dv_dy(:, :, 1:size(dv_dy,3)-r))), dim=1) &
         + uy_mean * sum(dw * (4.0_rk * 0.25_rk * (dw_dy(:, :, 1+r:size(dw_dy,3)) - dw_dy(:, :, 1:size(dw_dy,3)-r))), dim=1))
    s3 = sum(mean_du * sum(du * (0.5_rk * ddx) + dv * (0.5_rk * vdx) + dw * (0.5_rk * wdx), dim=1))
    s4 = sum(mean_dv * sum(du * (0.5_rk * ddy) + dv * (0.5_rk * vdy) + dw * (0.5_rk * wdy), dim=1))

    acc%ux_dS2_dXx_rx(ir) = acc%ux_dS2_dXx_rx(ir) + s1
    acc%uy_dS2_dXy_rx(ir) = acc%uy_dS2_dXy_rx(ir) + s2
    acc%dux_dS2_drx_rx(ir) = acc%dux_dS2_drx_rx(ir) + s3
    acc%duy_dS2_dry_rx(ir) = acc%duy_dS2_dry_rx(ir) + s4

    c = int(size(du,1), i8) * int(size(du,2), i8) * int(size(du,3), i8)
  end subroutine accumulate_one_rx

  subroutine accumulate_one_ry(ir, r, u, v, w, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, ctx, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(rk), intent(in) :: du_dx(:, :, :), du_dy(:, :, :), dv_dx(:, :, :), dv_dy(:, :, :), dw_dx(:, :, :), dw_dy(:, :, :)
    type(linear_transport_context_t), intent(in) :: ctx
    type(linear_transport_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: mean_du(:, :), mean_dv(:, :), ux_mean(:, :), uy_mean(:, :)
    real(rk), allocatable :: ddx(:, :, :), ddy(:, :, :), vdx(:, :, :), vdy(:, :, :), wdx(:, :, :), wdy(:, :, :)
    real(rk), allocatable :: tmp1(:, :), tmp2(:, :)
    real(rk) :: s1, s2, s3, s4

    allocate(du(size(u,1), size(u,2)-r, size(u,3)))
    allocate(dv(size(v,1), size(v,2)-r, size(v,3)))
    allocate(dw(size(w,1), size(w,2)-r, size(w,3)))
    allocate(mean_du(size(ctx%mean_u,1)-r, size(ctx%mean_u,2)))
    allocate(mean_dv(size(ctx%mean_v,1)-r, size(ctx%mean_v,2)))
    allocate(ux_mean(size(ctx%mean_u,1)-r, size(ctx%mean_u,2)))
    allocate(uy_mean(size(ctx%mean_v,1)-r, size(ctx%mean_v,2)))
    allocate(ddx(size(du_dx,1), size(du_dx,2)-r, size(du_dx,3)))
    allocate(ddy(size(du_dy,1), size(du_dy,2)-r, size(du_dy,3)))
    allocate(vdx(size(dv_dx,1), size(dv_dx,2)-r, size(dv_dx,3)))
    allocate(vdy(size(dv_dy,1), size(dv_dy,2)-r, size(dv_dy,3)))
    allocate(wdx(size(dw_dx,1), size(dw_dx,2)-r, size(dw_dx,3)))
    allocate(wdy(size(dw_dy,1), size(dw_dy,2)-r, size(dw_dy,3)))

    du = u(:, 1+r:size(u,2), :) - u(:, 1:size(u,2)-r, :)
    dv = v(:, 1+r:size(v,2), :) - v(:, 1:size(v,2)-r, :)
    dw = w(:, 1+r:size(w,2), :) - w(:, 1:size(w,2)-r, :)

    mean_du = ctx%mean_u(1+r:size(ctx%mean_u,1), :) - ctx%mean_u(1:size(ctx%mean_u,1)-r, :)
    mean_dv = ctx%mean_v(1+r:size(ctx%mean_v,1), :) - ctx%mean_v(1:size(ctx%mean_v,1)-r, :)
    ux_mean = 0.5_rk * (ctx%mean_u(1+r:size(ctx%mean_u,1), :) + ctx%mean_u(1:size(ctx%mean_u,1)-r, :))
    uy_mean = 0.5_rk * (ctx%mean_v(1+r:size(ctx%mean_v,1), :) + ctx%mean_v(1:size(ctx%mean_v,1)-r, :))

    ddx = du_dx(:, 1+r:size(du_dx,2), :) + du_dx(:, 1:size(du_dx,2)-r, :)
    ddy = du_dy(:, 1+r:size(du_dy,2), :) + du_dy(:, 1:size(du_dy,2)-r, :)
    vdx = dv_dx(:, 1+r:size(dv_dx,2), :) + dv_dx(:, 1:size(dv_dx,2)-r, :)
    vdy = dv_dy(:, 1+r:size(dv_dy,2), :) + dv_dy(:, 1:size(dv_dy,2)-r, :)
    wdx = dw_dx(:, 1+r:size(dw_dx,2), :) + dw_dx(:, 1:size(dw_dx,2)-r, :)
    wdy = dw_dy(:, 1+r:size(dw_dy,2), :) + dw_dy(:, 1:size(dw_dy,2)-r, :)

    tmp1 = sum(du * (du_dx(:, 1+r:size(du_dx,2), :) - du_dx(:, 1:size(du_dx,2)-r, :)) &
         + dv * (dv_dx(:, 1+r:size(dv_dx,2), :) - dv_dx(:, 1:size(dv_dx,2)-r, :)) &
         + dw * (dw_dx(:, 1+r:size(dw_dx,2), :) - dw_dx(:, 1:size(dw_dx,2)-r, :)), dim=1)
    s1 = sum(ux_mean * tmp1)

    tmp1 = sum(du * (du_dy(:, 1+r:size(du_dy,2), :) - du_dy(:, 1:size(du_dy,2)-r, :)) &
         + dv * (dv_dy(:, 1+r:size(dv_dy,2), :) - dv_dy(:, 1:size(dv_dy,2)-r, :)) &
         + dw * (dw_dy(:, 1+r:size(dw_dy,2), :) - dw_dy(:, 1:size(dw_dy,2)-r, :)), dim=1)
    s2 = sum(uy_mean * tmp1)

    tmp1 = sum(du * ddx + dv * vdx + dw * wdx, dim=1)
    tmp2 = sum(du * ddy + dv * vdy + dw * wdy, dim=1)
    s3 = 0.5_rk * sum(mean_du * tmp1)
    s4 = 0.5_rk * sum(mean_dv * tmp2)

    acc%ux_dS2_dXx_ry(ir) = acc%ux_dS2_dXx_ry(ir) + s1
    acc%uy_dS2_dXy_ry(ir) = acc%uy_dS2_dXy_ry(ir) + s2
    acc%dux_dS2_drx_ry(ir) = acc%dux_dS2_drx_ry(ir) + s3
    acc%duy_dS2_dry_ry(ir) = acc%duy_dS2_dry_ry(ir) + s4

    c = int(size(du,1), i8) * int(size(du,2), i8) * int(size(du,3), i8)
  end subroutine accumulate_one_ry

  subroutine finalize_linear_transport(acc, out)
    type(linear_transport_accumulator_t), intent(in) :: acc
    type(linear_transport_outputs_t), intent(out) :: out

    allocate(out%ux_dS2_dXx_rx(size(acc%ux_dS2_dXx_rx)))
    allocate(out%uy_dS2_dXy_rx(size(acc%uy_dS2_dXy_rx)))
    allocate(out%dux_dS2_drx_rx(size(acc%dux_dS2_drx_rx)))
    allocate(out%duy_dS2_dry_rx(size(acc%duy_dS2_dry_rx)))
    allocate(out%ux_dS2_dXx_ry(size(acc%ux_dS2_dXx_ry)))
    allocate(out%uy_dS2_dXy_ry(size(acc%uy_dS2_dXy_ry)))
    allocate(out%dux_dS2_drx_ry(size(acc%dux_dS2_drx_ry)))
    allocate(out%duy_dS2_dry_ry(size(acc%duy_dS2_dry_ry)))

    out%ux_dS2_dXx_rx = acc%ux_dS2_dXx_rx / real(acc%count_rx, rk)
    out%uy_dS2_dXy_rx = acc%uy_dS2_dXy_rx / real(acc%count_rx, rk)
    out%dux_dS2_drx_rx = acc%dux_dS2_drx_rx / real(acc%count_rx, rk)
    out%duy_dS2_dry_rx = acc%duy_dS2_dry_rx / real(acc%count_rx, rk)
    out%ux_dS2_dXx_ry = acc%ux_dS2_dXx_ry / real(acc%count_ry, rk)
    out%uy_dS2_dXy_ry = acc%uy_dS2_dXy_ry / real(acc%count_ry, rk)
    out%dux_dS2_drx_ry = acc%dux_dS2_drx_ry / real(acc%count_ry, rk)
    out%duy_dS2_dry_ry = acc%duy_dS2_dry_ry / real(acc%count_ry, rk)
  end subroutine finalize_linear_transport

end module mod_linear_transport_terms
