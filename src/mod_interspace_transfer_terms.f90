module mod_interspace_transfer_terms
  use mod_kinds, only: i4, i8, rk
  use mod_grid, only: grid_info_t
  implicit none
  private

  type, public :: interspace_transfer_accumulator_t
    real(rk), allocatable :: divXx_rx(:), divXy_rx(:), divX_rx(:)
    real(rk), allocatable :: divXx_ry(:), divXy_ry(:), divX_ry(:)
    integer(i8), allocatable :: count_rx(:), count_ry(:)
  end type

  type, public :: interspace_transfer_outputs_t
    real(rk), allocatable :: divXx_rx(:), divXy_rx(:), divX_rx(:)
    real(rk), allocatable :: divXx_ry(:), divXy_ry(:), divX_ry(:)
  end type

  public :: init_interspace_transfer_accumulator
  public :: accumulate_interspace_transfer_block
  public :: finalize_interspace_transfer

contains

  subroutine central_diff_3d_x(f, dx, df)
    real(rk), intent(in) :: f(:, :, :)
    real(rk), intent(in) :: dx
    real(rk), allocatable, intent(out) :: df(:, :, :)

    allocate(df(size(f, 1), size(f, 2) - 2, size(f, 3) - 2))
    df = (f(:, 2:size(f,2)-1, 3:size(f,3)) - f(:, 2:size(f,2)-1, 1:size(f,3)-2)) / (2.0_rk * dx)
  end subroutine central_diff_3d_x

  subroutine central_diff_3d_y(f, dy, df)
    real(rk), intent(in) :: f(:, :, :)
    real(rk), intent(in) :: dy
    real(rk), allocatable, intent(out) :: df(:, :, :)

    allocate(df(size(f, 1), size(f, 2) - 2, size(f, 3) - 2))
    df = (f(:, 3:size(f,2), 2:size(f,3)-1) - f(:, 1:size(f,2)-2, 2:size(f,3)-1)) / (2.0_rk * dy)
  end subroutine central_diff_3d_y

  subroutine init_interspace_transfer_accumulator(acc, nrx, nry)
    type(interspace_transfer_accumulator_t), intent(inout) :: acc
    integer(i4), intent(in) :: nrx, nry

    allocate(acc%divXx_rx(nrx), acc%divXy_rx(nrx), acc%divX_rx(nrx))
    allocate(acc%divXx_ry(nry), acc%divXy_ry(nry), acc%divX_ry(nry))
    allocate(acc%count_rx(nrx), acc%count_ry(nry))

    acc%divXx_rx = 0.0_rk
    acc%divXy_rx = 0.0_rk
    acc%divX_rx = 0.0_rk
    acc%divXx_ry = 0.0_rk
    acc%divXy_ry = 0.0_rk
    acc%divX_ry = 0.0_rk
    acc%count_rx = 0_i8
    acc%count_ry = 0_i8
  end subroutine init_interspace_transfer_accumulator

  subroutine accumulate_interspace_transfer_block(u, v, w, grid, rx_pts, ry_pts, acc)
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    type(grid_info_t), intent(in) :: grid
    integer(i4), intent(in) :: rx_pts(:), ry_pts(:)
    type(interspace_transfer_accumulator_t), intent(inout) :: acc

    real(rk), allocatable :: uc(:, :, :), vc(:, :, :), wc(:, :, :)
    real(rk), allocatable :: du_dx(:, :, :), du_dy(:, :, :)
    real(rk), allocatable :: dv_dx(:, :, :), dv_dy(:, :, :)
    real(rk), allocatable :: dw_dx(:, :, :), dw_dy(:, :, :)
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
  end subroutine accumulate_interspace_transfer_block

  subroutine accumulate_one_rx(ir, r, u, v, w, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(rk), intent(in) :: du_dx(:, :, :), du_dy(:, :, :)
    real(rk), intent(in) :: dv_dx(:, :, :), dv_dy(:, :, :)
    real(rk), intent(in) :: dw_dx(:, :, :), dw_dy(:, :, :)
    type(interspace_transfer_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: ux(:, :, :), uy(:, :, :)
    real(rk), allocatable :: dXx_u(:, :, :), dXx_v(:, :, :), dXx_w(:, :, :)
    real(rk), allocatable :: dXy_u(:, :, :), dXy_v(:, :, :), dXy_w(:, :, :)
    real(rk) :: s1, s2

    allocate(du(size(u,1), size(u,2), size(u,3)-r))
    allocate(dv(size(v,1), size(v,2), size(v,3)-r))
    allocate(dw(size(w,1), size(w,2), size(w,3)-r))
    allocate(ux(size(u,1), size(u,2), size(u,3)-r))
    allocate(uy(size(v,1), size(v,2), size(v,3)-r))
    allocate(dXx_u(size(u,1), size(u,2), size(u,3)-r))
    allocate(dXx_v(size(v,1), size(v,2), size(v,3)-r))
    allocate(dXx_w(size(w,1), size(w,2), size(w,3)-r))
    allocate(dXy_u(size(u,1), size(u,2), size(u,3)-r))
    allocate(dXy_v(size(v,1), size(v,2), size(v,3)-r))
    allocate(dXy_w(size(w,1), size(w,2), size(w,3)-r))

    du = u(:, :, 1+r:size(u,3)) - u(:, :, 1:size(u,3)-r)
    dv = v(:, :, 1+r:size(v,3)) - v(:, :, 1:size(v,3)-r)
    dw = w(:, :, 1+r:size(w,3)) - w(:, :, 1:size(w,3)-r)
    ux = 0.5_rk * (u(:, :, 1+r:size(u,3)) + u(:, :, 1:size(u,3)-r))
    uy = 0.5_rk * (v(:, :, 1+r:size(v,3)) + v(:, :, 1:size(v,3)-r))

    dXx_u = du_dx(:, :, 1+r:size(du_dx,3)) - du_dx(:, :, 1:size(du_dx,3)-r)
    dXx_v = dv_dx(:, :, 1+r:size(dv_dx,3)) - dv_dx(:, :, 1:size(dv_dx,3)-r)
    dXx_w = dw_dx(:, :, 1+r:size(dw_dx,3)) - dw_dx(:, :, 1:size(dw_dx,3)-r)
    dXy_u = du_dy(:, :, 1+r:size(du_dy,3)) - du_dy(:, :, 1:size(du_dy,3)-r)
    dXy_v = dv_dy(:, :, 1+r:size(dv_dy,3)) - dv_dy(:, :, 1:size(dv_dy,3)-r)
    dXy_w = dw_dy(:, :, 1+r:size(dw_dy,3)) - dw_dy(:, :, 1:size(dw_dy,3)-r)

    s1 = 2.0_rk * sum(ux * du * dXx_u + ux * dv * dXx_v + ux * dw * dXx_w)
    s2 = 2.0_rk * sum(uy * du * dXy_u + uy * dv * dXy_v + uy * dw * dXy_w)

    acc%divXx_rx(ir) = acc%divXx_rx(ir) + s1
    acc%divXy_rx(ir) = acc%divXy_rx(ir) + s2
    acc%divX_rx(ir) = acc%divX_rx(ir) + s1 + s2

    c = int(size(du,1), i8) * int(size(du,2), i8) * int(size(du,3), i8)
  end subroutine accumulate_one_rx

  subroutine accumulate_one_ry(ir, r, u, v, w, du_dx, du_dy, dv_dx, dv_dy, dw_dx, dw_dy, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(rk), intent(in) :: du_dx(:, :, :), du_dy(:, :, :)
    real(rk), intent(in) :: dv_dx(:, :, :), dv_dy(:, :, :)
    real(rk), intent(in) :: dw_dx(:, :, :), dw_dy(:, :, :)
    type(interspace_transfer_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: ux(:, :, :), uy(:, :, :)
    real(rk), allocatable :: dXx_u(:, :, :), dXx_v(:, :, :), dXx_w(:, :, :)
    real(rk), allocatable :: dXy_u(:, :, :), dXy_v(:, :, :), dXy_w(:, :, :)
    real(rk) :: s1, s2

    allocate(du(size(u,1), size(u,2)-r, size(u,3)))
    allocate(dv(size(v,1), size(v,2)-r, size(v,3)))
    allocate(dw(size(w,1), size(w,2)-r, size(w,3)))
    allocate(ux(size(u,1), size(u,2)-r, size(u,3)))
    allocate(uy(size(v,1), size(v,2)-r, size(v,3)))
    allocate(dXx_u(size(u,1), size(u,2)-r, size(u,3)))
    allocate(dXx_v(size(v,1), size(v,2)-r, size(v,3)))
    allocate(dXx_w(size(w,1), size(w,2)-r, size(w,3)))
    allocate(dXy_u(size(u,1), size(u,2)-r, size(u,3)))
    allocate(dXy_v(size(v,1), size(v,2)-r, size(v,3)))
    allocate(dXy_w(size(w,1), size(w,2)-r, size(w,3)))

    du = u(:, 1+r:size(u,2), :) - u(:, 1:size(u,2)-r, :)
    dv = v(:, 1+r:size(v,2), :) - v(:, 1:size(v,2)-r, :)
    dw = w(:, 1+r:size(w,2), :) - w(:, 1:size(w,2)-r, :)
    ux = 0.5_rk * (u(:, 1+r:size(u,2), :) + u(:, 1:size(u,2)-r, :))
    uy = 0.5_rk * (v(:, 1+r:size(v,2), :) + v(:, 1:size(v,2)-r, :))

    dXx_u = du_dx(:, 1+r:size(du_dx,2), :) - du_dx(:, 1:size(du_dx,2)-r, :)
    dXx_v = dv_dx(:, 1+r:size(dv_dx,2), :) - dv_dx(:, 1:size(dv_dx,2)-r, :)
    dXx_w = dw_dx(:, 1+r:size(dw_dx,2), :) - dw_dx(:, 1:size(dw_dx,2)-r, :)
    dXy_u = du_dy(:, 1+r:size(du_dy,2), :) - du_dy(:, 1:size(du_dy,2)-r, :)
    dXy_v = dv_dy(:, 1+r:size(dv_dy,2), :) - dv_dy(:, 1:size(dv_dy,2)-r, :)
    dXy_w = dw_dy(:, 1+r:size(dw_dy,2), :) - dw_dy(:, 1:size(dw_dy,2)-r, :)

    s1 = 2.0_rk * sum(ux * du * dXx_u + ux * dv * dXx_v + ux * dw * dXx_w)
    s2 = 2.0_rk * sum(uy * du * dXy_u + uy * dv * dXy_v + uy * dw * dXy_w)

    acc%divXx_ry(ir) = acc%divXx_ry(ir) + s1
    acc%divXy_ry(ir) = acc%divXy_ry(ir) + s2
    acc%divX_ry(ir) = acc%divX_ry(ir) + s1 + s2

    c = int(size(du,1), i8) * int(size(du,2), i8) * int(size(du,3), i8)
  end subroutine accumulate_one_ry

  subroutine finalize_interspace_transfer(acc, out)
    type(interspace_transfer_accumulator_t), intent(in) :: acc
    type(interspace_transfer_outputs_t), intent(out) :: out

    allocate(out%divXx_rx(size(acc%divXx_rx)), out%divXy_rx(size(acc%divXy_rx)), out%divX_rx(size(acc%divX_rx)))
    allocate(out%divXx_ry(size(acc%divXx_ry)), out%divXy_ry(size(acc%divXy_ry)), out%divX_ry(size(acc%divX_ry)))

    out%divXx_rx = acc%divXx_rx / real(acc%count_rx, rk)
    out%divXy_rx = acc%divXy_rx / real(acc%count_rx, rk)
    out%divX_rx = acc%divX_rx / real(acc%count_rx, rk)
    out%divXx_ry = acc%divXx_ry / real(acc%count_ry, rk)
    out%divXy_ry = acc%divXy_ry / real(acc%count_ry, rk)
    out%divX_ry = acc%divX_ry / real(acc%count_ry, rk)
  end subroutine finalize_interspace_transfer

end module mod_interspace_transfer_terms
