module mod_pdelta_terms
  use mod_kinds, only: i4, i8, rk
  use mod_grid, only: grid_info_t
  use mod_linear_transport_terms, only: linear_transport_context_t
  implicit none
  private

  type, public :: pdelta_accumulator_t
    real(rk), allocatable :: p111_rx(:), p121_rx(:), p212_rx(:), p222_rx(:), p313_rx(:), p323_rx(:), p_rx(:)
    real(rk), allocatable :: p111_ry(:), p121_ry(:), p212_ry(:), p222_ry(:), p313_ry(:), p323_ry(:), p_ry(:)
    integer(i8), allocatable :: count_rx(:), count_ry(:)
  end type

  type, public :: pdelta_outputs_t
    real(rk), allocatable :: p111_rx(:), p121_rx(:), p212_rx(:), p222_rx(:), p313_rx(:), p323_rx(:), p_rx(:)
    real(rk), allocatable :: p111_ry(:), p121_ry(:), p212_ry(:), p222_ry(:), p313_ry(:), p323_ry(:), p_ry(:)
  end type

  public :: init_pdelta_accumulator
  public :: accumulate_pdelta_block
  public :: finalize_pdelta

contains

  subroutine init_pdelta_accumulator(acc, nrx, nry)
    type(pdelta_accumulator_t), intent(inout) :: acc
    integer(i4), intent(in) :: nrx, nry

    allocate(acc%p111_rx(nrx), acc%p121_rx(nrx), acc%p212_rx(nrx), acc%p222_rx(nrx))
    allocate(acc%p313_rx(nrx), acc%p323_rx(nrx), acc%p_rx(nrx))
    allocate(acc%p111_ry(nry), acc%p121_ry(nry), acc%p212_ry(nry), acc%p222_ry(nry))
    allocate(acc%p313_ry(nry), acc%p323_ry(nry), acc%p_ry(nry))
    allocate(acc%count_rx(nrx), acc%count_ry(nry))

    acc%p111_rx = 0.0_rk
    acc%p121_rx = 0.0_rk
    acc%p212_rx = 0.0_rk
    acc%p222_rx = 0.0_rk
    acc%p313_rx = 0.0_rk
    acc%p323_rx = 0.0_rk
    acc%p_rx = 0.0_rk
    acc%p111_ry = 0.0_rk
    acc%p121_ry = 0.0_rk
    acc%p212_ry = 0.0_rk
    acc%p222_ry = 0.0_rk
    acc%p313_ry = 0.0_rk
    acc%p323_ry = 0.0_rk
    acc%p_ry = 0.0_rk
    acc%count_rx = 0_i8
    acc%count_ry = 0_i8
  end subroutine init_pdelta_accumulator

  subroutine accumulate_pdelta_block(u, v, w, grid, rx_pts, ry_pts, ctx, acc)
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    type(grid_info_t), intent(in) :: grid
    integer(i4), intent(in) :: rx_pts(:), ry_pts(:)
    type(linear_transport_context_t), intent(in) :: ctx
    type(pdelta_accumulator_t), intent(inout) :: acc

    integer(i4) :: ir, r
    integer(i8) :: c

    !$omp parallel do default(none) shared(rx_pts,u,v,w,grid,ctx,acc) private(ir,r,c) schedule(static)
    do ir = 1, size(rx_pts)
      r = rx_pts(ir)
      call accumulate_one_rx(ir, r, u, v, w, grid, ctx, acc, c)
      acc%count_rx(ir) = acc%count_rx(ir) + c
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(ry_pts,u,v,w,grid,ctx,acc) private(ir,r,c) schedule(static)
    do ir = 1, size(ry_pts)
      r = ry_pts(ir)
      call accumulate_one_ry(ir, r, u, v, w, grid, ctx, acc, c)
      acc%count_ry(ir) = acc%count_ry(ir) + c
    end do
    !$omp end parallel do
  end subroutine accumulate_pdelta_block

  subroutine accumulate_one_rx(ir, r, u, v, w, grid, ctx, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    type(grid_info_t), intent(in) :: grid
    type(linear_transport_context_t), intent(in) :: ctx
    type(pdelta_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: uc(:, :, :), vc(:, :, :), wc(:, :, :)
    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: gux(:, :), guy(:, :), gvx(:, :), gvy(:, :), gwx(:, :), gwy(:, :)
    real(rk) :: p111, p121, p212, p222, p313, p323
    integer(i4) :: nt, ny, nx

    nt = size(u, 1)
    ny = grid%ny - 2
    nx = grid%nx - 2

    allocate(uc(nt, ny, nx), vc(nt, ny, nx), wc(nt, ny, nx))
    uc = u(:, 2:grid%ny-1, 2:grid%nx-1)
    vc = v(:, 2:grid%ny-1, 2:grid%nx-1)
    wc = w(:, 2:grid%ny-1, 2:grid%nx-1)

    allocate(du(nt, ny, nx-r), dv(nt, ny, nx-r), dw(nt, ny, nx-r))
    du = uc(:, :, 1+r:nx) - uc(:, :, 1:nx-r)
    dv = vc(:, :, 1+r:nx) - vc(:, :, 1:nx-r)
    dw = wc(:, :, 1+r:nx) - wc(:, :, 1:nx-r)

    allocate(gux(ny, nx-r), guy(ny, nx-r), gvx(ny, nx-r), gvy(ny, nx-r), gwx(ny, nx-r), gwy(ny, nx-r))
    gux = 0.5_rk * (ctx%dmean_u_dx(:, 1+r:nx) + ctx%dmean_u_dx(:, 1:nx-r))
    guy = 0.5_rk * (ctx%dmean_u_dy(:, 1+r:nx) + ctx%dmean_u_dy(:, 1:nx-r))
    gvx = 0.5_rk * (ctx%dmean_v_dx(:, 1+r:nx) + ctx%dmean_v_dx(:, 1:nx-r))
    gvy = 0.5_rk * (ctx%dmean_v_dy(:, 1+r:nx) + ctx%dmean_v_dy(:, 1:nx-r))
    gwx = 0.5_rk * (ctx%dmean_w_dx(:, 1+r:nx) + ctx%dmean_w_dx(:, 1:nx-r))
    gwy = 0.5_rk * (ctx%dmean_w_dy(:, 1+r:nx) + ctx%dmean_w_dy(:, 1:nx-r))

    p111 = -2.0_rk * sum(du * du * spread(gux, 1, nt))
    p121 = -2.0_rk * sum(du * dv * spread(guy, 1, nt))
    p212 = -2.0_rk * sum(dv * du * spread(gvx, 1, nt))
    p222 = -2.0_rk * sum(dv * dv * spread(gvy, 1, nt))
    p313 = -2.0_rk * sum(dw * du * spread(gwx, 1, nt))
    p323 = -2.0_rk * sum(dw * dv * spread(gwy, 1, nt))

    acc%p111_rx(ir) = acc%p111_rx(ir) + p111
    acc%p121_rx(ir) = acc%p121_rx(ir) + p121
    acc%p212_rx(ir) = acc%p212_rx(ir) + p212
    acc%p222_rx(ir) = acc%p222_rx(ir) + p222
    acc%p313_rx(ir) = acc%p313_rx(ir) + p313
    acc%p323_rx(ir) = acc%p323_rx(ir) + p323
    acc%p_rx(ir) = acc%p_rx(ir) + p111 + p121 + p212 + p222 + p313 + p323

    c = int(nt, i8) * int(ny, i8) * int(nx-r, i8)
  end subroutine accumulate_one_rx

  subroutine accumulate_one_ry(ir, r, u, v, w, grid, ctx, acc, c)
    integer(i4), intent(in) :: ir, r
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    type(grid_info_t), intent(in) :: grid
    type(linear_transport_context_t), intent(in) :: ctx
    type(pdelta_accumulator_t), intent(inout) :: acc
    integer(i8), intent(out) :: c

    real(rk), allocatable :: uc(:, :, :), vc(:, :, :), wc(:, :, :)
    real(rk), allocatable :: du(:, :, :), dv(:, :, :), dw(:, :, :)
    real(rk), allocatable :: gux(:, :), guy(:, :), gvx(:, :), gvy(:, :), gwx(:, :), gwy(:, :)
    real(rk) :: p111, p121, p212, p222, p313, p323
    integer(i4) :: nt, ny, nx

    nt = size(u, 1)
    ny = grid%ny - 2
    nx = grid%nx - 2

    allocate(uc(nt, ny, nx), vc(nt, ny, nx), wc(nt, ny, nx))
    uc = u(:, 2:grid%ny-1, 2:grid%nx-1)
    vc = v(:, 2:grid%ny-1, 2:grid%nx-1)
    wc = w(:, 2:grid%ny-1, 2:grid%nx-1)

    allocate(du(nt, ny-r, nx), dv(nt, ny-r, nx), dw(nt, ny-r, nx))
    du = uc(:, 1+r:ny, :) - uc(:, 1:ny-r, :)
    dv = vc(:, 1+r:ny, :) - vc(:, 1:ny-r, :)
    dw = wc(:, 1+r:ny, :) - wc(:, 1:ny-r, :)

    allocate(gux(ny-r, nx), guy(ny-r, nx), gvx(ny-r, nx), gvy(ny-r, nx), gwx(ny-r, nx), gwy(ny-r, nx))
    gux = 0.5_rk * (ctx%dmean_u_dx(1+r:ny, :) + ctx%dmean_u_dx(1:ny-r, :))
    guy = 0.5_rk * (ctx%dmean_u_dy(1+r:ny, :) + ctx%dmean_u_dy(1:ny-r, :))
    gvx = 0.5_rk * (ctx%dmean_v_dx(1+r:ny, :) + ctx%dmean_v_dx(1:ny-r, :))
    gvy = 0.5_rk * (ctx%dmean_v_dy(1+r:ny, :) + ctx%dmean_v_dy(1:ny-r, :))
    gwx = 0.5_rk * (ctx%dmean_w_dx(1+r:ny, :) + ctx%dmean_w_dx(1:ny-r, :))
    gwy = 0.5_rk * (ctx%dmean_w_dy(1+r:ny, :) + ctx%dmean_w_dy(1:ny-r, :))

    p111 = -2.0_rk * sum(du * du * spread(gux, 1, nt))
    p121 = -2.0_rk * sum(du * dv * spread(guy, 1, nt))
    p212 = -2.0_rk * sum(dv * du * spread(gvx, 1, nt))
    p222 = -2.0_rk * sum(dv * dv * spread(gvy, 1, nt))
    p313 = -2.0_rk * sum(dw * du * spread(gwx, 1, nt))
    p323 = -2.0_rk * sum(dw * dv * spread(gwy, 1, nt))

    acc%p111_ry(ir) = acc%p111_ry(ir) + p111
    acc%p121_ry(ir) = acc%p121_ry(ir) + p121
    acc%p212_ry(ir) = acc%p212_ry(ir) + p212
    acc%p222_ry(ir) = acc%p222_ry(ir) + p222
    acc%p313_ry(ir) = acc%p313_ry(ir) + p313
    acc%p323_ry(ir) = acc%p323_ry(ir) + p323
    acc%p_ry(ir) = acc%p_ry(ir) + p111 + p121 + p212 + p222 + p313 + p323

    c = int(nt, i8) * int(ny-r, i8) * int(nx, i8)
  end subroutine accumulate_one_ry

  subroutine finalize_pdelta(acc, out)
    type(pdelta_accumulator_t), intent(in) :: acc
    type(pdelta_outputs_t), intent(out) :: out

    allocate(out%p111_rx(size(acc%p111_rx)), out%p121_rx(size(acc%p121_rx)))
    allocate(out%p212_rx(size(acc%p212_rx)), out%p222_rx(size(acc%p222_rx)))
    allocate(out%p313_rx(size(acc%p313_rx)), out%p323_rx(size(acc%p323_rx)), out%p_rx(size(acc%p_rx)))
    allocate(out%p111_ry(size(acc%p111_ry)), out%p121_ry(size(acc%p121_ry)))
    allocate(out%p212_ry(size(acc%p212_ry)), out%p222_ry(size(acc%p222_ry)))
    allocate(out%p313_ry(size(acc%p313_ry)), out%p323_ry(size(acc%p323_ry)), out%p_ry(size(acc%p_ry)))

    out%p111_rx = acc%p111_rx / real(acc%count_rx, rk)
    out%p121_rx = acc%p121_rx / real(acc%count_rx, rk)
    out%p212_rx = acc%p212_rx / real(acc%count_rx, rk)
    out%p222_rx = acc%p222_rx / real(acc%count_rx, rk)
    out%p313_rx = acc%p313_rx / real(acc%count_rx, rk)
    out%p323_rx = acc%p323_rx / real(acc%count_rx, rk)
    out%p_rx = acc%p_rx / real(acc%count_rx, rk)

    out%p111_ry = acc%p111_ry / real(acc%count_ry, rk)
    out%p121_ry = acc%p121_ry / real(acc%count_ry, rk)
    out%p212_ry = acc%p212_ry / real(acc%count_ry, rk)
    out%p222_ry = acc%p222_ry / real(acc%count_ry, rk)
    out%p313_ry = acc%p313_ry / real(acc%count_ry, rk)
    out%p323_ry = acc%p323_ry / real(acc%count_ry, rk)
    out%p_ry = acc%p_ry / real(acc%count_ry, rk)
  end subroutine finalize_pdelta

end module mod_pdelta_terms
