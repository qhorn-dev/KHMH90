module mod_sf_terms
  use mod_kinds, only: i4, i8, rk
  use mod_grid, only: valid_i_bounds
  use mod_stats_types, only: sf_accumulator_t, sf_outputs_t, &
                             init_transverse_accumulator, finalize_transverse_accumulator
  implicit none
  private

  public :: init_sf_accumulator
  public :: accumulate_sf_block
  public :: finalize_sf_accumulator

contains

  subroutine init_sf_accumulator(acc, nrx, nry, ny, nx)
    type(sf_accumulator_t), intent(inout) :: acc
    integer(i4), intent(in) :: nrx, nry, ny, nx

    call init_transverse_accumulator(acc%rx, nrx, ny)
    call init_transverse_accumulator(acc%ry, nry, nx)
  end subroutine init_sf_accumulator

  subroutine accumulate_sf_block(u, v, w, rx_pts, ry_pts, acc)
    real(rk), intent(in) :: u(:, :, :), v(:, :, :), w(:, :, :)
    integer(i4), intent(in) :: rx_pts(:), ry_pts(:)
    type(sf_accumulator_t), intent(inout) :: acc

    integer(i4) :: nt, ny, nx
    integer(i4) :: ir, r, t, iy, ix, i0, i1
    real(rk) :: du, dv, dw

    nt = size(u, 1)
    ny = size(u, 2)
    nx = size(u, 3)
    if (size(v, 1) /= nt .or. size(v, 2) /= ny .or. size(v, 3) /= nx) stop 'inconsistent v dimensions'
    if (size(w, 1) /= nt .or. size(w, 2) /= ny .or. size(w, 3) /= nx) stop 'inconsistent w dimensions'

    !$omp parallel do default(none) shared(u,v,w,rx_pts,acc,nt,ny,nx) &
    !$omp private(ir,r,t,iy,ix,i0,i1,du,dv,dw) schedule(static)
    do ir = 1, size(rx_pts)
      r = rx_pts(ir)
      call valid_i_bounds(nx, r, i0, i1)
      if (i1 < i0) cycle
      do iy = 1, ny
        do t = 1, nt
          do ix = i0, i1
            du = u(t, iy, ix + r) - u(t, iy, ix)
            dv = v(t, iy, ix + r) - v(t, iy, ix)
            dw = w(t, iy, ix + r) - w(t, iy, ix)
            acc%rx%sum_u2(ir, iy) = acc%rx%sum_u2(ir, iy) + du * du
            acc%rx%sum_v2(ir, iy) = acc%rx%sum_v2(ir, iy) + dv * dv
            acc%rx%sum_w2(ir, iy) = acc%rx%sum_w2(ir, iy) + dw * dw
            acc%rx%sum_u3(ir, iy) = acc%rx%sum_u3(ir, iy) + du * du * du
            acc%rx%sum_v3(ir, iy) = acc%rx%sum_v3(ir, iy) + dv * dv * dv
            acc%rx%sum_w3(ir, iy) = acc%rx%sum_w3(ir, iy) + dw * dw * dw
            acc%rx%count(ir, iy) = acc%rx%count(ir, iy) + 1_i8
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(u,v,w,ry_pts,acc,nt,ny,nx) &
    !$omp private(ir,r,t,iy,ix,i0,i1,du,dv,dw) schedule(static)
    do ir = 1, size(ry_pts)
      r = ry_pts(ir)
      call valid_i_bounds(ny, r, i0, i1)
      if (i1 < i0) cycle
      do ix = 1, nx
        do t = 1, nt
          do iy = i0, i1
            du = u(t, iy + r, ix) - u(t, iy, ix)
            dv = v(t, iy + r, ix) - v(t, iy, ix)
            dw = w(t, iy + r, ix) - w(t, iy, ix)
            acc%ry%sum_u2(ir, ix) = acc%ry%sum_u2(ir, ix) + du * du
            acc%ry%sum_v2(ir, ix) = acc%ry%sum_v2(ir, ix) + dv * dv
            acc%ry%sum_w2(ir, ix) = acc%ry%sum_w2(ir, ix) + dw * dw
            acc%ry%sum_u3(ir, ix) = acc%ry%sum_u3(ir, ix) + du * du * du
            acc%ry%sum_v3(ir, ix) = acc%ry%sum_v3(ir, ix) + dv * dv * dv
            acc%ry%sum_w3(ir, ix) = acc%ry%sum_w3(ir, ix) + dw * dw * dw
            acc%ry%count(ir, ix) = acc%ry%count(ir, ix) + 1_i8
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine accumulate_sf_block

  subroutine finalize_sf_accumulator(acc, out)
    type(sf_accumulator_t), intent(in) :: acc
    type(sf_outputs_t), intent(out) :: out

    call finalize_transverse_accumulator(acc%rx, out%rx)
    call finalize_transverse_accumulator(acc%ry, out%ry)
  end subroutine finalize_sf_accumulator

end module mod_sf_terms
