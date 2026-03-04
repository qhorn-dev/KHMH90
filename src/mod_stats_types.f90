module mod_stats_types
  use mod_kinds, only: i4, i8, rk
  implicit none
  private

  ! Generic storage for one family of terms that depends on separation r and
  ! keeps a transverse coordinate unresolved.
  type, public :: transverse_accumulator_t
    real(rk), allocatable :: sum_u2(:, :), sum_v2(:, :), sum_w2(:, :)
    real(rk), allocatable :: sum_u3(:, :), sum_v3(:, :), sum_w3(:, :)
    integer(i8), allocatable :: count(:, :)
  end type transverse_accumulator_t

  ! Final statistics for one separation direction.
  type, public :: transverse_outputs_t
    real(rk), allocatable :: s2u(:), s2v(:), s2w(:)
    real(rk), allocatable :: s3u(:), s3v(:), s3w(:)
    real(rk), allocatable :: s2u_t(:, :), s2v_t(:, :), s2w_t(:, :)
    real(rk), allocatable :: s3u_t(:, :), s3v_t(:, :), s3w_t(:, :)
  end type transverse_outputs_t

  ! Project-level container: rx keeps y unresolved, ry keeps x unresolved.
  type, public :: sf_accumulator_t
    type(transverse_accumulator_t) :: rx
    type(transverse_accumulator_t) :: ry
  end type sf_accumulator_t

  type, public :: sf_outputs_t
    type(transverse_outputs_t) :: rx
    type(transverse_outputs_t) :: ry
  end type sf_outputs_t

  public :: init_transverse_accumulator
  public :: finalize_transverse_accumulator

contains

  subroutine init_transverse_accumulator(acc, nr, nt)
    type(transverse_accumulator_t), intent(inout) :: acc
    integer(i4), intent(in) :: nr, nt

    allocate(acc%sum_u2(nr, nt), acc%sum_v2(nr, nt), acc%sum_w2(nr, nt))
    allocate(acc%sum_u3(nr, nt), acc%sum_v3(nr, nt), acc%sum_w3(nr, nt))
    allocate(acc%count(nr, nt))

    acc%sum_u2 = 0.0_rk
    acc%sum_v2 = 0.0_rk
    acc%sum_w2 = 0.0_rk
    acc%sum_u3 = 0.0_rk
    acc%sum_v3 = 0.0_rk
    acc%sum_w3 = 0.0_rk
    acc%count = 0_i8
  end subroutine init_transverse_accumulator

  subroutine finalize_transverse_accumulator(acc, out)
    type(transverse_accumulator_t), intent(in) :: acc
    type(transverse_outputs_t), intent(out) :: out

    integer(i4) :: nr, nt, ir, it
    integer(i8) :: ctot

    nr = size(acc%count, 1)
    nt = size(acc%count, 2)

    allocate(out%s2u(nr), out%s2v(nr), out%s2w(nr))
    allocate(out%s3u(nr), out%s3v(nr), out%s3w(nr))
    allocate(out%s2u_t(nr, nt), out%s2v_t(nr, nt), out%s2w_t(nr, nt))
    allocate(out%s3u_t(nr, nt), out%s3v_t(nr, nt), out%s3w_t(nr, nt))

    do ir = 1, nr
      do it = 1, nt
        out%s2u_t(ir, it) = acc%sum_u2(ir, it) / real(acc%count(ir, it), rk)
        out%s2v_t(ir, it) = acc%sum_v2(ir, it) / real(acc%count(ir, it), rk)
        out%s2w_t(ir, it) = acc%sum_w2(ir, it) / real(acc%count(ir, it), rk)
        out%s3u_t(ir, it) = acc%sum_u3(ir, it) / real(acc%count(ir, it), rk)
        out%s3v_t(ir, it) = acc%sum_v3(ir, it) / real(acc%count(ir, it), rk)
        out%s3w_t(ir, it) = acc%sum_w3(ir, it) / real(acc%count(ir, it), rk)
      end do

      ctot = sum(acc%count(ir, :))
      out%s2u(ir) = sum(acc%sum_u2(ir, :)) / real(ctot, rk)
      out%s2v(ir) = sum(acc%sum_v2(ir, :)) / real(ctot, rk)
      out%s2w(ir) = sum(acc%sum_w2(ir, :)) / real(ctot, rk)
      out%s3u(ir) = sum(acc%sum_u3(ir, :)) / real(ctot, rk)
      out%s3v(ir) = sum(acc%sum_v3(ir, :)) / real(ctot, rk)
      out%s3w(ir) = sum(acc%sum_w3(ir, :)) / real(ctot, rk)
    end do
  end subroutine finalize_transverse_accumulator

end module mod_stats_types
