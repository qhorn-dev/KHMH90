module mod_kinds
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  implicit none
  private

  integer, parameter, public :: i4 = int32
  integer, parameter, public :: i8 = int64
  integer, parameter, public :: rk = real64

end module mod_kinds
