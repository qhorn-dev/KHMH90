module mod_output_netcdf
  use mod_kinds, only: rk
  use mod_grid, only: grid_info_t
  use mod_stats_types, only: sf_outputs_t, transverse_outputs_t
  use mod_linear_transport_terms, only: linear_transport_outputs_t
  use mod_interscale_transfer_terms, only: interscale_transfer_outputs_t
  use mod_interspace_transfer_terms, only: interspace_transfer_outputs_t
  use mod_pdelta_terms, only: pdelta_outputs_t
  use mod_px_terms, only: px_outputs_t
  use netcdf
  implicit none
  private

  public :: write_sf_output

contains

  subroutine nc_check(status, where)
    integer, intent(in) :: status
    character(len=*), intent(in) :: where

    if (status /= nf90_noerr) then
      write(*,'(a)') 'NetCDF error in '//trim(where)//': '//trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

  subroutine define_var_1d(ncid, name, dimid, varid)
    integer, intent(in) :: ncid, dimid
    character(len=*), intent(in) :: name
    integer, intent(out) :: varid

    call nc_check( &
         nf90_def_var(ncid, trim(name), nf90_double, (/ dimid /), varid), &
         'def '//trim(name) &
    )
  end subroutine define_var_1d

  subroutine define_var_2d(ncid, name, dimid_a, dimid_b, varid)
    integer, intent(in) :: ncid, dimid_a, dimid_b
    character(len=*), intent(in) :: name
    integer, intent(out) :: varid

    call nc_check( &
         nf90_def_var(ncid, trim(name), nf90_double, (/ dimid_a, dimid_b /), varid), &
         'def '//trim(name) &
    )
  end subroutine define_var_2d

  subroutine define_direction_vars(ncid, prefix, suffix, dimid_r, dimid_t, vids)
    integer, intent(in) :: ncid, dimid_r, dimid_t
    character(len=*), intent(in) :: prefix, suffix
    integer, intent(out) :: vids(12)

    call define_var_1d(ncid, 'deltaux_2_'//trim(prefix), dimid_r, vids(1))
    call define_var_1d(ncid, 'deltauy_2_'//trim(prefix), dimid_r, vids(2))
    call define_var_1d(ncid, 'deltauz_2_'//trim(prefix), dimid_r, vids(3))
    call define_var_1d(ncid, 'deltaux_3_'//trim(prefix), dimid_r, vids(4))
    call define_var_1d(ncid, 'deltauy_3_'//trim(prefix), dimid_r, vids(5))
    call define_var_1d(ncid, 'deltauz_3_'//trim(prefix), dimid_r, vids(6))

    call define_var_2d(ncid, 'deltaux_2_'//trim(prefix)//'_'//trim(suffix), dimid_t, dimid_r, vids(7))
    call define_var_2d(ncid, 'deltauy_2_'//trim(prefix)//'_'//trim(suffix), dimid_t, dimid_r, vids(8))
    call define_var_2d(ncid, 'deltauz_2_'//trim(prefix)//'_'//trim(suffix), dimid_t, dimid_r, vids(9))
    call define_var_2d(ncid, 'deltaux_3_'//trim(prefix)//'_'//trim(suffix), dimid_t, dimid_r, vids(10))
    call define_var_2d(ncid, 'deltauy_3_'//trim(prefix)//'_'//trim(suffix), dimid_t, dimid_r, vids(11))
    call define_var_2d(ncid, 'deltauz_3_'//trim(prefix)//'_'//trim(suffix), dimid_t, dimid_r, vids(12))
  end subroutine define_direction_vars

  subroutine write_direction_data(ncid, vids, out)
    integer, intent(in) :: ncid, vids(12)
    type(transverse_outputs_t), intent(in) :: out

    call nc_check(nf90_put_var(ncid, vids(1), out%s2u), 'write s2u')
    call nc_check(nf90_put_var(ncid, vids(2), out%s2v), 'write s2v')
    call nc_check(nf90_put_var(ncid, vids(3), out%s2w), 'write s2w')
    call nc_check(nf90_put_var(ncid, vids(4), out%s3u), 'write s3u')
    call nc_check(nf90_put_var(ncid, vids(5), out%s3v), 'write s3v')
    call nc_check(nf90_put_var(ncid, vids(6), out%s3w), 'write s3w')

    call nc_check(nf90_put_var(ncid, vids(7), transpose(out%s2u_t)), 'write s2u_t')
    call nc_check(nf90_put_var(ncid, vids(8), transpose(out%s2v_t)), 'write s2v_t')
    call nc_check(nf90_put_var(ncid, vids(9), transpose(out%s2w_t)), 'write s2w_t')
    call nc_check(nf90_put_var(ncid, vids(10), transpose(out%s3u_t)), 'write s3u_t')
    call nc_check(nf90_put_var(ncid, vids(11), transpose(out%s3v_t)), 'write s3v_t')
    call nc_check(nf90_put_var(ncid, vids(12), transpose(out%s3w_t)), 'write s3w_t')
  end subroutine write_direction_data

  subroutine write_sf_output(path, grid, rx_m, ry_m, out, &
       lt_out, write_linear_transport, &
       ist_out, write_interscale_transfer, &
       xsp_out, write_interspace_transfer, &
       pd_out, write_pdelta, &
       px_out, write_px)
    character(len=*), intent(in) :: path
    type(grid_info_t), intent(in) :: grid
    real(rk), intent(in) :: rx_m(:), ry_m(:)
    type(sf_outputs_t), intent(in) :: out
    type(linear_transport_outputs_t), intent(in), optional :: lt_out
    logical, intent(in), optional :: write_linear_transport
    type(interscale_transfer_outputs_t), intent(in), optional :: ist_out
    logical, intent(in), optional :: write_interscale_transfer
    type(interspace_transfer_outputs_t), intent(in), optional :: xsp_out
    logical, intent(in), optional :: write_interspace_transfer
    type(pdelta_outputs_t), intent(in), optional :: pd_out
    logical, intent(in), optional :: write_pdelta
    type(px_outputs_t), intent(in), optional :: px_out
    logical, intent(in), optional :: write_px

    integer :: ncid, dimid_rx, dimid_ry, dimid_x, dimid_y
    integer :: vid_x, vid_y, vid_rx, vid_ry
    integer :: vids_rx(12), vids_ry(12)
    integer :: v_lt_rx(4), v_lt_ry(4)
    integer :: v_ist_rx(3), v_ist_ry(3)
    integer :: v_xsp_rx(3), v_xsp_ry(3)
    integer :: v_pd_rx(7), v_pd_ry(7)
    integer :: v_px_rx(7), v_px_ry(7)
    logical :: do_lt, do_ist, do_xsp, do_pd, do_px

    call nc_check(nf90_create(trim(path), nf90_clobber, ncid), 'create output')
    call nc_check(nf90_def_dim(ncid, 'dimrx', size(rx_m), dimid_rx), 'def dimrx')
    call nc_check(nf90_def_dim(ncid, 'dimry', size(ry_m), dimid_ry), 'def dimry')
    call nc_check(nf90_def_dim(ncid, 'dimx', size(grid%x), dimid_x), 'def dimx')
    call nc_check(nf90_def_dim(ncid, 'dimy', size(grid%y), dimid_y), 'def dimy')

    call define_var_1d(ncid, 'x', dimid_x, vid_x)
    call define_var_1d(ncid, 'y', dimid_y, vid_y)
    call define_var_1d(ncid, 'rx', dimid_rx, vid_rx)
    call define_var_1d(ncid, 'ry', dimid_ry, vid_ry)
    call define_direction_vars(ncid, 'rx', 'y', dimid_rx, dimid_y, vids_rx)
    call define_direction_vars(ncid, 'ry', 'x', dimid_ry, dimid_x, vids_ry)

    do_lt = .false.
    if (present(write_linear_transport)) do_lt = write_linear_transport
    do_ist = .false.
    if (present(write_interscale_transfer)) do_ist = write_interscale_transfer
    do_xsp = .false.
    if (present(write_interspace_transfer)) do_xsp = write_interspace_transfer
    do_pd = .false.
    if (present(write_pdelta)) do_pd = write_pdelta
    do_px = .false.
    if (present(write_px)) do_px = write_px

    if (do_lt) then
      call define_var_1d(ncid, 'uXx_mean_d_deltau_2_dXx_rx_meanxt', dimid_rx, v_lt_rx(1))
      call define_var_1d(ncid, 'uXy_mean_d_deltau_2_dXy_rx_meanxt', dimid_rx, v_lt_rx(2))
      call define_var_1d(ncid, 'deltaux_mean_d_deltau_2_drx_rx_meanxt', dimid_rx, v_lt_rx(3))
      call define_var_1d(ncid, 'deltauy_mean_d_deltau_2_dry_rx_meanxt', dimid_rx, v_lt_rx(4))
      call define_var_1d(ncid, 'uXx_mean_d_deltau_2_dXx_ry_meanxt', dimid_ry, v_lt_ry(1))
      call define_var_1d(ncid, 'uXy_mean_d_deltau_2_dXy_ry_meanxt', dimid_ry, v_lt_ry(2))
      call define_var_1d(ncid, 'deltaux_mean_d_deltau_2_drx_ry_meanxt', dimid_ry, v_lt_ry(3))
      call define_var_1d(ncid, 'deltauy_mean_d_deltau_2_dry_ry_meanxt', dimid_ry, v_lt_ry(4))
    end if

    if (do_ist) then
      call define_var_1d(ncid, 'divrx_deltaux_deltau_2_rx_meanxt', dimid_rx, v_ist_rx(1))
      call define_var_1d(ncid, 'divry_deltauy_deltau_2_rx_meanxt', dimid_rx, v_ist_rx(2))
      call define_var_1d(ncid, 'divr_deltau_deltau_2_rx_meanxt', dimid_rx, v_ist_rx(3))
      call define_var_1d(ncid, 'divrx_deltaux_deltau_2_ry_meanxt', dimid_ry, v_ist_ry(1))
      call define_var_1d(ncid, 'divry_deltauy_deltau_2_ry_meanxt', dimid_ry, v_ist_ry(2))
      call define_var_1d(ncid, 'divr_deltau_deltau_2_ry_meanxt', dimid_ry, v_ist_ry(3))
    end if

    if (do_xsp) then
      call define_var_1d(ncid, 'divXx_uXx_deltau_2_rx_meanxt', dimid_rx, v_xsp_rx(1))
      call define_var_1d(ncid, 'divXy_uXy_deltau_2_rx_meanxt', dimid_rx, v_xsp_rx(2))
      call define_var_1d(ncid, 'divX_uX_deltau_2_rx_meanxt', dimid_rx, v_xsp_rx(3))
      call define_var_1d(ncid, 'divXx_uXx_deltau_2_ry_meanxt', dimid_ry, v_xsp_ry(1))
      call define_var_1d(ncid, 'divXy_uXy_deltau_2_ry_meanxt', dimid_ry, v_xsp_ry(2))
      call define_var_1d(ncid, 'divX_uX_deltau_2_ry_meanxt', dimid_ry, v_xsp_ry(3))
    end if

    if (do_pd) then
      call define_p_like_vars(ncid, 'Pdelta', 'rx_meanxt', dimid_rx, v_pd_rx)
      call define_p_like_vars(ncid, 'Pdelta', 'ry_meanxt', dimid_ry, v_pd_ry)
    end if

    if (do_px) then
      call define_p_like_vars(ncid, 'PX', 'rx_meanxt', dimid_rx, v_px_rx)
      call define_p_like_vars(ncid, 'PX', 'ry_meanxt', dimid_ry, v_px_ry)
    end if

    call nc_check(nf90_enddef(ncid), 'enddef output')

    call nc_check(nf90_put_var(ncid, vid_x, grid%x), 'write x')
    call nc_check(nf90_put_var(ncid, vid_y, grid%y), 'write y')
    call nc_check(nf90_put_var(ncid, vid_rx, rx_m), 'write rx')
    call nc_check(nf90_put_var(ncid, vid_ry, ry_m), 'write ry')
    call write_direction_data(ncid, vids_rx, out%rx)
    call write_direction_data(ncid, vids_ry, out%ry)

    if (do_lt .and. present(lt_out)) then
      call nc_check(nf90_put_var(ncid, v_lt_rx(1), lt_out%ux_dS2_dXx_rx), 'write lt rx1')
      call nc_check(nf90_put_var(ncid, v_lt_rx(2), lt_out%uy_dS2_dXy_rx), 'write lt rx2')
      call nc_check(nf90_put_var(ncid, v_lt_rx(3), lt_out%dux_dS2_drx_rx), 'write lt rx3')
      call nc_check(nf90_put_var(ncid, v_lt_rx(4), lt_out%duy_dS2_dry_rx), 'write lt rx4')
      call nc_check(nf90_put_var(ncid, v_lt_ry(1), lt_out%ux_dS2_dXx_ry), 'write lt ry1')
      call nc_check(nf90_put_var(ncid, v_lt_ry(2), lt_out%uy_dS2_dXy_ry), 'write lt ry2')
      call nc_check(nf90_put_var(ncid, v_lt_ry(3), lt_out%dux_dS2_drx_ry), 'write lt ry3')
      call nc_check(nf90_put_var(ncid, v_lt_ry(4), lt_out%duy_dS2_dry_ry), 'write lt ry4')
    end if

    if (do_ist .and. present(ist_out)) then
      call nc_check(nf90_put_var(ncid, v_ist_rx(1), ist_out%divrx_rx), 'write ist rx1')
      call nc_check(nf90_put_var(ncid, v_ist_rx(2), ist_out%divry_rx), 'write ist rx2')
      call nc_check(nf90_put_var(ncid, v_ist_rx(3), ist_out%divr_rx), 'write ist rx3')
      call nc_check(nf90_put_var(ncid, v_ist_ry(1), ist_out%divrx_ry), 'write ist ry1')
      call nc_check(nf90_put_var(ncid, v_ist_ry(2), ist_out%divry_ry), 'write ist ry2')
      call nc_check(nf90_put_var(ncid, v_ist_ry(3), ist_out%divr_ry), 'write ist ry3')
    end if

    if (do_xsp .and. present(xsp_out)) then
      call nc_check(nf90_put_var(ncid, v_xsp_rx(1), xsp_out%divXx_rx), 'write xsp rx1')
      call nc_check(nf90_put_var(ncid, v_xsp_rx(2), xsp_out%divXy_rx), 'write xsp rx2')
      call nc_check(nf90_put_var(ncid, v_xsp_rx(3), xsp_out%divX_rx), 'write xsp rx3')
      call nc_check(nf90_put_var(ncid, v_xsp_ry(1), xsp_out%divXx_ry), 'write xsp ry1')
      call nc_check(nf90_put_var(ncid, v_xsp_ry(2), xsp_out%divXy_ry), 'write xsp ry2')
      call nc_check(nf90_put_var(ncid, v_xsp_ry(3), xsp_out%divX_ry), 'write xsp ry3')
    end if

    if (do_pd .and. present(pd_out)) then
      call write_p_like_data(ncid, v_pd_rx, pd_out%p111_rx, pd_out%p121_rx, pd_out%p212_rx, &
           pd_out%p222_rx, pd_out%p313_rx, pd_out%p323_rx, pd_out%p_rx, 'pd rx')
      call write_p_like_data(ncid, v_pd_ry, pd_out%p111_ry, pd_out%p121_ry, pd_out%p212_ry, &
           pd_out%p222_ry, pd_out%p313_ry, pd_out%p323_ry, pd_out%p_ry, 'pd ry')
    end if

    if (do_px .and. present(px_out)) then
      call write_p_like_data(ncid, v_px_rx, px_out%p111_rx, px_out%p121_rx, px_out%p212_rx, &
           px_out%p222_rx, px_out%p313_rx, px_out%p323_rx, px_out%p_rx, 'px rx')
      call write_p_like_data(ncid, v_px_ry, px_out%p111_ry, px_out%p121_ry, px_out%p212_ry, &
           px_out%p222_ry, px_out%p313_ry, px_out%p323_ry, px_out%p_ry, 'px ry')
    end if

    call nc_check(nf90_close(ncid), 'close output')
  end subroutine write_sf_output

  subroutine define_p_like_vars(ncid, prefix, suffix, dimid_r, vids)
    integer, intent(in) :: ncid, dimid_r
    character(len=*), intent(in) :: prefix, suffix
    integer, intent(out) :: vids(7)

    call define_var_1d(ncid, trim(prefix)//'_111_'//trim(suffix), dimid_r, vids(1))
    call define_var_1d(ncid, trim(prefix)//'_121_'//trim(suffix), dimid_r, vids(2))
    call define_var_1d(ncid, trim(prefix)//'_212_'//trim(suffix), dimid_r, vids(3))
    call define_var_1d(ncid, trim(prefix)//'_222_'//trim(suffix), dimid_r, vids(4))
    call define_var_1d(ncid, trim(prefix)//'_313_'//trim(suffix), dimid_r, vids(5))
    call define_var_1d(ncid, trim(prefix)//'_323_'//trim(suffix), dimid_r, vids(6))
    call define_var_1d(ncid, trim(prefix)//'_'//trim(suffix), dimid_r, vids(7))
  end subroutine define_p_like_vars

  subroutine write_p_like_data(ncid, vids, p111, p121, p212, p222, p313, p323, ptotal, label)
    integer, intent(in) :: ncid, vids(7)
    real(rk), intent(in) :: p111(:), p121(:), p212(:), p222(:), p313(:), p323(:), ptotal(:)
    character(len=*), intent(in) :: label

    call nc_check(nf90_put_var(ncid, vids(1), p111), 'write '//trim(label)//' 111')
    call nc_check(nf90_put_var(ncid, vids(2), p121), 'write '//trim(label)//' 121')
    call nc_check(nf90_put_var(ncid, vids(3), p212), 'write '//trim(label)//' 212')
    call nc_check(nf90_put_var(ncid, vids(4), p222), 'write '//trim(label)//' 222')
    call nc_check(nf90_put_var(ncid, vids(5), p313), 'write '//trim(label)//' 313')
    call nc_check(nf90_put_var(ncid, vids(6), p323), 'write '//trim(label)//' 323')
    call nc_check(nf90_put_var(ncid, vids(7), ptotal), 'write '//trim(label)//' total')
  end subroutine write_p_like_data

end module mod_output_netcdf
