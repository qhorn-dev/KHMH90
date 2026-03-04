program khmh_s2
  use mod_kinds, only: i4, rk
  use mod_config, only: s2_config_t, read_config
  use mod_grid, only: grid_info_t, build_r_values
  use mod_stats_types, only: sf_accumulator_t, sf_outputs_t
  use mod_sf_terms, only: init_sf_accumulator, accumulate_sf_block, finalize_sf_accumulator
  use mod_linear_transport_terms, only: linear_transport_context_t, linear_transport_accumulator_t, &
       linear_transport_outputs_t, init_linear_transport_context, init_linear_transport_accumulator, &
       accumulate_linear_transport_block, finalize_linear_transport
  use mod_interscale_transfer_terms, only: interscale_transfer_accumulator_t, interscale_transfer_outputs_t, &
       init_interscale_transfer_accumulator, accumulate_interscale_transfer_block, finalize_interscale_transfer
  use mod_interspace_transfer_terms, only: interspace_transfer_accumulator_t, interspace_transfer_outputs_t, &
       init_interspace_transfer_accumulator, accumulate_interspace_transfer_block, finalize_interspace_transfer
  use mod_pdelta_terms, only: pdelta_accumulator_t, pdelta_outputs_t, &
       init_pdelta_accumulator, accumulate_pdelta_block, finalize_pdelta
  use mod_px_terms, only: px_accumulator_t, px_outputs_t, init_px_accumulator, accumulate_px_block, finalize_px
  use mod_netcdf_io, only: load_stats_and_grid, open_run_reader, close_run_reader, read_run_block, run_reader_t
  use mod_output_netcdf, only: write_sf_output
#ifdef _OPENMP
  use omp_lib, only: omp_set_num_threads
#endif
  implicit none

  type(s2_config_t) :: cfg
  type(grid_info_t) :: grid
  type(run_reader_t) :: rr
  type(sf_accumulator_t) :: acc
  type(sf_outputs_t) :: out
  type(linear_transport_context_t) :: lt_ctx
  type(linear_transport_accumulator_t) :: lt_acc
  type(linear_transport_outputs_t) :: lt_out
  type(interscale_transfer_accumulator_t) :: ist_acc
  type(interscale_transfer_outputs_t) :: ist_out
  type(interspace_transfer_accumulator_t) :: xsp_acc
  type(interspace_transfer_outputs_t) :: xsp_out
  type(pdelta_accumulator_t) :: pd_acc
  type(pdelta_outputs_t) :: pd_out
  type(px_accumulator_t) :: px_acc
  type(px_outputs_t) :: px_out

  integer(i4), allocatable :: rx_pts(:), ry_pts(:)
  real(rk), allocatable :: mean_u(:, :), mean_v(:, :), mean_w(:, :)
  real(rk), allocatable :: u(:, :, :), v(:, :, :), w(:, :, :)
  real(rk), allocatable :: rx_m(:), ry_m(:)

  integer(i4) :: irun, t0, nt_block
  integer(i4) :: run_start, run_end
  integer(i4) :: r_max_x_eff, r_max_y_eff
  character(len=512) :: run_path

  call read_config('config/s2_config.nml', cfg)
#ifdef _OPENMP
  if (cfg%nthreads > 0) call omp_set_num_threads(cfg%nthreads)
#endif
  if (.not. any((/ cfg%calc_sf, cfg%calc_linear_transport, cfg%calc_interscale_transfer, &
       cfg%calc_interspace_transfer, cfg%calc_pdelta, cfg%calc_px /))) then
    stop 'No term family enabled in config.'
  end if
  call load_stats_and_grid(cfg%stats_file, grid, mean_u, mean_v, mean_w)

  r_max_x_eff = cfg%r_max_x_pts
  if (r_max_x_eff < 0) r_max_x_eff = grid%nx - 1
  r_max_x_eff = min(r_max_x_eff, grid%nx - 1)

  r_max_y_eff = cfg%r_max_y_pts
  if (r_max_y_eff < 0) r_max_y_eff = grid%ny - 1
  r_max_y_eff = min(r_max_y_eff, grid%ny - 1)

  call build_r_values(cfg%r_min_x_pts, r_max_x_eff, cfg%r_step_x_pts, rx_pts)
  call build_r_values(cfg%r_min_y_pts, r_max_y_eff, cfg%r_step_y_pts, ry_pts)
  call init_sf_accumulator(acc, size(rx_pts), size(ry_pts), grid%ny, grid%nx)
  if (cfg%calc_linear_transport .or. cfg%calc_pdelta .or. cfg%calc_px) then
    call init_linear_transport_context(lt_ctx, grid, mean_u, mean_v, mean_w)
  end if
  if (cfg%calc_linear_transport) then
    call init_linear_transport_accumulator(lt_acc, size(rx_pts), size(ry_pts))
  end if
  if (cfg%calc_interscale_transfer) then
    call init_interscale_transfer_accumulator(ist_acc, size(rx_pts), size(ry_pts))
  end if
  if (cfg%calc_interspace_transfer) then
    call init_interspace_transfer_accumulator(xsp_acc, size(rx_pts), size(ry_pts))
  end if
  if (cfg%calc_pdelta) then
    call init_pdelta_accumulator(pd_acc, size(rx_pts), size(ry_pts))
  end if
  if (cfg%calc_px) then
    call init_px_accumulator(px_acc, size(rx_pts), size(ry_pts))
  end if

  run_start = cfg%run_start
  run_end = cfg%run_end

  do irun = run_start, run_end
    write(run_path, '(a,i2.2,a)') trim(cfg%input_prefix), irun, trim(cfg%input_suffix)
    call open_run_reader(trim(run_path), rr)

    t0 = 1
    do while (t0 <= rr%nt)
      nt_block = min(cfg%time_block_size, rr%nt - t0 + 1)
      call read_run_block(rr, t0, nt_block, u, v, w)

      u = u - spread(mean_u, dim=1, ncopies=nt_block)
      v = v - spread(mean_v, dim=1, ncopies=nt_block)
      w = w - spread(mean_w, dim=1, ncopies=nt_block)

      call accumulate_sf_block(u, v, w, rx_pts, ry_pts, acc)
      if (cfg%calc_linear_transport) then
        call accumulate_linear_transport_block(u, v, w, grid, rx_pts, ry_pts, lt_ctx, lt_acc)
      end if
      if (cfg%calc_interscale_transfer) then
        call accumulate_interscale_transfer_block(u, v, w, grid, rx_pts, ry_pts, ist_acc)
      end if
      if (cfg%calc_interspace_transfer) then
        call accumulate_interspace_transfer_block(u, v, w, grid, rx_pts, ry_pts, xsp_acc)
      end if
      if (cfg%calc_pdelta) then
        call accumulate_pdelta_block(u, v, w, grid, rx_pts, ry_pts, lt_ctx, pd_acc)
      end if
      if (cfg%calc_px) then
        call accumulate_px_block(u, v, w, grid, rx_pts, ry_pts, lt_ctx, px_acc)
      end if

      deallocate(u, v, w)
      t0 = t0 + nt_block
    end do

    call close_run_reader(rr)
  end do

  call finalize_sf_accumulator(acc, out)
  if (cfg%calc_linear_transport) then
    call finalize_linear_transport(lt_acc, lt_out)
  end if
  if (cfg%calc_interscale_transfer) then
    call finalize_interscale_transfer(ist_acc, ist_out)
  end if
  if (cfg%calc_interspace_transfer) then
    call finalize_interspace_transfer(xsp_acc, xsp_out)
  end if
  if (cfg%calc_pdelta) then
    call finalize_pdelta(pd_acc, pd_out)
  end if
  if (cfg%calc_px) then
    call finalize_px(px_acc, px_out)
  end if

  allocate(rx_m(size(rx_pts)), ry_m(size(ry_pts)))
  rx_m = real(rx_pts, rk) * grid%dx
  ry_m = real(ry_pts, rk) * grid%dy

  call write_sf_output(cfg%output_file, grid, rx_m, ry_m, out, lt_out, cfg%calc_linear_transport, &
       ist_out, cfg%calc_interscale_transfer, xsp_out, cfg%calc_interspace_transfer, &
       pd_out, cfg%calc_pdelta, px_out, cfg%calc_px)

  write(*,'(a)') 'Done.'
end program khmh_s2
