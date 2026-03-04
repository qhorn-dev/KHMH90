module mod_config
  use mod_kinds, only: i4
  implicit none
  private

  type, public :: s2_config_t
    character(len=256) :: input_prefix = 'inputs/2bars_P1_run'
    character(len=256) :: input_suffix = '.nc'
    character(len=256) :: stats_file = 'inputs/2bars_P1_stats.nc'
    character(len=256) :: output_file = 's2_output.nc'

    integer(i4) :: run_start = 1
    integer(i4) :: run_end = 1

    integer(i4) :: r_min_x_pts = 1
    integer(i4) :: r_max_x_pts = -1
    integer(i4) :: r_step_x_pts = 1

    integer(i4) :: r_min_y_pts = 1
    integer(i4) :: r_max_y_pts = -1
    integer(i4) :: r_step_y_pts = 1

    integer(i4) :: time_block_size = 16
    integer(i4) :: nthreads = 0

    logical :: all_terms = .false.
    logical :: calc_sf = .true.
    logical :: calc_linear_transport = .false.
    logical :: calc_interscale_transfer = .false.
    logical :: calc_interspace_transfer = .false.
    logical :: calc_pdelta = .false.
    logical :: calc_px = .false.
  end type s2_config_t

  public :: read_config

contains

  subroutine read_config(path, cfg)
    character(len=*), intent(in) :: path
    type(s2_config_t), intent(inout) :: cfg

    character(len=256) :: input_prefix, input_suffix, stats_file, output_file
    integer(i4) :: run_start, run_end
    integer(i4) :: r_min_x_pts, r_max_x_pts, r_step_x_pts
    integer(i4) :: r_min_y_pts, r_max_y_pts, r_step_y_pts
    integer(i4) :: time_block_size, nthreads
    logical :: all_terms
    logical :: calc_sf, calc_linear_transport, calc_interscale_transfer
    logical :: calc_interspace_transfer, calc_pdelta, calc_px
    integer :: u, ios

    namelist /sf/ input_prefix, input_suffix, stats_file, output_file, &
                  run_start, run_end, &
                  r_min_x_pts, r_max_x_pts, r_step_x_pts, &
                  r_min_y_pts, r_max_y_pts, r_step_y_pts, &
                  time_block_size, nthreads, &
                  all_terms, calc_sf, calc_linear_transport, calc_interscale_transfer, &
                  calc_interspace_transfer, calc_pdelta, calc_px

    input_prefix = cfg%input_prefix
    input_suffix = cfg%input_suffix
    stats_file = cfg%stats_file
    output_file = cfg%output_file
    run_start = cfg%run_start
    run_end = cfg%run_end
    r_min_x_pts = cfg%r_min_x_pts
    r_max_x_pts = cfg%r_max_x_pts
    r_step_x_pts = cfg%r_step_x_pts
    r_min_y_pts = cfg%r_min_y_pts
    r_max_y_pts = cfg%r_max_y_pts
    r_step_y_pts = cfg%r_step_y_pts
    time_block_size = cfg%time_block_size
    nthreads = cfg%nthreads
    all_terms = cfg%all_terms
    calc_sf = cfg%calc_sf
    calc_linear_transport = cfg%calc_linear_transport
    calc_interscale_transfer = cfg%calc_interscale_transfer
    calc_interspace_transfer = cfg%calc_interspace_transfer
    calc_pdelta = cfg%calc_pdelta
    calc_px = cfg%calc_px

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') 'Config file not found. Using built-in defaults: '//trim(path)
      return
    end if

    read(u, nml=sf, iostat=ios)
    close(u)
    if (ios /= 0) then
      write(*,'(a)') 'Failed reading namelist /sf/. Using defaults for missing fields.'
    end if

    cfg%input_prefix = input_prefix
    cfg%input_suffix = input_suffix
    cfg%stats_file = stats_file
    cfg%output_file = output_file
    cfg%run_start = run_start
    cfg%run_end = run_end
    cfg%r_min_x_pts = r_min_x_pts
    cfg%r_max_x_pts = r_max_x_pts
    cfg%r_step_x_pts = r_step_x_pts
    cfg%r_min_y_pts = r_min_y_pts
    cfg%r_max_y_pts = r_max_y_pts
    cfg%r_step_y_pts = r_step_y_pts
    cfg%time_block_size = time_block_size
    cfg%nthreads = nthreads

    if (all_terms) then
      calc_sf = .true.
      calc_linear_transport = .true.
      calc_interscale_transfer = .true.
      calc_interspace_transfer = .true.
      calc_pdelta = .true.
      calc_px = .true.
    end if

    cfg%all_terms = all_terms
    cfg%calc_sf = calc_sf
    cfg%calc_linear_transport = calc_linear_transport
    cfg%calc_interscale_transfer = calc_interscale_transfer
    cfg%calc_interspace_transfer = calc_interspace_transfer
    cfg%calc_pdelta = calc_pdelta
    cfg%calc_px = calc_px
  end subroutine read_config

end module mod_config
