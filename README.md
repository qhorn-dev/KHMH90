# KHMH

Fortran code to compute structure functions and several KHMH budget terms from experimental NetCDF files.

## Prerequisites

- `cmake`
- a Fortran compiler (e.g. gfortran)
- the `netcdf-fortran` library
- optional: OpenMP

## Compilation

Standard build with NetCDF:

```bash
cmake -S . -B build -DKHMH_USE_NETCDF=ON
cmake --build build -j
```

The generated executable is:

```bash
./build/khmh_s2
```

## Running

By default, the program reads:

- the configuration file [s2_config.nml](config/s2_config.nml)
- the input NetCDF files specified in that configuration

Basic run:

```bash
./build/khmh_s2
```

If OpenMP is available, the number of threads can be set directly in the configuration with `nthreads`.

Example:

```fortran
nthreads = 8
```

If `nthreads = 0`, OpenMP uses its default behavior which is usueally set by the environment variable:

```bash
export OMP_NUM_THREADS=8
./build/khmh_s2
```

## Configuration

Configuration is handled in [s2_config.nml](config/s2_config.nml).

Example:

```fortran
&sf
  input_prefix = 'inputs/2bars_P1_run'
  input_suffix = '.nc'
  stats_file = 'inputs/2bars_P1_stats.nc'
  output_file = 's2_output.nc'

  run_start = 1
  run_end = 1

  r_min_x_pts = 1
  r_max_x_pts = -1
  r_step_x_pts = 1

  r_min_y_pts = 1
  r_max_y_pts = -1
  r_step_y_pts = 1

  time_block_size = 16
  nthreads = 0

  all_terms = .false.
  calc_sf = .true.
  calc_linear_transport = .true.
  calc_interscale_transfer = .true.
  calc_interspace_transfer = .false.
  calc_pdelta = .false.
  calc_px = .false.
/
```

### Input and output parameters

- `input_prefix`: prefix of the run files
- `input_suffix`: suffix of the run files
- `stats_file`: NetCDF file containing the grid and mean fields
- `output_file`: output NetCDF file

Run file names are built as:

```text
input_prefix + NN + input_suffix
```

with `NN` on 2 digits, for example `run01`, `run02`, etc.

### Run range

- `run_start`
- `run_end`

Example:

```fortran
run_start = 1
run_end = 4
```

This processes runs `01` to `04`.

### Separations

- `r_min_x_pts`, `r_max_x_pts`, `r_step_x_pts`: separations in `x`, in grid points
- `r_min_y_pts`, `r_max_y_pts`, `r_step_y_pts`: separations in `y`, in grid points

Convention:

- `r_max_x_pts = -1` means: use all possible separations up to `nx-1`
- `r_max_y_pts = -1` means: use all possible separations up to `ny-1`

### Time streaming

- `time_block_size`: number of time instants read per block

A larger value means:

- fewer I/O calls
- more memory usage

That's why this value has to be adapted to the RAM you have. You can for example proceed with `time_block_size=100` and lower it if you get some out of memory issue.

### Threads

- `nthreads`: number of OpenMP threads

Convention:

- `nthreads = 0`: default OpenMP behavior
- `nthreads > 0`: the code explicitly sets this number of threads at startup

### Term families

- `all_terms`
- `calc_sf`
- `calc_linear_transport`
- `calc_interscale_transfer`
- `calc_interspace_transfer`
- `calc_pdelta`
- `calc_px`

Convention:

- if `all_terms = .true.`, all families are enabled
- otherwise, each family is controlled individually

## Mathematical meaning of the switches

The code works with velocity fluctuations:

\[
u' = u - \overline{u}, \qquad
v' = v - \overline{v}, \qquad
w' = w - \overline{w}
\]

and with velocity increments between two points separated by `r_x` or `r_y`.

For a separation in `x`:

\[
\delta u_x = u'(x+r_x,y,t)-u'(x,y,t)
\]
\[
\delta u_y = v'(x+r_x,y,t)-v'(x,y,t)
\]
\[
\delta u_z = w'(x+r_x,y,t)-w'(x,y,t)
\]

For a separation in `y`:

\[
\delta u_x = u'(x,y+r_y,t)-u'(x,y,t)
\]
\[
\delta u_y = v'(x,y+r_y,t)-v'(x,y,t)
\]
\[
\delta u_z = w'(x,y+r_y,t)-w'(x,y,t)
\]

The code produces:

- outputs averaged over the separation direction and time
- outputs not averaged in the transverse direction

### `calc_sf`

Second- and third-order structure functions.

For a separation in `x`, for example:

\[
S_{2,x}^{(u)}(r_x,y)=\left\langle (\delta u_x)^2 \right\rangle_{x,t}
\]
\[
S_{3,x}^{(u)}(r_x,y)=\left\langle (\delta u_x)^3 \right\rangle_{x,t}
\]

and similarly for the `v` and `w` components.

The code also outputs the versions averaged over the transverse coordinate:

\[
S_{2,x}^{(u)}(r_x)=\left\langle S_{2,x}^{(u)}(r_x,y) \right\rangle_y
\]
\[
S_{3,x}^{(u)}(r_x)=\left\langle S_{3,x}^{(u)}(r_x,y) \right\rangle_y
\]

Similarly for `r_y`:

\[
S_{2,y}^{(u)}(r_y,x)=\left\langle (\delta u_x)^2 \right\rangle_{y,t}
\]
\[
S_{3,y}^{(u)}(r_y,x)=\left\langle (\delta u_x)^3 \right\rangle_{y,t}
\]

### `calc_linear_transport`

Linear terms of the form:

\[
\overline{\mathbf{U}}_X \cdot \nabla_X |\delta \mathbf{u}|^2
\qquad \text{and} \qquad
\overline{\delta \mathbf{U}} \cdot \nabla_r |\delta \mathbf{u}|^2
\]

Associated output variables:

- `uXx_mean_d_deltau_2_dXx_*`
- `uXy_mean_d_deltau_2_dXy_*`
- `deltaux_mean_d_deltau_2_drx_*`
- `deltauy_mean_d_deltau_2_dry_*`

### `calc_interscale_transfer`

Interscale transfer in separation space:

\[
\nabla_r \cdot \left( \delta \mathbf{u} \, |\delta \mathbf{u}|^2 \right)
\]

The code outputs:

- `divrx_deltaux_deltau_2_*`
- `divry_deltauy_deltau_2_*`
- `divr_deltau_deltau_2_*`

### `calc_interspace_transfer`

Transport in physical center space:

\[
\nabla_X \cdot \left( \mathbf{U}_X \, |\delta \mathbf{u}|^2 \right)
\]

The code outputs:

- `divXx_uXx_deltau_2_*`
- `divXy_uXy_deltau_2_*`
- `divX_uX_deltau_2_*`

### `calc_pdelta`

Production term associated with mean gradients in separation space:

\[
P_{\delta}
=
- \delta u_i \, \delta u_j \,
\frac{\partial \overline{\delta U_i}}{\partial r_j}
\]

Resolved subterms:

- `Pdelta_111`
- `Pdelta_121`
- `Pdelta_212`
- `Pdelta_222`
- `Pdelta_313`
- `Pdelta_323`

and their sum:

\[
P_{\delta}
=
P_{\delta,111}+P_{\delta,121}+P_{\delta,212}+P_{\delta,222}+P_{\delta,313}+P_{\delta,323}
\]

### `calc_px`

Production term associated with gradients in center coordinates:

\[
P_X
=
- \delta u_i \, U_{X,j} \,
\frac{\partial \overline{\delta U_i}}{\partial X_j}
\]

Resolved subterms:

- `PX_111`
- `PX_121`
- `PX_212`
- `PX_222`
- `PX_313`
- `PX_323`

and their sum:

\[
P_X
=
P_{X,111}+P_{X,121}+P_{X,212}+P_{X,222}+P_{X,313}+P_{X,323}
\]

## Outputs

The code writes a NetCDF file containing:

- coordinates `x`, `y`, `rx`, `ry`
- second- and third-order structure functions
- averaged versions:
  - `SF(rx)`
  - `SF(ry)`
- transverse non-averaged versions:
  - `SF(rx,y)`
  - `SF(ry,x)`
- optional term families enabled in the namelist:
  - `linear_transport`
  - `interscale_transfer`
  - `interspace_transfer`
  - `Pdelta`
  - `PX`

