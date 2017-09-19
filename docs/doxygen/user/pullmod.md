Extensibility of pull code {#page_pullcodeextension}
==========================

Historically, bias potentials have been applied by patching \Gromacs at the pull code.

This approach introduces an immediate difficulty in staying up to date with other
developments in \Gromacs and leaves many problems for the researcher, such as how to
get user input into the additional code and how to transfer data in and out of the
extension code.

Case study
----------

2015 \Gromacs fork adding an ensemble biasing scheme after Roux, dx.doi.org/10.1021/jp3110369

Commit ac9ce76

Summary of changes:

    +11  -0  M docs/user-guide/mdp-options.rst
    +4   -0  M src/gromacs/fileio/tpxio.cpp
    +45  -3  M src/gromacs/gmxana/gmx_mindist.cpp
    +2   -2  M src/gromacs/gmxlib/names.cpp
    +7   -1  M src/gromacs/gmxpreprocess/readpull.cpp
    +3   -2  M src/gromacs/legacyheaders/types/enums.h
    +298 -0  A src/gromacs/pulling/csvparser.cpp
    +48  -0  A src/gromacs/pulling/csvparser.h
    +75  -27 M src/gromacs/pulling/pull.cpp
    +1   -0  M src/gromacs/pulling/pullutil.cpp
    +113 -0  A src/gromacs/pulling/roux.cpp
    +9   -0  A src/gromacs/pulling/roux.h
    +1   -1  M src/programs/mdrun/md.cpp

First, note that over 300 lines are added for processing comma-separated-value data for the extension code.
This sort of dependency certainly should not need to be satisfied in the `libgromacs` code.

Globals modified

`gmxlib/names.cpp`

* add "roux" to `const char *epull_names[]`
* add "distance-reference" to `const char *epullg_names[]`

`/legacyheaders/types/enums.h`

* add `epullROUX`
* add `epullgDISTREF`

Functions modified

`fileio/tpxio.cpp`

* `static void do_pull_coord()`
  * call `gmx_fio_do_int(fio, pcrd->group[2])` for `pcrd->eGeom == epullgDIRRELATIVE`

`gmxana/gmx_mindist.cpp`

* `void dist_plot()`
  * add `int rectmat` call parameter
  * handle rectangular distance matrix
* `int gmx_mindist()`
  * additional data member and CLI argument
  * handle new call signature for `dist_plot()`

`gmxpreprocess/readpull.cpp`

* `static void init_pull_coord()`
  * extend scope of calculation to apply for `pcrd->eGeom == epullgDISTREF`
* `char **read_pullparams()`
  * set `ngroup = 3` when `pcrd->eGeom == epullgDISTREF`
* `void make_pull_coords()`
  * Add a line of status output when `pcrd->eGeom == epullgDISTREF`

`pulling/pull.cpp`

* `static void low_get_pull_coord_dr()`
  * add `dvec xref2` to call signature
  * change how `dr[m]` and `dr2` are set when `pcrd->params.eGeom == epullgDISTREF`
  * extend scope of calculation to also apply when `pcrd->params.eGeom == epullgDISTREF`
* `static void get_pull_coord_dr()`
  * add a small data member
  * handle new call signature for `low_get_pull_coord_dr()`
* `static void get_pull_coord_distance()`
  * extend scope of calculation to also apply when `pcrd->params.eGeom == epullgDISTREF`
* `static double get_pull_coord_deviation()`
  * extend scope of calculations to an additional case `epullgDISTREF`
* `static void do_constraint()`
  * add a small data member
  * adapt to new call signature of `low_get_pull_coord_dr()`
  * extend scope of calculations to an additional case `epullgDISTREF`
* `static void calc_pull_coord_force()`
  * add parameter to call signature
  * extend scope of calculations to apply to `pcrd->params.eGeom == epullgDISTREF`
  * add handling for case `epullROUX`
    * set `dev`, `pcrd->f_scal`, `*V`, and `*dVdl`
    * call `getRouxForce(dev, coord_ind, k)` and `getRouxPotential(dev, coord_ind, k)`
* `void set_pull_coord_reference_value()`
  * adapt to new call signature of `calc_pull_coord_force()`
* `static void do_pull_pot_coord()`
  * adapt to new call signature of `calc_pull_coord_force()`
* `struct pull_t * init_pull()`
  * extend scope of calculations to an additional case `epullgDISTREF`
  
Functions added

* `double getRouxForce(double dev, int coord_ind, double K)`
* `double getRouxPotential(double dev, int coord_ind, double K)`

Other

* `md.cpp` uses `step_rel` instead of `step` to determine whether this is a neighbor search step.


Current pull code implementation
--------------------------------

Relevant types:

* \ref t_pull_coord
* \ref pull_t
* \ref pull_coord_work_t

Dependents

* `t_mdebin` tracks pull as an energy data provider

Hooks:

* `register_external_pull_potential()` only checks consistency. It doesn't actually enable pulling implementations.
  That is hard-coded and must be supported by `init_pull()`, `pull_potential()`, `do_potential()`, `apply_external_pull_coord_force()`
  
Code flow:

`Mdrunner::mdrunner()`

conditionally set `inputrec->pull_work = init_pull(...)`

`finish_pull(inputrec->pull_work)` after call to integrator

`do_md()`

call `pull_print_output(ir->pull_work, step, t)` in simulation loop.

`do_md()` calls `do_force()` calls `do_force_cutsVERLET()` calls `pull_potential_wrapper()` calls `pull_potential()`


Strategy
--------

Reimplement `init_pull()` and refactor `pull_t`. Allow client code to provide alternatives.
