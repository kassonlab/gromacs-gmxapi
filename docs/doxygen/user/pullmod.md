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
* Additional input: environment variable provides a filename that is opened and
  parsed by additional code in the extension.

Old implementation
------------------

- do_md()
    - do_force()
      - do_force_cutsVERLET()
        - pull_potential_wrapper()
          - pull_potential()
            - do_pull_pot_coord()
              - calc_pull_coord_force()
                - getRouxForce()
                - getRouxPotential()

```
static void calc_pull_coord_force(pull_coord_work_t *pcrd, int coord_ind,
                              double dev, real lambda,
                              real *V, tensor vir, real *dVdl)
{   //...
    switch (pcrd->params.eType)
    {   //...
        case epullROUX:
            dev = pcrd->value;
            /* We want to pass the value of the pull coordinate, not the
             * difference between the pull coordinate and pull-init (which
             * would be dev = pcrd->value - pcrd->value_ref)
             */
            pcrd->f_scal = getRouxForce(dev, coord_ind, k);
            *V          += getRouxPotential(dev, coord_ind, k);
            *dVdl       += *V/k * dkdl;
            break;
    }
    //...
}
```

Current pull code implementation
--------------------------------

Relevant types:

* t_pull_coord
* pull_t
* pull_coord_work_t

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

- `gmx_mdrun()`
  - `Mdrunner:mdrunner()`
    - `init_pull()` -> `pull_t`

- `do_md()` accesses the `pull_t` member of the input record to call...
  - `do_force()`
    - `do_force_cutsVERLET()`
      - `pull_potential_wrapper()`
        - `pull_potential()`
          - `do_pull_pot_coord()`
            - `calc_pull_coord_scalar_force_and_potential(pcrd, dev, lambda, V, dVdl);`
            - `calc_pull_coord_vector_force(pcrd);`


Strategy
--------

Reimplement `init_pull()` and refactor `pull_t`. Allow client code to provide alternatives.

Refactor `do_pull_pot_coord()` to use member functions of `pull_t` argument.

Where control flow is currently guided by enum values, wrap the function to dynamically dispatch to functor or legacy free functions.

A custom `PullPotential` class may (re)implement `calc_pull_coord_scalar_force_and_potential()` or `calc_pull_coord_vector_force();` or some subset of their internals.

Note that `pull_potential()` is the lowest public function in this hierarchy and is passed a pull_t from the inputrec.

Usage
-----

Instantiate a modified `do_md()` integrator with a builder for the plugin to replace `init_pull()`

Extension code inherits from base `do_pull_pot_coord()` provider that `pull_potential()` will call. Override `calc_...()` method (or submethod).


Python client code:

```{.py}
import numpy
import csv
import roux
import gmx

with open('rawdata.csv', 'r') as datafile:
    histogram = numpy.asarray(csv.reader(datafile), dtype=float)

puller = roux.PullPotential()
puller.histogram = histogram
puller.sigma = 1.0
puller.k = 30
puller.nbins = 50

system = gmx.system.from_tpr('topol.tpr')
# Note: input record does not indicate the custom code
system.md.addPotential(puller)

# In the simple case, let nranks == num_ensemble_members
with gmx.context.MpiEnsemble(system) as session:
    # get an mpi4py COMM_WORLD handle
    comm = session.communicator
    size = comm.Get_size()
    rank = comm.Get_rank()
    for iteration in range(10):
        session.run()

        # gather and accumulate statistics
        recvbuf = None
        if rank == 0:
            recvbuf_shape = [size] + puller.data.shape
            recvbuf = numpy.empty(recvbuf_shape, dtype=float)
        comm.Gather(puller.data, recvbuf, root=0)

        # perform analysis
        new_histogram = roux.analyze(recvbuf)

        # broadcast updated histogram
        comm.Bcast(new_histogram, root=0)
        puller.histogram = new_histogram
        if rank == 0:
            with open('rawdata.csv', 'w') as datafile:
                csv.writer(datafile).writerows(puller.histogram)

```

sample C++ client code

```{.cpp}

RouxPuller::Builder rouxSetup;
rouxSetup.addHistogram(arraydata);
rouxSetup.setDimensions(...);
std::shared_ptr<RouxPuller> myRoux = rouxSetup.build();

gmx::MdrunnerBuilder simulation;
simulation.addPullPotential(myRoux);
gmx::Mdrunner session = myContext(simulation);

for (auto i=0; i < n_iter; i++)
{
    session->run();
    myRoux->update();
}
```

Implementing a PullPotential derived class.

```{.cpp}

class RouxPuller : public PullPotential
{
public:
    /// Something we can use to set parameters from external code.
    class Builder;
    
    /// Provide interface for runner to set parameters
    struct pull_t* init_pull(FILE *fplog,
                             const pull_params_t *pull_params,
                             const t_inputrec *ir,
                             int nfile,
                             const t_filenm fnm[],
                             const gmx_mtop_t *mtop,
                             t_commrec *cr,
                             const gmx_output_env_t *oenv,
                             real lambda,
                             gmx_bool bOutFile,
                             unsigned long Flags) override;
    
    // No need to override
    /* static void do_pull_pot_coord(struct pull_t *pull, int coord_ind, t_pbc *pbc,
                                     double t, real lambda,
                                     real *V, tensor vir, real *dVdl) override;
     */
    
    // No need to override
    // calc_pull_coord_vector_force(pcrd) override;

    /// Implement potential and force
    gmxapi::Status calc_pull_coord_scalar_force_and_potential(struct pull_t *pull,
                                                              int coord_ind,
                                                              const t_pbc *pbc,
                                                              double t,
                                                              real lambda,
                                                              real *V,
                                                              real *dVdl
                                                              ) override;
    
    // No need to override
    // void finish_pull() override;
    
    /// Some sort of external interface
    gmxapi::Status update() override;
    
    // custom member functions.
    // ...
    
private:
    /// custom Potential parameters
    double sigma_;
    double k_;
    /// Histogram data for reference
    std::vector< std::vector<double> > histogram_;
};

struct pull_t* RouxPuller::init_pull(FILE *fplog,
                             const pull_params_t *pull_params,
                             const t_inputrec *ir,
                             int nfile,
                             const t_filenm fnm[],
                             const gmx_mtop_t *mtop,
                             t_commrec *cr,
                             const gmx_output_env_t *oenv,
                             real lambda,
                             gmx_bool bOutFile,
                             unsigned long Flags)
{
    assert(!histogram_.empty());
    this->setCommunicator(cr);
    this->setInputRecord(ir);
    return this;
}

void RouxPuller::calc_pull_coord_scalar_force_and_potential(int coord_ind,
                                                              double t,
                                                              real lambda,
                                                              real *V,
                                                              real *dVdl
                                                              )
{
    pull_coord_work_t *pcrd = this->coord[coord_ind];
    // Calculate distances and set V, dVdl
    dev = pcrd->value;
    /* We want to pass the value of the pull coordinate, not the
     * difference between the pull coordinate and pull-init (which
     * would be dev = pcrd->value - pcrd->value_ref)
     */
    pcrd->f_scal = getRouxForce(dev, coord_ind, k);
    *V          += getRouxPotential(dev, coord_ind, k);
    *dVdl       += *V/k * dkdl;
}

PYBIND11_MODULE(roux_, m) {
    m.doc() = somedocstring;

    pybind11::class_< RouxPuller, PullPotential, std::shared_ptr<RouxPuller> > roux_plugin(m, "PullPotential", "Implements Roux biasing potential.");
    roux_plugin.def_property("data", &RouxPuller::getData);
    // The rest is inherited from PullPotential
    
    pybind11::class_< RouxBuilder > roux_builder(m, "Builder", "Set up potential");
    roux_builder.def(pybind11::init());
    roux_builder.def_property("histogram", &RouxBuilder::addHistogram);
    roux_builder.def_property("dimensions", &RouxBuilder::setDimensions);
    roux_builder.def_property("sigma", &RouxBuilder::setSigma);
    roux_builder.def_property("k", &RouxBuilder::setK);
    roux_builder.def_property("nbins", &RouxBuilder::setNBins);
    roux_builder.def("build", &RouxBuilder::build);
}

```

Python wrapper

```{.py}

# roux.py
import gmx
import roux_.PullPotential
import roux_.Builder

class PullPotential(gmx.PullPotential):
    ...
    def bind(self, system):
        builder = roux_.Builder()
        builder.histogram = self.histogram
        builder.sigma = self.sigma
        builder.k = self.k
        builder.nbins = self.nbins
        self.api_object = builder.build();
        ...
     ...
     @property
     def histogram(self):
        ...
     @property
     def sigma(self):
        ...
     ...

def analyze(data):
    ...
    return histogram

```

...
