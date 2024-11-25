.. toctree::
   :maxdepth: 2
   :caption: Contents:

========
Examples
========

Here we include three examples of how to compile and run three basic thornado
applications:

1. A sine wave advection problem with GRHD using only thornado

2. A sine wave advection problem with GRHD using thornado and AMReX

3. A <problem> with :math:`\mathcal{O}\left(v\right)` neutrino transport
   using only thornado

To use thornado, a user will need some dependencies:

HDF5 (from your package manager, or download the
`source code <https://www.hdfgroup.org/download-hdf5/source-code/>`_,
or get it via `github <https://github.com/HDFGroup/hdf5>`_);

LAPACK (from your package manager, or download the
`source code <https://www.netlib.org/lapack/>`_,
or get it via `github <https://github.com/Reference-LAPACK/lapack>`_);

BLAS (from your package manager, or download the
`source code <https://www.netlib.org/blas/>`_)
(BLAS is included in the LAPACK github repository).

Additionally, if the user wants to use AMReX, they will need to clone their
`github <https://github.com/AMReX-Codes/amrex>`_ repository.

------------------
Running with AMReX
------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Evolve GRHD sine wave advection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0. Download dependencies: HDF5, LAPACK, BLAS, open-mpi (or mpich)

    * On Mac (assuming you are using `homebrew <https://brew.sh/>`_)::

        $ brew install hdf5 lapack openblas open-mpi

    * On Ubuntu::

        $ sudo apt-get install libhdf5-dev liblapack-dev libblas-dev libopenmpi-dev

1. Clone thornado::

    $ git clone https://github.com/endeve/thornado

2. Export environment variable::

    $ export THORNADO_DIR=/path/to/thornado

3. Clone amrex::

    $ git clone https://github.com/AMReX-Codes/amrex

4. Export environment variable::

    $ export AMREX_DIR=/path/to/amrex

5. Navigate to appropriate directory::

    $ cd ${THORNADO_DIR}/SandBox/AMReX/dgExperiments_Euler_Relativistic_IDEAL

6. Compile code::

    $ make DIM=1 DEBUG=FALSE -jN
where ``N`` is the desired number of threads (e.g., 20)

7. Execute::

    $ mpiexec -n N ./main1d.gnu.MPI.ex Advection1D_Gaussian_test.inputs > out
where ``N`` is the desired number of MPI processes.

-------------------------
Example AMReX inputs file
-------------------------

.. literalinclude:: grhd.Example.inputs
