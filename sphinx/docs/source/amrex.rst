==================
Running with AMReX
==================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

------------------------------------------
Evolve GRHD sine wave advection with AMReX
------------------------------------------

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

-------------------
Example inputs file
-------------------

.. literalinclude:: grhd.Example.inputs
