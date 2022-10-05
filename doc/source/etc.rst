Support for Linking with External Codes
=============================================

Using the ChIMES Calculator with LAMMPS
******************************************

We are currently working toward ChIMES calculator implementation in `LAMMPS <https://lammps.sandia.gov>`_ as a USER package. In the interim, the following provides a guide to implementing the ChIMES calculator as a LAMMPS pairstyle.



Quick start
^^^^^^^^^^^^^^^^

Provided a system with a C++11-compatible compiler and an MPI compatible compiler are available, LAMMPS can be downloaded, installed, linked to ChIMES, and compiled all at once by navigating to ``etc/lmp``, adding Intel compilers to your path  and executing ``./install.sh``. Once complete, the installation can be tested by navigating to ``etc/lmp/tests`` and running the example via ``../exe/lmp_mpi_chimes -i in.lammps``. 

As with installation of the ChIMES Calculator itself, if you are on a HPC using module files, you may need to load them first. Module files are already configured for a handful of HPC - inspect the contents of modfiles to see if
yours is listed. If it is (e.g., LLNL-LC.mod), execute ``export hosttype=LLNL-LC; ./install.sh`` to install. Otherwise, load the appropriate modules by hand before running the
install script.

Note that Intel oneapi compilers (which are now free) can be used to properly configure your enviroment for all Intel capabilities (e.g., icc, mpiicpc, mkl, etc.) - simply locate and execute the setvars.sh script within your Intel installation.


-----

Detailed Compilation Overview
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. Note::

    This example assumes users have downloaded the 29 Oct 2020 release of LAMMPS (stable version as of 10/29/20), which can be downloaded `here <https://lammps.sandia.gov/download.html>`_. 
    
    

To integrate the ChIMES calculator in LAMMPS, locate the following files, and place them in the following destination among the LAMMPS source code:

========================    ================    ================    ==============
File                        Location            Destination         Description
========================    ================    ================    ==============
``chimesFF.{h,cpp}``        ``chimesFF/src``    ``src/MANYBODY``    ChIMES calculator files
``pair_chimes.{h,cpp}``     ``etc/lmp/src``     ``src/MANYBODY``    ChIMES pair_style definition files
``pair.{h,cpp}``            ``etc/lmp/etc``     ``src``             Updated LAMMPS pair files (new ev_tally definition added)
``Makefile.mpi_chimes``     ``etc/lmp/etc``     ``src/MAKE``        Makefile for compiling with ChIMES support
========================    ================    ================    ==============


Following, compile from the base LAMMPS directory with:

.. code-block:: shell

    make yes-manybody
    make mpi_chimes

Note that a successful compilation should produce an executable named ``lmp_mpi_chimes``.

.. Tip::

        If you are using an intel compiler, either delete the ``pair_list.*`` files that appear in the src folder following the ``make yes-manybody`` command, or add ``-restrict`` to ``CCFLAGS`` in ``MAKE/Makefile.mpi_chimes``. Note that the presently provided ``Makefile.mpi_chimes`` utilizes the latter approach.


Running
^^^^^^^^^^^^^^^^

To run a simulation using ChIMES parameters, a block like the following is needed in the main LAMMPS input file (i.e. ``in.lammps``):

.. code-block:: text

    pair_style	chimesFF
    pair_coeff	* *   some_standard_chimes_parameter_file.txt 

Note that the following must also be set in the main LAMMPS input file, to use ChIMES:

.. code-block:: text

    units       real		
    newton      on 		
    atom_style  atomic		
    atom_modify sort 0 0.0	


.. Warning::

    1. Implementation assumes outer cutoffs for (n+1)-body interactions are always :math:`\le` those for n-body interactions
    2. This capability is still under testing - please `let us know <https://groups.google.com/g/chimes_software>`_ if you observe strange behavior
    3. Assumes user wants single-atom energies to be added to the system energy. If you don't want to, zero the energy offsets in the parameter file
