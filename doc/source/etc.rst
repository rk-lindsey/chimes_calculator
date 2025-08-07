.. _page-etc:

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

**Important Considerations:**

Atom Type Mass Matching
"""""""""""""""""""""

The element masses defined in the ChIMES parameter file must match the masses defined in the LAMMPS data file. The matching is done by comparing masses to within 0.001 atomic mass units. The order of atom types in the files does not matter - only the mass values need to match.

For example, if your ChIMES parameter file defines elements with masses 12.011 and 15.999, then your LAMMPS data file must have atom types with matching masses, regardless of the order they appear in the files.

**ChIMES Parameter File Example:**
The parameter file defines atom types with specific masses:

.. code-block:: text

    ATOM TYPES: 2

    # TYPEIDX #     # ATM_TYP #     # ATMCHRG #     # ATMMASS #
    0               C               0               12.011
    1               O               0               15.999

**LAMMPS Data File Example:**
The data file must have atom types with matching masses (order doesn't matter):

.. code-block:: text
  
    2 atom types

    Masses

    1 12.0107  # C (matches ChIMES mass 12.011 within tolerance)
    2 15.9994  # O (matches ChIMES mass 15.999 within tolerance)

**Important:** If no mass matches are found between the ChIMES parameter file and LAMMPS data file, the simulation will terminate with an error, as ChIMES cannot be used for any interactions.

Hybrid Overlay Usage
"""""""""""""""""""

ChIMES can be used simultaneously with other potentials using LAMMPS' hybrid/overlay pair style. This allows you to combine ChIMES with additional force fields for specific interactions.

**Example of hybrid/overlay usage:**

.. code-block:: text

    pair_style      hybrid/overlay chimesFF momb 9.0 0.75 20.0 lj/cut 10.0
    pair_coeff      * * chimesFF ${param_file}
    pair_coeff      1 1 momb 0.0 1.0 1.0 418.26 2.904
    pair_coeff      1 2 lj/cut 0.25   3.5
    pair_coeff      2 2 lj/cut 0.25   3.5

In this example, ChIMES is combined with MOMB (Many-body van der Waals) and Lennard-Jones potentials for different atom type interactions. 

Models with D2 Dispersion Correction
"""""""""""""""""""""""""""""""""""

Using hybrid ChIMES and MOMB is specifically for adding D2 dispersion correction at the time of using LAMMPS. Whether to include MOMB and the parameters used is specific to the ChIMES parameterization and should not be added arbitrarily.

**Important:** When using MOMB with ChIMES, you must include the ``make yes-extra-pair`` command in the install.sh script when compiling LAMMPS to enable the MOMB potential support.

.. code-block:: text

    # Compile

    cd build/${lammps}/src
    make yes-manybody
    make yes-extra-pair

.. Warning::

    1. Implementation assumes outer cutoffs for (n+1)-body interactions are always :math:`\le` those for n-body interactions
    2. This capability is still under testing - please `let us know <https://groups.google.com/g/chimes_software>`_ if you observe strange behavior
    3. Assumes user wants single-atom energies to be added to the system energy. If you don't want to, zero the energy offsets in the parameter file

