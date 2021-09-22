Support for Linking with External Codes
=============================================

Using the ChIMES Calculator with LAMMPS
******************************************

We are currently working to get toward ChIMES calculator implementation in `LAMMPS <https://lammps.sandia.gov>`_ as a USER package. In the interim, the following provides a guide to implementing the ChIMES calculator as a LAMMPS pairstyle.

.. Note::

    This example assumes users have downloaded the 29 Oct 2020 version of LAMMPS (stable version as of 9/21/2021), which can be downloaded `here <https://lammps.sandia.gov/download.html>`_. 

Compiling
^^^^^^^^^^^^^^^^

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

        If you are using an intel compiler, either delete the ``pair_list.*`` files that appear in the src folder following the ``make yes-manybody`` command, or add ``-restrict`` to ``CCFLAGS`` in ``MAKE/Makefile.mpi_chimes``. 


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
    atom_modify map array

.. Warning::

    1. Implementation assumes the max 2-body cutoff is always larger than the max (n>2)-body cutoffs; this is required for compatibility with LAMMPS
    2. This capability is still under testing - please `let us know <https://groups.google.com/g/chimes_software>`_ if you observe strange behavior
    3. Assumes user wants single-atom energies to be added to the system energy. If you don't want to, zero the energy offsets in the parameter file
    4. See comments at the top of chimesFF.{h.cpp} and pair_chimes.{h,cpp} for additional warnings
