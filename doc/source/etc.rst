.. _etc:

ChIMES in LAMMPS
=============================================

We are currently working toward ChIMES calculator implementation in `LAMMPS <https://lammps.sandia.gov>`_ as a USER package. In the interim, we have provided a capability for compiling the `stable_29Aug2024_update1 <https://github.com/lammps/lammps.git build/>`_ version of LAMMPS with ChIMES as a pairstyle in this package.



Quick start
^^^^^^^^^^^^^^^^

.. Note ::

    To compile LAMMPS with ChIMES, you must have C++11 and MPI compilers avillable. As with installation of the ChIMES Calculator itself, if you are on a HPC using module files, you may need to load them first. Module files are already configured for a handful of HPC - inspect the contents of modfiles to see if yours is listed. If it is (e.g., LLNL-LC.mod), execute ``export hosttype=LLNL-LC`` Otherwise, load the appropriate modules by hand before running the install script.

    Note that Intel oneapi compilers (which are now free) can be used to properly configure your enviroment for all Intel capabilities (e.g., icc, mpiicpc, mkl, etc.) - simply locate and execute the setvars.sh script within your Intel installation.



Navigate to ``etc/lmp`` and execute  ``./install.sh`` to install. Once complete, the installation can be tested by navigating to ``etc/lmp/tests`` and either running the entire test suite via ``./run_tests.sh``, which takes roughly 15 minutes, or running an individual test by entering an example folder and executing ``../../exe/lmp_mpi_chimes -i in.lammps``. 

The install script compiles with LAMMPS packages ``manybody`` and ``extra-pair``. Additional packages can be included by adding appropriate commands to the ``etc/lmp/install.sh`` script. For example, to add the MOLECULE package, one would add the line highlighted in yellow below:


.. code-block :: 
    :lineno-start: 114
    :emphasize-lines: 6
    
    # Compile
    
    cd build/${lammps}/src
    make yes-manybody
    make yes-extra-pair
    make yes-molecule
    make -j 4 mpi_chimes
    cd -

    


.. Tip ::

    Additional flags can be specified during installation to enable features such as model tabulation and configuration fingerprint generation. For more details, see the :ref:`utils` page.
    

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


Examples
^^^^^^^^^^^^^^^^
