.. _page-getting_started:

Getting Started
=============================================

.. _sec-obtaining:

Obtaining the code 
****************************************

LLNL Employees:
######################

Please see :ref:`the special Bitbucket instructions <page-getting_started-LLNL-BB>`


GitHub Users:
######################

Please see :ref:`the special Github instructions <page-getting_started-Open-GH>`



---------------


.. _sec-compiling:

Compiling and running the code
****************************************

As described above, the ``chimes_calculator`` comprises library tools for evaluating ChIMES interactions. Further information about the model architecture can be found in :ref:`ChIMES Calculator <page-chimesFF>`. ``chimesFF`` (expert mode) is the source file that implements the core logic of calculating ChIMES interactions given the atom cluster. It is focused on individual atom clusters  and does not perform neighbor clustering. 

An example of using chimesFF is provided in :ref:`ChIMES Calculator Serial Interface <sec-ser-use-examples-api>`, which demonstrates how to use the ``chimes_calculator`` for an overall system evaluation. This implementation is intended to be instructional and is not optimized for performance.

These examples can be compiled by navigating to a given example subdirectory (e.g. ``chimesFF/examples/cpp/``) and typing ``make``.

You can compile and install the entire software suite in different computing environments. All require C, C++ and Fortran compilers:

 * Standard Installation (Local Machine or Basic Environment)

      *  If your environment is correctly configured, you can simply execute ``./install.sh`` in the outer most directory.


 * Installation on an HPC cluster

      * You need to load pre-configured module files first. Module files are already configured for a handful of HPC in ``/modfiles/`` ; verify if yours is listed. If it is (e.g., LLNL-LC.mod), execute ``export hosttype=LLNL-LC; ./install.sh`` to install. Otherwise, load the appropriate modules by hand before running the install script.

 * CMake-Based Installation

      * You can execute the appropriate CMake commands by simply running ``./install.sh`` from the base       ``chimes_calculator`` directory. If this option is used, a list of generated executables/library files and their respective install locations can be found in the generated ``build/install_manifest.txt`` file. 

.. Note :: 

   Consider submitting module files and corresponding install.sh changes as a PR, for your HPC! 


Sample ChIMES parameter and input files are provided in the ``serial_interface/tests/force_fields`` and ``serial_interface/tests/configurations`` directories, allowing compiled executables to be tested via, e.g.:

.. code-block:: bash
    
    serial_interface/examples/cpp/chimescalc \
    serial_interface/tests/force_fields/published_params.liqC.2b.cubic.txt \
    serial_interface/tests/configurations/liqC.2.5gcc_6000K.OUTCAR_#000.xyz | tee my_test.log 

Compiling and running LAMMPS with ChIMES
****************************************

We have implemented the CHIMES in the `LAMMPS <https://lammps.sandia.gov/>`_ and further information regarding Compiling and running can be found in :ref:`The ChIMES Calculator with LAMMPS <page-etc.rst>`. In here, the ChIMES force field is compiled into LAMMPS as a custom pair style, allowing it to be used like any other interatomic potential within LAMMPS input scripts.
    

For additional details on using, integrating, compiling, and contributing, see:

* :ref:`The ChIMES Calculator <page-chimesFF>`
* :ref:`The ChIMES Calculator Serial Interface <page-serial_interface>`
* :ref:`The ChIMES Calculator with LAMMPS <page-etc.rst>`
* :ref:`Contributing <page-contributing>`
