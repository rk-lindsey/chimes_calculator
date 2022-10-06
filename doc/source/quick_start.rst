ChIMES Calculator Quickstart guide
=============================================

For more detailed instructions, see the :ref:`Getting Started <page-getting_started>` page.

Obtain a copy
######################

1. Create a fork of `the code <https://github.com/rk-lindsey/chimes_calculator>`_
2. Clone a copy of the code to your computer or high-performance computer (HPC)


Installing
######################

If your environment is correctly configured, you can install by simply executing ``./install.sh``.

If you are on a HPC using module files, you may need to load them first. Module files are already configured for a handful of HPC - inspect the contents of modfiles to see if
yours is listed. If it is (e.g., ``LLNL-LC.mod``), execute ``export hosttype=LLNL-LC; ./install.sh`` to install. Otherwise, load the appropriate modules by hand before running the
install script.

.. Note :: 

   Consider submitting module files and corresponding install.sh changes as a pull request, for your HPC!
   
Running
######################


You can test your installation by running an example job, e.g., by executing the following in your base chimes_calculator directory:


.. code-block:: bash
    
    serial_interface/examples/cpp/chimescalc \
    serial_interface/tests/force_fields/published_params.liqC.2b.cubic.txt \
    serial_interface/tests/configurations/liqC.2.5gcc_6000K.OUTCAR_#000.xyz | tee my_test.log 
 
