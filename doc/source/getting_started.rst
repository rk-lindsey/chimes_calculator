Getting Started
=============================================

* :ref:`Obtaining the code     <sec-obtaining>`
* :ref:`Compiling and running  <sec-compiling>`

---------------

.. _sec-obtaining:

Obtaining the code 
****************************************

LLNL Employees:
######################

Note that the ChIMES calculator is stored in a `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_calculator/browse>`_. The following steps must be completed to obtain/access:

* Send an e-mail to lindsey11 titled "ChIMES Calculator BB Access Request" to be granted access. 
* If the code is being copied to the LC: Set up SSH RSA keys

    * `Generate your keys <https://www.ssh.com/ssh/keygen/>`_. Use the default location and do not leave the pass-phrase blank.
    * Copy the contents of your key to the LLNL-hosted Bitbucket and add them to your account `(accessible here) <https://mybitbucket.llnl.gov/plugins/servlet/ssh/account/keys>`_.

Once access is granted and, if applicable, keys have been configured, the ChIMES calculator can be obtained three ways, described below:

1. Forking (recommended for developers/contributors)
2. Cloning 
3. Downloading as an archive

.. Note::

    All contributions to the ChIMES calculator repository proceed through a fork/pull request mechanism


1. Forking
^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_calculator/browse>`_ page, click the fork button (fourth button down from the top left, which looks like a two-prong trident). This will create a copy of ``chimes_calculator`` in your LLNL Bitbucket account. You can clone, edit, and commit to this repository as you see fit, and make any contributions to the original "parent" repository via the pull request mechanism. For additional details, see :ref:`Contributing <page-contributing>`. 


2. Cloning
^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_calculator/browse>`_ page, click the clone button (first button from the top left, which looks like a monitor with an arrow). Copy the ``chimes_calculator`` repository URL that appears. In a terminal, navigate to the desired download location, an type: ``git clone <the full copied URL> .``


2. Archive Download
^^^^^^^^^^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_calculator/browse>`_ page, click the three dots next to "master" under the heading "Source," and select download from the drop down menu



External Researchers:
######################

The ChIMES calculator is available on `Github <https://github.com/rk-lindsey/chimes_calculator>`_ and can be obtained three ways, described below:

1. Forking (recommended for developers/contributors)
2. Cloning 
3. Downloading as an archive

.. Note::

    All contributions to the ChIMES calculator repository proceed through a fork/pull request mechanism

.. Tip::

    For day-to-day use, we recommend using the most recent release branch. For developers, we recommend the `develop <https://github.com/rk-lindsey/chimes_calculator/tree/develop>`_ branch. For additional details, see :ref:`Releases <page-releases>` and :ref:`Contributing <page-contributing>`

1. Forking
^^^^^^^^^^

On the `Github repository <https://github.com/rk-lindsey/chimes_calculator>`_ page, click the fork button. This will create a copy of ``chimes_calculator`` in your Github account. You can clone, edit, and commit to this repository as you see fit, and make any contributions to the original "parent" repository via the pull request mechanism. For additional details, see :ref:`Contributing <page-contributing>`. 


2. Cloning
^^^^^^^^^^

On the `Github repository <https://github.com/rk-lindsey/chimes_calculator>`_ page, click the "code" button and select "clone". Copy the ``chimes_calculator`` repository URL that appears. In a terminal, navigate to the desired download location, an type: ``git clone <the full copied URL> .``


2. Archive Download
^^^^^^^^^^^^^^^^^^^

On the `Github repository <https://github.com/rk-lindsey/chimes_calculator>`_ page, click the "code" button and select "Download zip".




---------------


.. _sec-compiling:

Compiling and running the code
****************************************

As described above, the ``chimes_calculator`` comprises library tools for evaluating ChIMES interactions. However, the repository contains several usage examples (see, e.g. :ref:`ChIMES Calculator <sec-use-examples-api>` and :ref:`ChIMES Calculator Serial Interface <sec-ser-use-examples-api>` examples.). These examples can be compiled by navigating to a given example sub directory (e.g. ``chimes FF/examples/cpp/``) and typing ``make``, or executing the appropriate CMake commands by simply running ``./install.sh`` from the base ``chimes_calculator`` directory. If the latter option is used, a list of generated executables/library files and their respective install locations can be found in the generated ``build/install_manifest.txt`` file. Note that C++, C, *and* Fortran compilers are all required to use the ``./install.sh`` approach. 

Sample ChIMES parameter and input files are provided in the ``serial_interface/tests/force_fields`` and ``serial_interface/tests/configurations`` directories, allowing compiled executables to be tested via, e.g.:

.. code-block:: bash
    
    serial_interface/examples/cpp/chimescalc \
    serial_interface/tests/force_fields/published_params.liqC.2b.cubic.txt \
    serial_interface/tests/configurations/liqC.2.5gcc_6000K.OUTCAR_#000.xyz | tee my_test.log 
    

For additional details on using, integrating, and compiling, and contributing, see:

* :ref:`The ChIMES Calculator <page-chimesFF>`
* :ref:`The ChIMES Calculator Serial Interface <page-serial_interface>`
* :ref:`Contributing <page-contributing>`
