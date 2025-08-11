.. _page-contributing:


Contributing to  ChIMES
==========================

The ChIMES calculator is an open source project, and we welcome contributions, e.g. bug fixes, updates to the documentation, extensions, etc. 

Contributions are made through the fork/pull request (PR) mechanism and generally, PRs should start from and target the develop branch. Additionally, PRs should include an attached test suite log file (see below).

Running the test suite
************************

To run the ChIMES calculator tests, simply navigate to ``serial_interface/tests/`` an run ``./run_tests.sh | tee run_tests.log``. 

.. Note::

    The ``run_tests.sh`` shell script assumes that a binary named ``python3.7`` exists in the users $PATH. If it does not exist, users can set the ``PYTH3`` variable near the top of ``run_tests.sh``  
    
.. Tip::

    The above command (i.e. ``./run_tests.sh | tee run_tests.log``) should be used generating a test suite log file for a PR, but if one desires quickers tests for debugging purposes, the test suite can be run as ``./run_tests.sh SHORT | tee run_tests.log``, which reduces the number of test calculations by a factor of roughly ten.
    
Making a test
************************

To add a new test, please navigate to ``serial_interface/tests/`` and follow the existing structure:

- Add your configuration file to the ``configurations/`` directory.
- Add the corresponding parameter file to the ``force_fields/`` directory.
- Add the expected output to the ``expected_output/`` directory.
- Update ``test_list.dat`` to include your new test.

Once everything is in place, running ``run_tests.sh`` should execute your added test.

LAMMPS contributions
************************

Currently, `LAMMPS <https://www.lammps.org/#gsc.tab=0>`_ version ``stable_29Aug2024_update1`` is supported. A set of tests is available in the ``etc/lmp/tests/`` directory. 
LAMMPS interacts with the following scripts: ``etc/lmp/src/pair_chimes.{h,cpp}``, ``etc/lmp/etc/pair.{h,cpp}`` and ``chimesFF/src/chimesFF.{h,cpp}`` Any contributing LAMMPS PRs should include an attached test suite log file and a test for the new functionality.

To add a new LAMMPS test, please navigate to ``etc/lmp/tests`` and follow the existing structure:

- Make a directory with a descriptive name for your test. e.g. ``test_<my_feature>``.
- In the directory, there should be a file named ``in.lammps``, a file named ``data.lammps``, and a parameter file.
- You should also make an ``expected_output`` directory that contains a file called ``out.lammps`` and ``log.lammps`` holiding the expected output.
.. Note::

    The ``in.lammps`` should follow exactly the format of the existing ``in.lammps`` files in other tests.
- Update the ``etc/tests/create.sh`` script so that it generates the expected output, follow the existing examples.
- Update ``etc/tests/run_tests.sh`` to include your new test directory in the list of tests to run.
.. Note::

    Note that all tests are made to run on a single node with 36 processors.

   


For additional questions and concerns, we can be contacted through our `Google group <https://groups.google.com/g/chimes_software>`_.



