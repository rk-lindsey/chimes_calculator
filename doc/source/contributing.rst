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
    

For additional questions and concerns, we can be contacted through our `Google group <https://groups.google.com/g/chimes_software>`_.



