.. _page-getting_started-LLNL-BB:

Obtaining the code: LLNL Employees
=============================================

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

