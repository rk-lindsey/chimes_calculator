.. chimes_calculator documentation master file, created by
   sphinx-quickstart on Mon Oct  5 15:20:18 2020.
   
.. Lines starting with two dots are comments if no special commands are found



ChIMES Calculator Documentation
=============================================

The **Ch**\ ebyshev **I**\ nteraction **M**\ odel for **E**\ fficient **S**\ imulation (ChIMES) is a machine-learned interatomic potential targeting chemistry in condensed phase systems. ChIMES models are able to approach quantum-accuracy through a systematically improvable explicitly many-bodied basis comprised of linear combinations of Chebyshev polynomials. Though originally developed to enable description of organic molecular materials, ChIMES has successfully been applied to systems spanning ambient water to molten carbon, and leveraged as correction for density functional based tight binding simulations. 

The ChIMES calculator comprises a flexible tool set for evaluating ChIMES interactions (e.g. in simulations, single point calculations, etc). Users have the option of directly embedding the ChIMES calculator within their codes (e.g. see ''The ChIMES Calculator,'' for advanced users), or evaluating interactions through the beginner-friendly serial interface, each of which have Python, C++, C, and Fortran API's. Files necessary for linking to popular simulation codes are being continually added with ancillary support.

The ChIMES Calculator was originally developed by R.K. Lindsey, N. Goldman, and L.E. Fried at Lawrence Livermore National Laboratory with funding from the US Department of Energy (DOE). Its development is now led by the Lindsey Lab in the Department of Chemical Engineering at the University of Michigan, Ann Arbor. The ChIMES Calculator is open source, distributed under the terms of the LGPL v3.0 License.

.. toctree::
   :maxdepth: 1
   :hidden:

   quick_start
   getting_started
   releases
   chimesFF
   serial_interface
   etc
   force_fields
   units
   utils
   citing
   contributing
   legal
   contact

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


---------------