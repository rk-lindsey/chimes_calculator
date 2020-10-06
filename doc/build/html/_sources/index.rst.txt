.. chimes_calculator documentation master file, created by
   sphinx-quickstart on Mon Oct  5 15:20:18 2020.
   
.. Lines starting with two dots are comments if no special commands are found



ChIMES Calculator Documentation
=============================================

The **Ch**\ ebyshev **I**\ nteraction **M**\ odel for **E**\ fficient **S**\ imulation (ChIMES) is a machine-learned interatomic potential targeting chemistry in condensed phase systems.ChIMES models are able to approach quantum-accuracy through a systamtically improvable explicitly many-bodied basis comprised of linear combinations of Chebyshev polynomials. Though originally developed as to enable description of organic molecular materials, ChIMES has successfuly been applied to systems spanning ambient water to PuH, and leveraged as correction for density functional based tight binding simulations. 

The ChIMES calculator comprises a flexible toolset for running evaluating ChIMES interations (e.g. in simulations, single point calculations, etc). Users have the option of directly imbedding the ChIMES calculator within thier codes (e.g. see ''The ChIMES Calculator,'' for advanced users), or evaluating interactions through the beginner-friendly Serial Interface, each of which have Python, C++, C, and Fortran API's. Files necessary for linking to popular simulation codes including LAMMPS and DFTBplus are available as well. For mor information see the links below.

The ChIMES Calculator is developed at Lawrence Livermore National Laboratory with funding from the US Department of Energy (DOE), and is open source, distributed freely under the terms of the GNU Public License (GPL).

.. toctree::
   :maxdepth: 4
   
   chimesFF
   serial_interface
   etc
   force_fields
   citing
   units
   
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
