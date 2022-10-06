.. chimes_calculator documentation master file, created by
   sphinx-quickstart on Mon Oct  5 15:20:18 2020.
   
.. Lines starting with two dots are comments if no special commands are found



ChIMES Calculator Documentation
=============================================

The **Ch**\ ebyshev **I**\ nteraction **M**\ odel for **E**\ fficient **S**\ imulation (ChIMES) is a machine-learned interatomic potential targeting chemistry in condensed phase systems. ChIMES models are able to approach quantum-accuracy through a systematically improvable explicitly many-bodied basis comprised of linear combinations of Chebyshev polynomials. Though originally developed to enable description of organic molecular materials, ChIMES has successfully been applied to systems spanning ambient water to molten carbon, and leveraged as correction for density functional based tight binding simulations. 

The ChIMES calculator comprises a flexible tool set for evaluating ChIMES interactions (e.g. in simulations, single point calculations, etc). Users have the option of directly embedding the ChIMES calculator within their codes (e.g. see ''The ChIMES Calculator,'' for advanced users), or evaluating interactions through the beginner-friendly serial interface, each of which have Python, C++, C, and Fortran API's. Files necessary for linking to popular simulation codes are being continually added with ancillary support. For more information see the links below.

The ChIMES Calculator is developed at Lawrence Livermore National Laboratory with funding from the US Department of Energy (DOE), and is open source, distributed under the terms of the LGPL v3.0 License.

`Note: This documentation is under still construction.`

.. toctree::
   :maxdepth: 1
   
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

This work was produced under the auspices of the U.S. Department of Energy by
Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

This work was prepared as an account of work sponsored by an agency of the
United States Government. Neither the United States Government nor Lawrence
Livermore National Security, LLC, nor any of their employees makes any warranty,
expressed or implied, or assumes any legal liability or responsibility for the
accuracy, completeness, or usefulness of any information, apparatus, product, or
process disclosed, or represents that its use would not infringe privately owned
rights. Reference herein to any specific commercial product, process, or service
by trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or Lawrence Livermore National Security, LLC. The views and
opinions of authors expressed herein do not necessarily state or reflect those
of the United States Government or Lawrence Livermore National Security, LLC,
and shall not be used for advertising or product endorsement purposes.
