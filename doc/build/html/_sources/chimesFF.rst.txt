The ChIMES Calculator
=====================

Overview
********

ChIMES is a reactive explicitly many-bodied machine learned interatomic potential (ML-IAP) for which interactions are 
computed on the basis of atom clusters. For example, the total ChIMES energy is given as:

.. math::
   :nowrap:

   \begin{eqnarray}
      E_{n_\mathrm{B}} 
      = \sum_{i_1}        ^{n_a} {}^{1}\!E_{i_1}
      + \sum_{i_1>i_2}    ^{n_a} {}^{2}\!E_{i_1i_2} 
      + \sum_{i_1>i_2>i_3}^{n_a} {}^{3}\!E_{i_1i_2i_3} 
      + \dots
      + \sum_{i_1>i_2\dots >i_{n\mathrm{B}-1}>i_{n\mathrm{B}}}^{n_a} {}^{n}\!E_{i_1i_2\dots i_n}
   
   \end{eqnarray}

   
where :math:`E_{n_\mathrm{B}}` is the total ChIMES system energy, :math:`n_{\mathrm{B}}` is the maximum bodiedness, 
:math:`{}^{n}\!E_{i_1i_2\dots i_n}` is the :math:`n`-body ChIMES energy for 
a given set of :math:`n` atoms with indices :math:`i = {i_1, i_2, \dots , i_n}`, and :math:`n_a` is the total number of atoms in the system. In the 
ChIMES framework, single-body energies are constant values and :math:`n`-body energies are constructed from the product of
polynomials of transformed atom pair distances. Thus, a 2-body interaction would involve a single pair, :math:`ij`, while a 
three-body interaction would involve three pairs, :math:`ij, ik,` and :math:`jk`, a 4-body interaction would involve :math:`4\choose 2` pairs, 
and so on. Currently, the ChIMES calculator supports up to 4-body interactions.


For further details of the ChIMES ML-IAP, the reader is referred to the following. For a complete set of ChIMES references, see :ref:`Citing ChIMES <page-citing>`.

#. R.K. Lindsey*, L.E. Fried, N. Goldman, S. Bastea, `J. Chem. Phys.`           **153** ??     (2020). (`link <http://google.com>`_)
#. R.K. Lindsey*, N. Goldman, L.E. Fried, S. Bastea, `J. Chem. Phys.`           **153** 054103 (2020). (`link <https://doi.org/10.1063/5.0012840>`_)
#. R.K. Lindsey*, L.E. Fried, N. Goldman,            `J. Chem. Theory Comput.`  **15**  436    (2019). (`link <https://doi.org/10.1021/acs.jctc.8b00831>`_)
#. R.K. Lindsey*, L.E. Fried, N. Goldman,            `J. Chem. Theory Comput.`, **13**  6222   (2017). (`link <https://doi.org/10.1021/acs.jctc.7b00867>`_)


Sections
********

* :ref:`ChIMES Calculator <sec-chimes-calc>`
* :ref:`C API <sec-c-api>`
* :ref:`Fortran API <sec-fortran-api>`
* :ref:`Python API <sec-python-api>`
* :ref:`Implementation Examples <sec-use-examples-api>`

---------------

.. _sec-chimes-calc:

The ChIMES Calculator
*********************

The ChIMES Calculator source files are located in chimesFF/src. To use in a C++ code, simply include ``chimesFF.h`` in the target 
code and instantiate a ``chimesFF`` object. Any such code must at least include the following operations, in order:
      
     .. code-block:: cpp

	 int my_rank = 0;
	 chimesFF my_chimesFF_object;	   // Instantiate
	 my_chimesFF_object.init(my_rank); // Set MPI rank (replace with zero if used in serial code)
	 my_chimesFF_object.read_parameters("my_parameter_file"); 
	 
Note that the ChIMES calculator ``chimesFF`` class provides users with the following functions:	 

=========== =================  =================
Return Type Name               Arguments and Description
=========== =================  =================
void        init               ======	===	    
			       Type	Description
			       ======	===
			       int      MPI rank
			       ======   ===
			       Set the MPI rank. With the exception of error messages, 
			       the ChIMES calculator will only print output for rank 0.
			       
void        read_parameters    ======   ===
                               Type     Description
			       ======   ===
                               string   Parameter file
			       ======   ===
			       Read the chimes parameter file.

void        set_atomtypes      ==============  ===
			       Type            Description
			       ==============  ===
			       vector<string>  List of atom types defined by parameter file (updated by function)
			       ==============  ===
			       Update the input vector with atom types in the parameter file.

double      max_cutoff_2B      ======	===
			       Type	Description
			       ======	===
			       bool	Flag: If true, prints largest 2-body cutoff
			       ======	===
			       Returns the maximum 2-body outer cutoff distance.
			       
double      max_cutoff_3B      ======	===
			       Type	Description
			       ======	===
			       bool	Flag: If true, prints largest 3-body cutoff
			       ======	===
			       Returns the maximum 3-body outer cutoff distance.
			       
double      max_cutoff_4B      ======	===
			       Type	Description
			       ======	===
			       bool	Flag: If true, prints largest 4-body cutoff
			       ======	===
			       Returns the maximum 4-body outer cutoff distance.
			       
void        compute_1B         ======	===
			       Type	Description
			       ======	===
			       int	Atom type index
			       double   Energy (updated)
			       ======	===
			       Update energy with the single atom contribution.
			       
void        compute_2B         ==========================   ===
			       Type	Description
			       ==========================   ===
			       double                       Distance between two atoms, i and j
			       vector<double>               Distance vector components for each atom
			       vector<int>                  Type indices for atoms i and j 
			       vector<vector<double* > >    Force pointer ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)])
			       vector<double*>              Stress tensor pointer ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz])
			       double                       Energy
			       ==========================   ===
			       Update the force pointer, stress tensor pointer, and energy with the two-atom contribution.
			       
void        compute_3B         ==========================   ===
			       Type	Description
			       ==========================   ===
			       vector<double>               Distance between three atoms, i, j, and k
			       vector<vector<double> >      Distance vector components for each atom
			       vector<int>                  Type indices for atoms i, j and k
			       vector<vector<double* > >    Force pointer ([atom index (out of 3)][component index (i.e. fx=0, fy=1, fz=3)])
			       vector<double*>              Stress tensor pointer ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz])
			       double                       Energy
			       ==========================   ===
			       Update the force pointer, stress tensor pointer, and energy with the three-atom contribution.
			       
			       
void        compute_4B         ==========================   ===
			       Type	Description
			       ==========================   ===
			       vector<double>		    Distance between four atoms, i, j, k and l
			       vector<vector<double> >      Distance vector components for each atom
			       vector<int>		    Type indices for atoms i, j, k  and l
			       vector<vector<double* > >    Force pointer ([atom index (out of 4)][component index (i.e. fx=0, fy=1, fz=3)])
			       vector<double*>  	    Stress tensor pointer ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz])
			       double			    Energy
			       ==========================   ===
			       Update the force pointer, stress tensor pointer, and energy with the four-atom contribution.
			       
=========== =================  =================


---------------


.. _sec-c-api:

The C API
^^^^^^^^^

The C API (``wrapper-C*``) is located in ``chimesFF/api``. This wrapper provides C style name mangling and creates a 
set of C-style wrapper functions. The latter are needed because ``chimesFF`` natively uses vectors, which are not 
supported in most other languages.


---------------


.. _sec-fortran-api:

The Fortran API
^^^^^^^^^^^^^^^

The Fortran API (``wrapper-F*``) is located in ``chimesFF/api``. This wrapper enables access to ``chimesFF`` functions 
through the C API and handles other details like differences in array storage order. 


---------------


.. _sec-python-api:

The Python API
^^^^^^^^^^^^^^

The Python API (``wrapper_py*``) is located in ``chimesFF/api``. Like the Fortran API, this wrapper enables access to 
``chimesFF`` functions through the C API, via ctypes. 



---------------

.. _sec-use-examples-api:

Implementation Examples
^^^^^^^^^^^^^^^^^^^^^^^

The following codes demonstrates how ``chimesFF{h,cpp}`` can be used to obtain the
stress tensor, energy, and per-atom forces for a given system configuration using C, C++ 
Fortran, and python. See the ``main.*`` files in each corresponding subdirectory of ``chimesFF/examples``
for further implementation details. Note that sample system configurations (i.e. ``*xyz`` files) and 
parameter files can be found in ``serial_interface/test/configurations`` and ``serial_interface/test/force_fields``, respectively. For user
generated tests, note that ``*.xyz`` files must be fore orthorhombic systems, with lattice vectors provided in the comment line,
e.g. lx 0.0 0.0 0.0 ly 0.0 0.0 0.0 lz.

Disclaimer: These codes are for demonstrative purposes only and come with no guarantees.


* **C Example:** The ``main`` function of this example includes the C API, ``wrapper-C.{h,cpp}``, which creates a global static pointer to a ``chimesFF`` object. 
  The ``chimesFF`` pointer object is set up, i.e. by ``set_chimes()``, and used for access to ``chimesFF`` member functions, etc.
     
   * Compile with: ``make all``
   * Test with: ``./test_wrapper-C <parameter file> <xyz file>``
   * Additional notes: 
   
      * ``*.xyz`` files must not contain any information beyond atom type and x-, y-, and z- coordinate on coordinate lines.
      * This implementation does NOT use ghost atoms/layering thus the input system MUST have box lengths greater than two times the largest outer cutoff, or results will not be correct.
      
* **C++ Example:** The ``main`` function of this example creates an instance of ``serial_chimes_interface`` (i.e. a class inheriting ``chimesFF``, 
  which computes energy, per-atom forces, and stress tensor for an overall system). For additional details, see [**REFERENCE THE SERIAL CHIMES INTERFACE**] 
   
   * Compile with: ``make test-CPP``
   * Test with: ``./test-CPP <parameter file> <xyz file> <nlayers>``
   * Additional notes: 
   
   	* ``<nlayers>`` is an integer definiting how many layers of ghost atoms should be used. Note that ``<nlayers>`` should be large enough to prevent self-interaction across the periodic boundary.
	
* **Fortran Example:** Similar to the C example, this ``main`` function establishes a pointer to a ``chimesFF`` object via ``f_set_chimes()``. 
  The ``f_set_chimes()`` function call is defined in ``wrapper-F.F90,`` a wrapper for the C API ``wrapper-C.cpp`` (i.e which facilitates C-style access to 
  ``chimesFF`` member functions, etc). Actual linking is achieved at compilation. See the ``Makefile`` for details. 
  
   * Compile with: ``make all``
   * Test with: ``./test_wrapper-F <parameter file> <xyz file>``
   * Additional notes: 
   
      * ``*.xyz`` files must not contain any information beyond atom type and x-, y-, and z- coordinate on coordinate lines.
      * This implementation does NOT use ghost atoms/layering thus the input system MUST have box lengths greater than two times the largest outer cutoff, or results will not be correct.
      
* **Python Example:** This example accesses ``chimesFF`` functions through ``wrapper_py.py``, a ctypes-based python API for access to the C API functions 
  (i.e. through ``wrapper-C.cpp``). Once ``wrapper_py.py`` is imported, it is associated with a compiled C API library file, i.e. ``libwrapper-C.so`` and 
  can be used to access ``chimesFF`` member functions. 

   * Compile libwrapper-C.so with: ``make all``
   * Test with: python3 main.py <parameter file> <coordinate file>
   * Additional notes: 
   
      * Requires ``libwrapper-C.so`` in the same directory, which is generated via ``make all``
      * Expects to be run with Python version 3.X
      * This implementation does NOT use ghost atoms/layering thus the input system MUST have box lengths greater than two times the largest outer cutoff, or results will not be correct.

---------------
