.. _page-chimesFF:

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


For further details of the ChIMES ML-IAP equations, the reader is referred to the following. For a complete set of ChIMES references, see :ref:`Citing ChIMES <page-citing>`.

#. R.K. Lindsey*, L.E. Fried, N. Goldman,            `J. Chem. Theory Comput.`, **13**  6222   (2017). (`link <https://doi.org/10.1021/acs.jctc.7b00867>`_)
#. R.K. Lindsey*, L.E. Fried, N. Goldman,            `J. Chem. Theory Comput.`  **15**  436    (2019). (`link <https://doi.org/10.1021/acs.jctc.8b00831>`_)
#. R.K. Lindsey*, N. Goldman, L.E. Fried, S. Bastea, `J. Chem. Phys.`           **153** 054103 (2020). (`link <https://doi.org/10.1063/5.0012840>`_)
#. R.K. Lindsey*, L.E. Fried, N. Goldman, S. Bastea, `J. Chem. Phys.`           **153** 134117 (2020). (`link <https://doi.org/10.1063/5.0021965>`_)

`Corresponding authors are indicated with an asterisk (*).`

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

The ChIMES Calculator source files are located in ``chimesFF/src``. To use in a C++ code, simply ``#include "chimescalc.h"`` in the target code and instantiate a ``chimesFF`` object. As described in greater detail below, ``chimesFF`` objects take information on individual atom clusters and provide the corresponding ChIMES energy, stress tensor, and forces.  Any such code must at least include the following operations, in order:

    .. code-block:: cpp

       int my_rank = 0;
       chimesFF my_chimesFF_object;      // Instantiate
       my_chimesFF_object.init(my_rank); // Set MPI rank (replace with zero if used in serial code)
       my_chimesFF_object.read_parameters("my_parameter_file");

Note that the ChIMES calculator ``chimesFF`` class provides users with the following functions:

=========== =================  =================
Return Type Name               Arguments and Description
=========== =================  =================
void        init               ======   ===
                               Type     Description
                               ======   ===
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

double      max_cutoff_2B      ======    ===
                               Type      Description
                               ======    ===
                               bool      Flag: If true, prints largest 2-body cutoff
                               ======    ===

                               Returns the maximum 2-body outer cutoff distance.

double      max_cutoff_3B      ======    ===
                               Type      Description
                               ======    ===
                               bool      Flag: If true, prints largest 3-body cutoff
                               ======    ===

                               Returns the maximum 3-body outer cutoff distance.

double      max_cutoff_4B      ======    ===
                               Type      Description
                               ======    ===
                               bool      Flag: If true, prints largest 4-body cutoff
                               ======    ===

                               Returns the maximum 4-body outer cutoff distance.

void        compute_1B         ======    ===
                               Type      Description
                               ======    ===
                               int       Atom type index
                               double    Energy (updated)
                               ======    ===

                               Update energy with the single atom contribution.

void        compute_2B         ==========================   ===
                               Type                         Description
                               ==========================   ===
                               double                       Distance between two atoms, i and j
                               vector<double>               Distance vector components for each atom
                               vector<int>                  Type indices for atoms i and j
                               vector<vector<double* > >    Force pointer ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                               vector<double*>              Stress tensor pointer ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                               double                       Energy (updated by function)
                               ==========================   ===

                               Update the force pointer, stress tensor pointer, and energy with the two-atom contribution.

void        compute_3B         ==========================   ===
                               Type                         Description
                               ==========================   ===
                               vector<double>               Distances between three atoms, ij, ik, and jk
                               vector<vector<double> >      Distance vector components for each atom
                               vector<int>                  Type indices for atoms i, j and k
                               vector<vector<double* > >    Force pointer ([atom index (out of 3)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                               vector<double*>              Stress tensor pointer ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                               double                       Energy (updated by function)
                               ==========================   ===

                               Update the force pointer, stress tensor pointer, and energy with the three-atom contribution.

void        compute_4B         ==========================   ===
                               Type                         Description
                               ==========================   ===
                               vector<double>               Distance between four atoms, ij, ik, il, jk, jl, and kl
                               vector<vector<double> >      Distance vector components for each atom
                               vector<int>                  Type indices for atoms i, j, k  and l
                               vector<vector<double* > >    Force pointer ([atom index (out of 4)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                               vector<double*>              Stress tensor pointer ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                               double                       Energy (updated by function)
                               ==========================   ===

                               Update the force pointer, stress tensor pointer, and energy with the four-atom contribution.

=========== =================  =================



---------------


.. _sec-c-api:

The C API
^^^^^^^^^

The C API (``chimescalc_C*``) is located in ``chimesFF/api``. This wrapper provides C style name mangling and creates a
set of C-style wrapper functions. The latter are needed for compatibility with std::vector which is heavily used in ``chimesFF`` and not supported in most other languages. Any C code attempting to use the ChIMES Calculator should ``#include "chimescalc_C.h"``
and at least include the following operations, in order:

    .. code-block:: cpp

       int my_rank = 0;
       set_chimes();         // Instantiate
       init_chimes(my_rank); // Set MPI rank (replace with zero if used in serial code)
       chimes_read_params("my_parameter_file");

For additional information on compiling, see :ref:`Implementation Examples <sec-use-examples-api>`.

Note that the ChIMES calculator ``chimescalc_C`` API provides users with the following functions:

=========== ================================  =================
Return Type Name                              Arguments and Description
=========== ================================  =================
void        set_chimes                        No arguments. Instantiates a pointer to a ``chimesFF`` object.

void        init_chimes                       ======   ===
                                              Type     Description
                                              ======   ===
                                              int      MPI rank
                                              ======   ===

                                              Set the MPI rank. With the exception of error messages,
                                              the ChIMES calculator will only print output for rank 0.

void        chimes_read_params                ======   ===
                                              Type     Description
                                              ======   ===
                                              char*    Parameter file
                                              ======   ===

                                              Read the chimes parameter file.

int         get_chimes_2b_order               No arguments. Returns the two body order set by the parameter file.
int         get_chimes_3b_order               No arguments. Returns the three body order set by the parameter file.
int         get_chimes_4b_order               No arguments. Returns the four body order set by the parameter file.

double      get_chimes_max_2b_cutoff          No arguments. Returns the two body maximum outer cutoff set by the parameter file.
double      get_chimes_max_3b_cutoff          No arguments. Returns the three body maximum outer cutoff set by the parameter file.
double      get_chimes_max_4b_cutoff          No arguments. Returns the four body maximum outer cutoff set by the parameter file.


void        chimes_compute_2b_props           ============  ===
                                              Type          Description
                                              ============  ===
                                              double        Distance between two atoms, i and j
                                              double array  Distance vector components for each atom
                                              char*  array  Atom types for atoms i and j
                                              double array  Forces for atoms i and j ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                              double array  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                              double*       Energy (updated by function)
                                              ============  ===

                                              Update the force, stress tensor, and energy with the two-atom contribution.

void        chimes_compute_3b_props           ============  ===
                                              Type          Description
                                              ============  ===
                                              double array  Distances between three atoms, ij, ik, and jk
                                              double array  Distance vector components for each atom [atom][x, y, or z component]
                                              char*  array  Atom types for atoms i, j and k
                                              double array  Forces for atoms i, j, and k ([atom index (out of 3)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                              double array  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                              double*       Energy (updated by function)
                                              ============  ===

                                              Update the force, stress tensor, and energy with the three-atom contribution.

void        chimes_compute_4b_props           ============  ===
                                              Type          Description
                                              ============  ===
                                              double array  Distances between four atoms, ij, ik, il, jk, jl, and kl
                                              double array  Distance vector components for each atom [atom][x, y, or z component]
                                              char*  array  Atom types for atoms i, j, k  and l
                                              double array  Forces for atoms i, j, k, and l ([atom index (out of 4)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                              double array  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                              double*       Energy (updated by function)
                                              ============  ===

                                              Update the force, stress tensor, and energy with the four-atom contribution.


void        chimes_compute_2b_props_fromf90   ============  ===
                                              Type          Description
                                              ============  ===
                                              double*       Distance between two atoms, i and j
                                              double array  Distance vector components for each atom
                                              char*         Type for atom i
                                              char*         Type for atom j
                                              double array  Forces for atoms i and j ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                              double array  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                              double*       Energy (updated by function)
                                              ============  ===

                                              For calls from a Fortran code. Update the force, stress tensor, and energy with the two-atom contribution.

void        chimes_compute_3b_props_fromf90   ============  ===
                                              Type          Description
                                              ============  ===
                                              double        Distances between three atoms, ij, ik, and jk
                                              double array  Distance vector components for each atom [atom][x, y, or z component]
                                              char*         Type for atom i
                                              char*         Type for atom j
                                              char*         Type for atom k
                                              double array  Forces for atoms i, j, and k ([atom index (out of 3)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                              double array  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                              double*       Energy (updated by function)
                                              ============  ===

                                              For calls from a Fortran code. Update the force, stress tensor, and energy with the three-atom contribution.

void        chimes_compute_4b_props_fromf90   ============  ===
                                              Type          Description
                                              ============  ===
                                              double        Distances between four atoms, ij, ik, il, jk, jl, and kl
                                              double array  Distance vector components for each atom [atom][x, y, or z component]
                                              char*         Type for atom i
                                              char*         Type for atom j
                                              char*         Type for atom k
                                              char*         Type for atom l
                                              double array  Forces for atoms i, j, k, and l ([atom index (out of 4)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                              double array  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                              double*       Energy (updated by function)
                                              ============  ===

                                              For calls from a Fortran code. Update the force, stress tensor, and energy with the four-atom contribution.

=========== ================================  =================









---------------


.. _sec-fortran-api:

The Fortran API
^^^^^^^^^^^^^^^

The Fortran API (``chimescalc_F*``) is located in ``chimesFF/api``. This wrapper enables access to ``chimesFF`` functions
through the C API and handles other details like differences in array storage order.


Any Fortran code attempting to use the ChIMES Calculator should ``use chimescalc`` and at least include the following
operations, in order:

    .. code-block:: fortran

       integer(C_int) :: my_rank
       call f_set_chimes()         ! Instantiate
       call f_init_chimes(my_rank) ! Set MPI rank (replace with zero if used in serial code)
       call f_chimes_read_params(string2Cstring("my_parameter_file"))

For additional information on compiling, see :ref:`Implementation Examples <sec-use-examples-api>`.

Note that the ChIMES calculator ``chimescalc_F`` API provides users with the following functions:

=========== ==================================  =================
Return Type Name                                Arguments and Description
=========== ==================================  =================
none        f_chimes_compute_2b_props_fromf90   ==============   ===
                                                Type             Description
                                                ==============   ===
                                                C_double         Distance between two atoms, i and j
                                                C_double array   Distance vector components for each atom
                                                C_char           Type for atom i
                                                C_char           Type for atom j
                                                C_double array   Forces for atoms i and j ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                                C_double array   Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                                C_double         Energy (updated by function)
                                                ==============   ===

                                                Update the force, stress tensor, and energy with the two-atom contribution.

none        f_chimes_compute_3b_props_fromf90   ==============   ===
                                                Type             Description
                                                ==============   ===
                                                C_double array   Distances between three atoms, ij, ik, and jk
                                                C_double array   Distance vector components for each atom
                                                C_char           Type for atom i
                                                C_char           Type for atom j
                                                C_char           Type for atom k
                                                C_double array   Forces for atoms i, j, and k ([atom index (out of 3)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                                C_double array   Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                                C_double         Energy (updated by function)
                                                ==============   ===

                                                Update the force, stress tensor, and energy with the three-atom contribution.

none        f_chimes_compute_4b_props_fromf90   ==============   ===
                                                Type             Description
                                                ==============   ===
                                                C_double array   Distances between four atoms, ij, ik, il, jk, jl, and
                                                C_double array   Distance vector components for each atom
                                                C_char           Type for atom i
                                                C_char           Type for atom j
                                                C_char           Type for atom k
                                                C_char           Type for atom l
                                                C_double array   Forces for atoms i, j, k, and l ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                                C_double array   Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                                C_double         Energy (updated by function)
                                                ==============   ===

                                                Update the force, stress tensor, and energy with the four-atom contribution.

none        f_set_chimes                        No arguments. Instantiates a pointer to a ``chimesFF`` object.

none        f_init_chimes                       ======   ===
                                                Type     Description
                                                ======   ===
                                                int      MPI rank
                                                ======   ===

                                                Set the MPI rank. With the exception of error messages,
                                                the ChIMES calculator will only print output for rank 0.

none        f_chimes_read_params                ======   ===
                                                Type     Description
                                                ======   ===
                                                C_char   Parameter file
                                                ======   ===

                                                Read the chimes parameter file.


C_int       f_get_chimes_2b_order               No arguments. Returns the two body order set by the parameter file.
C_int       f_get_chimes_3b_order               No arguments. Returns the three body order set by the parameter file.
C_int       f_get_chimes_4b_order               No arguments. Returns the four body order set by the parameter file.

C_double    f_get_chimes_max_2b_cutoff          No arguments. Returns the two body maximum outer cutoff.
C_double    f_get_chimes_max_3b_cutoff          No arguments. Returns the three body maximum outer cutoff.
C_double    f_get_chimes_max_4b_cutoff          No arguments. Returns the four body maximum outer cutoff.

C_string    string2Cstring                      ======   ===
                                                Type     Description
                                                ======   ===
                                                string   Any text
                                                ======   ===

                                                Converts a Fortran string to a C_string

=========== ==================================  =================

---------------


.. _sec-python-api:

The Python API
^^^^^^^^^^^^^^

The Python API (``chimescalc_py*``) is located in ``chimesFF/api``. Like the Fortran API, this wrapper enables access to
``chimesFF`` functions through the C API, via ctypes.

Any python code attempting to use the ChIMES Calculator should ``import chimescalc_py`` and at least include the following
operations, in order:

    .. code-block:: python

       chimescalc_py.chimes_wrapper = chimescalc_py.init_chimes_wrapper("lib-C_wrapper-serial_interface.so") # Associate the wrapper with a compiled C API library file
       chimescalc_py.set_chimes()  # Instantiate
       chimescalc_py.init_chimes() # If run with MPI, an integer MPI rank can be passed to this function. By default, assumes rank = 0
       chimescalc_py.read_params("my_parameter_file")


For additional information on compiling (i.e. generation of ``lib-C_wrapper-serial_interface.so``), see :ref:`Implementation Examples <sec-use-examples-api>`.

Note that the ChIMES calculator ``chimescalc_py`` API provides users with the following functions:


=========== ==================================  =================
Return Type Name                                Arguments and Description
=========== ==================================  =================
ctypes      init_chimes_wrapper                 ==============   ===
                                                Type             Description
                                                ==============   ===
                                                str              C-wrapper library name (i.e. "lib-C_wrapper-serial_interface.so")
                                                ==============   ===

none        set_chimes                          No arguments. Instantiates a pointer to a ``chimesFF`` object.

none        init_chimes                         ==============   ===
                                                Type             Description
                                                ==============   ===
                                                int              MPI rank (optional parameter)
                                                ==============   ===

                                                Set the MPI rank. With the exception of error messages,
                                                the ChIMES calculator will only print output for rank 0.

none        read_params                         ==============   ===
                                                Type             Description
                                                ==============   ===
                                                str              Parameter file
                                                ==============   ===

float       get_chimes_max_2b_cutoff            No arguments. Returns the two body order set by the parameter file.
float       get_chimes_max_2b_cutoff            No arguments. Returns the three body order set by the parameter file.
float       get_chimes_max_2b_cutoff            No arguments. Returns the four body order set by the parameter file.

int         get_chimes_2b_order                 No arguments. Returns the two body maximum outer cutoff.
int         get_chimes_3b_order                 No arguments. Returns the three body maximum outer cutoff.
int         get_chimes_4b_order                 No arguments. Returns the four body maximum outer cutoff.

none        chimes_compute_2b_props             ==========  ===
                                                Type        Description
                                                ==========  ===
                                                float       Distances between atoms i and j
                                                float list  Distance vector components for each atom
                                                str list    Types for atom i and j
                                                float list  Forces for atoms i, and j ([atom index (out of 2)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                                float list  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                                float       Energy (updated by function)
                                                ==========  ===

                                                Update the force, stress tensor, and energy with the two-atom contribution.


none        chimes_compute_3b_props             ==========  ===
                                                Type        Description
                                                ==========  ===
                                                float list  Distances between three atoms, ij, ik, and jk
                                                float list  Distance vector components for each atom
                                                str list    Types for atom i, j, and k
                                                float list  Forces for atoms i, j, and k ([atom index (out of 3)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                                float list  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                                float       Energy (updated by function)
                                                ==========  ===

                                                Update the force, stress tensor, and energy with the three-atom contribution.

none        chimes_compute_4b_props              ==========  ===
                                                 Type        Description
                                                 ==========  ===
                                                 float list  Distances between four atoms, ij, ik, il, jk, jl, and kl
                                                 float list  Distance vector components for each atom
                                                 str list    Types for atom i, j, k, and l
                                                 float list  Forces for atoms i, j, k, and l ([atom index (out of 4)][component index (i.e. fx=0, fy=1, fz=3)]) (contents updated by function)
                                                 float list  Stress tensor ([s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz]) (contents updated by function)
                                                 float       Energy (updated by function)
                                                 ==========  ===

                                                Update the force, stress tensor, and energy with the four-atom contribution.

=========== ==================================  =================


---------------

.. _sec-use-examples-api:

Implementation Examples
^^^^^^^^^^^^^^^^^^^^^^^

The following codes demonstrates how ``chimesFF.{h,cpp}`` can be used to obtain the overall stress tensor, energy, and per-atom forces for a given system configuration using C, C++ Fortran, and Python. See the ``main.*`` files in each corresponding subdirectory of ``chimesFF/examples`` for further implementation details. Note that sample system configurations (i.e. ``*xyz`` files) and parameter files can be found in ``serial_interface/test/configurations`` and ``serial_interface/test/force_fields``, respectively. For user generated tests, note that ``*.xyz`` files must provide lattice vectors in the comment line, e.g. lx 0.0 0.0 0.0 ly 0.0 0.0 0.0 lz. Click :ref:`here <page-units>` for an overview of ChIMES units.

.. Note::

    All implementation examples are intended to be run on Unix-based systems (e.g. Linux, OSX).


.. Warning::

    These codes are for demonstrative purposes only and come with no guarantees.


.. Note::

    All example executables can be compiled at once via ``./install.sh`` from the ``chimes_calculator`` base directory, and similarly uninstalled via ``./uninstall.sh``. However, the examples below compile via the user-generated Makefiles located in each ``examples`` subdirectory, for demonstrative purposes.


* **C Example:** The ``main`` function of this example includes the C API, ``chimescalc_C.{h,cpp}``, which creates a global static pointer to a ``chimesFF`` object.
  The ``chimesFF`` pointer object is set up, i.e. by ``set_chimes()``, and used for access to ``chimesFF`` member functions, etc.

   * Navigate to ``chimesFF/examples/c``
   * Compile with: ``make all``
   * Test with: ``./C_wrapper-direct_interface <parameter file> <xyz file>``
   * Additional notes:

      * ``*.xyz`` files must not contain any information beyond atom type and x-, y-, and z- coordinate on coordinate lines.
      * This implementation does NOT use ghost atoms/layering thus the input system MUST have box lengths greater than two times the largest outer cutoff, or results will not be correct.

* **C++ Example:** The ``main`` function of this example creates an instance of ``serial_chimes_interface`` (i.e. a class inheriting ``chimesFF``,
  which computes energy, per-atom forces, and stress tensor for an overall system). For additional details, see :ref:`The ChIMES Calculator Serial Interface <page-serial_interface>`

   * Navigate to ``chimesFF/examples/cpp``
   * Compile with: ``make all``
   * Test with: ``./CPP-interface <parameter file> <xyz file> ``

* **Fortran Example:** Similar to the C example, this ``main`` function establishes a pointer to a ``chimesFF`` object via ``f_set_chimes()``.
  The ``f_set_chimes()`` function call is defined in ``chimescalc_F.f90,`` a wrapper for the C API ``chimescalc_C.cpp`` (i.e which facilitates C-style access to
  ``chimesFF`` member functions, etc). Actual linking is achieved at compilation. See the ``Makefile`` for details.

   * Navigate to ``chimesFF/examples/fortran``
   * Compile with: ``make all``
   * Test with: ``./fortran_wrapper-direct_interface <parameter file> <xyz file>``
   * Additional notes:

      * ``*.xyz`` files must not contain any information beyond atom type and x-, y-, and z- coordinate on coordinate lines.
      * This implementation does NOT use ghost atoms/layering thus the input system MUST have box lengths greater than two times the largest outer cutoff, or results will not be correct.

* **Python Example:** This example accesses ``chimesFF`` functions through `chimescalc_py.py``, a ctypes-based python API for access to the C API functions
  (i.e. through ``chimescalc_C.cpp``). Once ``chimescalc_py.py`` is imported, it is associated with a compiled C API library file, i.e. ``lib-C_wrapper-direct_interface.so`` and  can be used to access ``chimesFF`` member functions.

   * Navigate to ``chimesFF/examples/python``
   * Compile lib-C_wrapper-direct_interface.so with: ``make all``
   * Test with: python main.py <parameter file> <coordinate file>
   * Additional notes:

      * Requires ``lib-C_wrapper-direct_interface.so`` in the same directory, which is generated via ``make all``
      * Expects to be run with Python version 3.X
      * This implementation does NOT use ghost atoms/layering thus the input system MUST have box lengths greater than two times the largest outer cutoff, or results will not be correct.
