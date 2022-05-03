.. _page-citing:


Citing ChIMES
=============

- :ref:`Reference Key for ChIMES Methods         <sec-methods>`
- :ref:`Reference Key for ChIMES Parameter Sets  <sec-params>`
- :ref:`Reference Key Definitions                <sec-defs>`

---------------


.. _sec-methods:

----------------------------------------
Reference Key for Methods/Applications
----------------------------------------

Key definitions are given :ref:`below <sec-defs>`.

=============================  =================
Method                         Reference Key
=============================  =================
2+3-body ChIMES                #. Carbon-1
ChIMES+DFTB                    #. PuH-DFTB
                               #. DNTF-DFTB
                               #. TiH-DFTB
                               #. QMD-DFTB
Iterative Refinement           #. CO-1
Carbon Condensation            #. CO-1
2+3+4-body ChIMES              #. CO-2
Distributed LASSO              #. CO-2
Active Learning                #. CO-2                 
ChIMES+MSST                    #. HN-1
                               #. DNTF-DFTB
=============================  =================

..                                #. DNTF-1


----------------

.. _sec-params:

---------------------------------
Reference Key for Parameter Sets
---------------------------------

Parameter set and key name are interchangeable. Key definitions are given :ref:`below <sec-defs>`.

=============  ====================================  ==========   ============================================   =========
KEY            Material                              Bodiedness   `T` (K)/ :math:`\rho` (gcc) Range              Comments
=============  ====================================  ==========   ============================================   =========
Carbon-1       Molten Carbon                         2            5000/2.43                                      N/A      
Carbon-1       Molten Carbon                         2            5000/2.43                                      N/A      
Carbon-1       Molten Carbon                         2+3          5000/2.43                                      N/A      
Carbon-1       Molten Carbon                         2+3          6000/2.25-3.00                                 N/A      
Water-1        Water                                 2+3          298/1.00                                       N/A
PuH-DFTB       Pu/H                                  2+3          0-300/N/A                                      DFTB :math:`E_{\mathrm{rep}}` 
CO-1           Carbon Monoxide (1:1)                 2+3          6500-9350/2.5                                  N/A
CO-2           Carbon Monoxide (1:1)                 2+3+4        2400/1.79                                      N/A
HN-1           Hydrazoic Acid H/N                    2+3+4        300-4500/1-2                                   N/A
DNTF-DFTB      3,4-bis(4-nitrofurazan-3-yl)furoxan   2+3          300-9000/1.86-3.4                              DFTB correction, Not applicable to other atom type ratios
TiH-DFTB       Ti/H                                  2+3          N/A/5.5                                      	 DFTB :math:`E_{\mathrm{rep}}` 
QMD-DFTB       C/N/O/H (based on QM database)        2+3          0/ambient                                      DFTB correction
=============  ====================================  ==========   ============================================   =========

.. DNTF-1 DNTF    3,4-bis(4-nitrofurazan-3-yl)furoxan   2+3+4        300-?/?-?                                      N/A 

---------

.. _sec-defs:

---------------------------------
Reference Key Definitions
---------------------------------

Corresponding authors are indicated with an asterisk (*).

==============   ==========================================================   ==============
Key              Link                                                         Definition
==============   ==========================================================   ==============
Carbon-1         (`link <https://doi.org/10.1021/acs.jctc.7b00867>`_)         R.K. Lindsey*, L.E. Fried, N. Goldman, `J. Chem. Theory Comput.`, **13**  6222   (2017).
PuH-DFTB         (`link <https://doi.org/10.1021/acs.jctc.8b00165>`_)         N. Goldman*, B. Aradi, R.K. Lindsey, L.E. Fried, `J. Chem. Theory Comput.` **14** 2652 (2018).
Water-1          (`link <https://doi.org/10.1021/acs.jctc.8b00831>`_)         R.K. Lindsey*, L.E. Fried, N. Goldman, `J. Chem. Theory Comput.`  **15**  436    (2019).
CO-1             (`link <https://doi.org/10.1063/5.0012840>`_)                R.K. Lindsey*, N. Goldman, L.E. Fried, S. Bastea, `J. Chem. Phys.` **153** 054103 (2020).
CO-2             (`link <https://doi.org/10.1063/5.0021965>`_)                R.K. Lindsey*, L.E. Fried, N. Goldman, S. Bastea, `J. Chem. Phys.` **153** 134117 (2020).
COND-1           (`link <https://doi.org/10.1038/s41467-019-14034-z>`_)       M.R. Armstrong*, R.K. Lindsey*, N. Goldman, M.H. Nielsen, E. Stavrou, L.E. Fried, J.M. Zaug, S. Bastea*, `Nat, Commun.` **11** 353 (2020).
HN-1             (`link <https://doi.org/10.1063/5.0029011>`_)                H. Pham*, R.K. Lindsey, L.E. Fried, N. Goldman, `J. Chem. Phys.` **153** 224102 (2020).
DNTF-DFTB        (`link <https://doi.org/10.26434/chemrxiv.14043839.v1>`_)    R.K. Lindsey*, S. Bastea*, N. Goldman, L. Fried, `J. Chem. Phys.` **154** 164115 (2021).
TiH-DFTB         (`link <https://doi.org/10.1021/acs.jctc.1c00172>`_)         N. Goldman*, K. Kweon, R.K. Lindsey, L.E. Fried, T.W. Heo, B. Sadigh, P. Soderlind, A. Landa, A. Perron, J. Jeffries, B. Wood, `J. Chem. Theory Comput.` **17** 4435 (2021).
QMD-DFTB	 (`link <https://doi.org/10.1021/acs.jpclett.2c00453>`_)     C.H. Pham*, R.K. Lindsey, L.E. Fried, N. Goldman, `J. Phys. Chem. Lett.` **13** 2934 (2022).
==============   ==========================================================   ==============

.. .. DNTF-1 DNTF      (`link <https://doi.org/?????????????????>`_)                R.K. Lindsey*, S. Bastea*, N. Goldman, L. Fried, ???.
