.. _basics.ed:

******************************************
Introduction to Exact-diagonalization (ED)
******************************************

1. Second quantization
======================

1.1. Many-body Operators
------------------------

The exact-diagonalization (ED) algorithms implemented in edrixs are based on the second quantization language. In this framework, any many-body operator can be written as

.. math::

    \hat{O}=\sum_{\alpha,\beta}E_{\alpha,\beta}\hat{d}^{\dagger}_{\alpha}\hat{d}_{\beta} + \sum_{\alpha,\beta,\gamma,\delta} U_{\alpha,\beta,\gamma,
    \delta}\hat{d}^{\dagger}_{\alpha}\hat{d}^{\dagger}_{\beta}\hat{d}_{\gamma}\hat{d}_{\delta}

where, the first part are two-fermion terms and the second part are the four-fermion terms, with the matrix

.. math::

    E_{\alpha,\beta}=<\alpha|\hat{E}|\beta>

and the rank-4 tensor

.. math::

    U_{\alpha,\beta,\gamma,\delta}=<\alpha|<\beta|\hat{U}|\gamma>|\delta>

in the single-particle basis :math:`|\alpha>,|\beta>,|\gamma>,|\delta>`. Thus, the most important thing to perform calculations with edrixs is to define a single particle basis with particular orbital ordering and then write all the matrix or 4-rank tensor in the same single-particle basis.

1.2. Fock basis
---------------

The Fock basis can then be defined in the given single particle basis with specific orbital ordering. A Fock basis describes which orbital is occupied or empty. It can be written as a binary number with digital 1 representing occupied status and 0 empty, for example,

.. math::

    |110100>

where, from the left to right, the first, second and fourth orbital is occupied, while other orbitals are empty.


