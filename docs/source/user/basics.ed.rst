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

    |110101>

where, from the left to right, the first, second, fourth and sixth orbital is occupied, while other orbitals are empty.

A fermion operator acts on it and get a new Fock state with proper sign, for example, annihilate one electron on the sixth orbital and create one electron on the third orbital,

.. math::

    \hat{d}_{6}&|110101> = (-1)^3 |110100>\\
    \hat{d}^{\dagger}_{3}\hat{d}_{6}&|110101> = (-1)^2(-1)^3|111100>

where, in :math:`(-1)^2` and :math:`(-1)^3`, 2 (3) means that it passes two (three) occupied orbitals when going from the first to the targeted orbital. Based on these operations, the matrix elements of any many-body operator can then be obtained.


2. Common orbital basis and ordering in edrixs
==============================================

2.1. Spin basis :math:`|s_z>`
----------------------------

In edrixs, we use a following convention of the spin ordering,

.. math::

    |\uparrow>, |\downarrow>, |\uparrow>, |\downarrow>, ..., |\uparrow>, |\downarrow>

2.2. Complex spherical harmonics :math:`|l_z, s_z>`
--------------------------------------------------

For an orbital angular momentum :math:`l`, the orbital ordering are:

.. math::

    l=0: &|0, \uparrow>, |0, \downarrow>\\
    l=1: &|-1, \uparrow>, |-1, \downarrow>, |0, \uparrow>, |0, \downarrow>, |+1, \uparrow>, |+1, \downarrow>\\
    l=2: &|-2, \uparrow>, |-2, \downarrow>, |-1, \uparrow>, |-1, \downarrow>, |0, \uparrow>, |0, \downarrow>, |+1, \uparrow>, |+1, \downarrow>,\\
         &|+2, \uparrow>, |+2, \downarrow>\\
    l=3: &|-3, \uparrow>, |-3, \downarrow>, |-2, \uparrow>, |-2, \downarrow>, |-1, \uparrow>, |-1, \downarrow>, |0, \uparrow>, |0, \downarrow>, |+1, \uparrow>, |+1, \downarrow>,\\
         &|+2, \uparrow>, |+2, \downarrow>, |+3, \uparrow>, |+3, \downarrow>

2.3. Real harmonics
-------------------

.. math::

    p: &|p_x, \uparrow>, |p_x, \downarrow>, |p_y, \uparrow>, |p_y, \downarrow>, |p_z, \uparrow>, |p_z, \downarrow>\\
    t_{2g}: &|d_{xz}, \uparrow>, |d_{xz}, \downarrow>, |d_{yz}, \uparrow>, |d_{yz}, \downarrow>, |d_{xy}, \uparrow>, |d_{xy}, \downarrow>\\
    d: &|d_{3z^2-r^2, \uparrow>, |d_{3z^2-r^2, \downarrow}, |d_{xz}, \uparrow>, |d_{xz}, \downarrow>, |d_{yz}, \uparrow>, |d_{yz}, \downarrow>,\\
       &|d_{x^2-y^2}, \uparrow>, |d_{x^2-y^2}, \downarrow>, |d_{xy}, \uparrow>, |d_{xy}, \downarrow>\\
    f: &|f_{z^3}, \uparrow>, |f_{z^3}, \downarrow>, |f_{xz^2}, \uparrow>, |f_{xz^2}, \downarrow>, |f_{yz^2}, \uparrow>, |f_{yz^2}, \downarrow>\\
       &|f_{z(x^2-y^2)}, \uparrow>, |f_{z(x^2-y^2)}, \downarrow>, |f_{xyz}, \uparrow>, |f_{xyz}, \downarrow>,\\
       &|f_{x(x^2-3y^2)}, \uparrow>, |f_{x(x^2-3y^2)}, \downarrow>, |f_{y(3x^2-y^2)}, \uparrow>, |f_{y(3x^2-y^2)}, \uparrow>

2.4. Real cubic harmonics
-------------------------

The real cubic harmonics basis is the irreducible representation of the cubic point group :math:`O_h`,

.. math::

    T_{1u}: &\{p_x, p_y, p_z\},\\
            &\{-\frac{\sqrt{6}}{4}f_{xz^2}+\frac{\sqrt{10}}{4}f_{x(x^2-3y^2)}, -\frac{\sqrt{6}}{4}f_{yz^2}-\frac{\sqrt{10}}{4}f_{y(3x^2-y^2)},f_{z^3}\}\\
    T_{2u}: &\{-\frac{\sqrt{10}}{4}f_{xz^2}-\frac{\sqrt{6}}{4}f_{x(x^2-3y^2)}, \frac{\sqrt{10}}{4}f_{yz^2}-\frac{\sqrt{6}}{4}f_{y(3x^2-y^2)},f_{z(x^2-y^2)}\}




