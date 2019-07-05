__all__ = ['HR', 'KVec', 'SymKVec', 'UniKVec']

import numpy as np


class HR():
    """
    Class for post-processing of the `Wannier90 <http://www.wannier.org/>`_
    tight-binding (TB) Hamiltonian in real space.
    """

    def __init__(self, nwann, nrpt, irpt0, rpts, deg_rpt, hr):

        """
        Parameters
        ----------
        nwann: int
            Number of Wannier orbitals.
        nrpt: int
            Number of :math:`r` points.
        irpt0: int
            Index of the origin :math:`r`-point :math:`(0, 0, 0)`.
        rpts: 3 elements of 1d float array
            The fractional coordinates of :math:`r`-points wrt. primitive lattice vectors.
        deg_rpt: 1d int array
            Degenerancy of :math:`r`-points.
        hr: 3d complex array
            TB Hamiltonian :math:`H(r)` from Wannier90.
        """
        self.nwann = nwann
        self.nrpt = nrpt
        self.irpt0 = irpt0
        self.deg_rpt = deg_rpt
        self.rpts = rpts
        self.hr = hr

    @staticmethod
    def from_file(fname='wannier90_hr.dat'):
        """
        Generate a TB Hamiltonian from Wannier90 output file ``case_hr.dat``.

        Parameters
        ----------
        fname: string
            The file that contains the Wannier90 output file ``case_hr.dat``.

        Returns
        -------
        HR: HR object
            An instance of HR object.
        """

        with open(fname, 'r') as f:
            # skip the header (1 line)
            f.readline()
            nwann = int(f.readline().strip())
            nrpt = int(f.readline().strip())
            # the degeneracy of R points
            nline = nrpt // 15 + 1
            tmp = []
            for i in range(nline):
                tmp.extend(f.readline().strip().split())
            tmp = [np.int(item) for item in tmp]
            deg_rpt = np.array(tmp, dtype=np.int)
            # read hr for each r-point
            rpts = np.zeros((nrpt, 3), dtype=np.int)
            hr = np.zeros((nrpt, nwann, nwann), dtype=np.complex128)
            for i in range(nrpt):
                for j in range(nwann):
                    for k in range(nwann):
                        rx, ry, rz, hr_i, hr_j, hr_real, hr_imag = f.readline().strip().split()
                        rpts[i, :] = int(rx), int(ry), int(rz)
                        if int(rx) == 0 and int(ry) == 0 and int(rz) == 0:
                            irpt0 = i
                        hr[i, k, j] = np.float64(hr_real) + np.float64(hr_imag) * 1j
            # construct the HR instance
            return HR(nwann, nrpt, irpt0, rpts, deg_rpt, hr)

    @staticmethod
    def copy_hr(other):
        """
        Copy an instance of HR.

        Parameters
        ----------
        other: HR object
            A HR object to be copied.

        Returns
        -------
        HR: HR object
            Return a new HR object.
        """

        return HR(other.nwann,
                  other.nrpt,
                  other.irpt0,
                  np.copy(other.rpts),
                  np.copy(other.deg_rpt),
                  np.copy(other.hr))

    def get_hr0(self, ispin=0):
        """
        Return the on-site term :math:`H(r=0)`.

        Parameters
        ----------
        ispin: int
            How to include spin degree of freedom:

            - 0: do not include spin degree of freedom manually

            - 1: include spin degree of freedom manually,

            spin order: up, dn, up, dn, ..., up, dn

            - 2: include spin degree of freedom manually,

            spin order: up, up, up, ..., dn, dn, dn, ...

        Returns
        -------
        hr0: 2d complex array
            The on-site Hamiltonian.
        """

        if ispin == 1:
            norbs = 2 * self.nwann
            hr0 = np.zeros((norbs, norbs), dtype=np.complex128)
            hr0[0:norbs:2, 0:norbs:2] = self.hr[self.irpt0, :, :]
            hr0[1:norbs:2, 1:norbs:2] = self.hr[self.irpt0, :, :]
            return hr0
        elif ispin == 2:
            norbs = 2 * self.nwann
            hr0 = np.zeros((norbs, norbs), dtype=np.complex128)
            hr0[0:self.nwann, 0:self.nwann] = self.hr[self.irpt0, :, :]
            hr0[self.nwann:norbs, self.nwann:norbs] = self.hr[self.irpt0, :, :]
            return hr0
        else:
            return self.hr[self.irpt0, :, :]

    def get_hr(self, ispin):
        """
        Return TB Hamiltonian of :math:`H(r)`.

        Parameters
        ----------
        ispin: int
            How to include spin degree of freedom:

            - 0: do not include spin degree of freedom manually

            - 1: include spin degree of freedom manually,

            spin order: up, dn, up, dn, ..., up, dn

            - 2: include spin degree of freedom manually,

            spin order: up, up, up, ..., dn, dn, dn, ...

        Returns
        -------
        hr: 3d complex array
            The Hamiltonian :math:`H(r)`.
        """

        # with spin, spin order: up dn up dn up dn
        if ispin == 1:
            norbs = 2 * self.nwann
            hr_spin = np.zeros((self.nrpt, norbs, norbs), dtype=np.complex128)
            hr_spin[:, 0:norbs:2, 0:norbs:2] = self.hr
            hr_spin[:, 1:norbs:2, 1:norbs:2] = self.hr
            return hr_spin
        # with spin, spin order: up up up dn dn dn
        elif ispin == 2:
            norbs = 2 * self.nwann
            hr_spin = np.zeros((self.nrpt, norbs, norbs), dtype=np.complex128)
            hr_spin[:, 0:self.nwann, 0:self.nwann] = self.hr
            hr_spin[:, self.nwann:norbs, self.nwann:norbs] = self.hr
            return hr_spin
        # without spin
        else:
            return self.hr[:, :, :]


class KVec():
    """
    Define :math:`k` points in BZ, high symmetry line or uniform grid.
    The coordinates are fractional w.r.t. the primitive reciprocal lattice vectors.
    """

    def __init__(self, kpt_type='uni', kbase=None, with_twopi=False, nkpt=None,
                 kvec=None, weights=None):
        """
        Parameters
        ----------
        kpt_type: string
            The type of :math:`k` points, 'uni' or 'sym'.
        kbase: a :math:`3 \\times 3` float array
            The basis vectors of the primitive reciprocal space.
        with_twopi: logical
            Whether the basis vector ``kbase`` including the :math:`2\\pi` factor.
        nkpt: int
            Number of :math:`k` points.
        kvec: 2d float array
            The fractional coordinates of :math:`k` points.
        weights: 1d float array
            The weights of :math:`k` points.
        """
        self.nkpt = nkpt
        self.kbase = np.array(kbase, dtype=np.float64)
        self.kvec = np.array(kvec, dtype=np.float64)
        self.weights = np.array(weights, dtype=np.float64)
        self.kpt_type = kpt_type
        self.with_twopi = with_twopi

    def set_base(self, kbase, with_twopi=False):
        """
        Set the basis of the primitive reciprocal.

        Parameters
        ----------
        kbase: a :math:`3 \\times 3` float array
            The basis with respect to the global axis.
        with_twopi: logical
            Whether the basis vector ``kbase`` including the :math:`2\\pi` factor.
        """

        self.kbase = np.array(kbase, dtype=np.float64)
        self.with_twopi = with_twopi

    def kvec_from_file(self, fname, read_weights=False):
        """
        Read the fractional coordinates and weights of :math:`k` points from a file.

        Parameters
        ----------
        fname: string
            File name.
        read_weights: logical
            Whether to read the weights.
        """

        k_tmp = []
        w_tmp = []
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip().split()
                if line != []:
                    k_tmp.append(line[0:3])
                    if read_weights:
                        w_tmp.append(line[3])
            self.kvec = np.array(k_tmp, dtype=np.float64)
            self.nkpt = len(k_tmp)
            if read_weights:
                self.weights = np.array(w_tmp, dtype=np.float64)
            else:
                self.weights = np.ones(self.nkpt) / self.nkpt

    def write_kvec(self, fname, write_weights=False):
        """
        Write the fractional coordinates of :math:`k` points to a file.

        Parameters
        ----------
        fname: string
            File name.
        write_weights: logical
            Whether to write the weights of each :math:`k` point.
        """
        with open(fname, 'w') as f:
            for i in range(len(self.kvec)):
                if write_weights:
                    fmt = "{:20.10f}"*4 + str("\n")
                    line = fmt.format(self.kvec[i, 0], self.kvec[i, 1], self.kvec[i, 2],
                                      self.weights[i])
                else:
                    fmt = "{:20.10f}"*3 + str("\n")
                    line = fmt.format(self.kvec[i, 0], self.kvec[i, 1], self.kvec[i, 2])
                f.write(line)


class SymKVec(KVec):
    """
    Class for defining :math:`k` points in high symmetry line, derived from :class:`KVec`.

    """

    def __init__(self, kbase=None, with_twopi=False, hsymkpt=None, klen=None):
        """
        Parameters
        ----------
        kbase: a :math:`3 \\times 3` float array
            Basis of the primitive reciprocal lattice.
        with_twopi: logical
            Whether the basis vector ``kbase`` including the :math:`2\\pi` factor.
        hsymkpt: float array
            Starting and end :math:`k` points along high symmetry lines.
        klen: float array
            Length of segments of :math:`k` points line.
        """
        self.klen = np.array(klen, dtype=np.float64)
        self.hsymkpt = np.array(hsymkpt, dtype=np.float64)
        KVec.__init__(self, kpt_type='sym', kbase=kbase, with_twopi=with_twopi)

    def get_klen(self):
        """
        Calculate the length of :math:`k` points segments.
        """

        self.klen = np.zeros(self.nkpt, dtype=np.float64)
        self.klen[0] = 0.0
        prev_kpt = self.kvec[0]

        for i in range(1, self.nkpt):
            curr_kpt = self.kvec[i, :]
            tmp_kpt = curr_kpt - prev_kpt
            kx = np.dot(tmp_kpt, self.kbase[:, 0])
            ky = np.dot(tmp_kpt, self.kbase[:, 1])
            kz = np.dot(tmp_kpt, self.kbase[:, 2])
            self.klen[i] = self.klen[i - 1] + np.sqrt(np.dot((kx, ky, kz), (kx, ky, kz)))
            prev_kpt = curr_kpt

    def from_hsymkpt_npt(self, nkpt_per_path=20):
        """
        Given the number of points per each segment,
        calculate the franctional coordinates of high symmetry :math:`k` points.

        Parameters
        ----------
        nkpt_per_path: int
            Number of :math:`k` points per each segment.
        """

        self.nkpt = nkpt_per_path * (len(self.hsymkpt) - 1)
        self.kvec = np.zeros((self.nkpt, 3), dtype=np.float64)
        self.weights = np.ones(self.nkpt) / self.nkpt
        for i in range(1, len(self.hsymkpt)):
            kpt_prev = self.hsymkpt[i - 1, :]
            kpt_curr = self.hsymkpt[i, :]
            for j in range(nkpt_per_path):
                ikpt = (i - 1) * nkpt_per_path + j
                self.kvec[ikpt, :] = (float(j) / float(nkpt_per_path - 1) *
                                      (kpt_curr - kpt_prev) + kpt_prev)

    def from_hsymkpt_step(self, step):
        """
        Given a step,
        calculate the franctional coordinates of high symmetry :math:`k` points.

        Parameters
        ----------
        step: float
            Step size.
        """

        kvec = []
        self.hsym_dis = np.zeros(len(self.hsymkpt), dtype=np.float64)
        self.hsym_dis[0] = 0.0
        for i in range(0, len(self.hsymkpt) - 1):
            kpt_prev = self.hsymkpt[i, :]
            kpt_curr = self.hsymkpt[i + 1, :]
            tmp = np.dot(self.kbase.transpose(), kpt_curr - kpt_prev)
            dis = np.sqrt(np.dot(tmp, tmp))
            self.hsym_dis[i + 1] = self.hsym_dis[i] + dis
            pts = np.arange(0, dis, step) / dis
            for ipt in pts:
                kvec.append(ipt * (kpt_curr - kpt_prev) + kpt_prev)
        self.kvec = np.array(kvec, dtype=np.float64)
        self.nkpt = len(self.kvec)
        self.weights = np.ones(self.nkpt) / self.nkpt


class UniKVec(KVec):
    """
    Class for defining uniform :math:`k` points grid, derived from :class:`KVec`.
    """

    def __init__(self, kbase=None, with_twopi=False, grid=None):
        """
        Parameters
        ----------
        kbase: a :math:`3 \\times 3` float array
            Basis of the primitive reciprocal lattice.
        with_twopi: logical
            Whether the basis vector ``kbase`` including the :math:`2\\pi` factor.
        grid: 3-elements tuple
            Three numbers defining a uniform grid, for example: :math:`11 \\times 11 \\times 11`.
        """
        self.grid = grid
        KVec.__init__(self, kpt_type='uni', kbase=kbase, with_twopi=with_twopi)

    def from_grid(self, shift_delta=0.0):
        """
        Return uniform :math:`k` points.

        Parameters
        ----------
        shift_delta: float
            A small shift.
        """

        nx, ny, nz = self.grid
        self.nkpt = nx * ny * nz
        self.kvec = np.zeros((self.nkpt, 3), dtype=np.float64)
        self.weights = np.ones(self.nkpt) / self.nkpt
        ikpt = 0
        for i in range(nx):
            if nx == 1:
                kx = 0.0
            else:
                kx = float(i) / float(nx)
            for j in range(ny):
                if ny == 1:
                    ky = 0.0
                else:
                    ky = float(j) / float(ny)
                for k in range(nz):
                    if nz == 1:
                        kz = 0.0
                    else:
                        kz = float(k) / float(nz)
                    ikpt = ikpt + 1
                    self.kvec[ikpt-1, :] = kx + shift_delta, ky + shift_delta, kz + shift_delta
