__all__ = ['get_spectra_from_poles', 'plot_spectrum', 'plot_rixs_map']

import numpy as np
import matplotlib.pyplot as plt
from .utils import boltz_dist
from .iostream import read_poles_from_file


def get_spectra_from_poles(poles_dict, omega_mesh, gamma_mesh, temperature):
    """
    Given dict of poles, calculate XAS or RIXS spectra.

    poles_dict: dict
        Dict containing information of poles, which are calculated from
        xas_fsolver and rixs_fsolver.
        This dict is constructed by :func:`iostream.read_poles_from_file`.
    omega_mesh: 1d float array
        Energy grid.
    gamma_mesh: 1d float array
        Life-time broadening.
    temperature: float number
        Temperature (K) for boltzmann distribution.

    Returns
    -------
    spectra: 1d float array
        The calculated XAS or RIXS spectra.

    See also
    --------
    iostream.read_poles_from_file: read XAS or RIXS poles files.
    """
    nom = len(omega_mesh)
    spectra = np.zeros(nom, dtype=np.float64)
    gs_dist = boltz_dist(poles_dict['eigval'], temperature)
    ngs = len(poles_dict['eigval'])
    for i in range(ngs):
        tmp_vec = np.zeros(nom, dtype=np.complex)
        neff = poles_dict['npoles'][i]
        alpha = poles_dict['alpha'][i]
        beta = poles_dict['beta'][i]
        eigval = poles_dict['eigval'][i]
        norm = poles_dict['norm'][i]
        for j in range(neff - 1, 0, -1):
            tmp_vec = beta[j - 1]**2 / (omega_mesh + 1j * gamma_mesh + eigval
                                        - alpha[j] - tmp_vec)
        tmp_vec = 1.0 / (omega_mesh + 1j * gamma_mesh +
                         eigval - alpha[0] - tmp_vec)
        spectra[:] += -1.0 / np.pi * np.imag(tmp_vec) * norm * gs_dist[i]

    return spectra


def plot_spectrum(file_list, omega_mesh, gamma_mesh, T=1.0, fname='spectrum.dat',
                  om_shift=0.0, fmt_float='{:.15f}'):
    """
    Reading poles :math:`\\alpha` and :math:`\\beta`, and calculate
    the spectrum using continued fraction formula,

    .. math::
        I(\\omega_{i}) =-\\frac{1}{\\pi}\\text{Im} \\left[ \\frac{1}{x - \\alpha_{0} -
        \\frac{\\beta_{1}^2}{x-\\alpha_{1} - \\frac{\\beta_{2}^2}{x-\\alpha_{2} - ...}} }\\right],

    where, :math:`x = \\omega_{i}+i\\Gamma_{i} + E_{g}`.

    Parameters
    ----------
    file_list: list of string
        Name of poles file.
    omega_mesh: 1d float array
        The frequency mesh.
    gamma_mesh: 1d float array
        The broadening factor, in general, it is frequency dependent.
    T: float (default: 1.0K)
        Temperature (K).
    fname: str (default: 'spectrum.dat')
        File name to store spectrum.
    om_shift: float (default: 0.0)
        Energy shift.
    fmt_float: str (default: '{:.15f}')
        Format for printing float numbers.
    """

    pole_dict = read_poles_from_file(file_list)
    spectrum = get_spectra_from_poles(pole_dict, omega_mesh, gamma_mesh, T)

    space = "  "
    fmt_string = (fmt_float + space) * 2 + '\n'
    f = open(fname, 'w')
    for i in range(len(omega_mesh)):
        f.write(fmt_string.format(omega_mesh[i] + om_shift, spectrum[i]))
    f.close()


def plot_rixs_map(rixs_data, ominc_mesh, eloss_mesh, fname='rixs_map.pdf'):
    """
    Given 2d RIXS data, plot a RIXS map and save it to a pdf file.

    Parameters
    ----------
    rixs_data: 2d float array
        Calculated RIXS data as a function of incident energy and energy loss.
    ominc_mesh: 1d float array
        Incident energy mesh.
    eloss_mesh: 1d float array
        Energy loss mesh.
    fname: string
        File name to save RIXS map.
    """

    fig, ax = plt.subplots()
    m, n = np.array(rixs_data).shape
    if len(ominc_mesh) == m and len(eloss_mesh) == n:
        plt.imshow(
            rixs_data, extent=[min(eloss_mesh), max(eloss_mesh), min(ominc_mesh), max(ominc_mesh)],
            origin='lower', aspect='auto', cmap='rainbow', interpolation='gaussian')
        plt.xlabel(r'Energy loss (eV)')
        plt.ylabel(r'Energy of incident photon (eV)')
    elif len(eloss_mesh) == m and len(ominc_mesh) == n:
        plt.imshow(
            rixs_data, extent=[min(ominc_mesh), max(ominc_mesh), min(eloss_mesh), max(eloss_mesh)],
            origin='lower', aspect='auto', cmap='rainbow', interpolation='gaussian')
        plt.ylabel(r'Energy loss (eV)')
        plt.xlabel(r'Energy of incident photon (eV)')
    else:
        raise Exception("Dimension of rixs_data is not consistent with ominc_mesh or eloss_mesh")
    plt.savefig(fname)
