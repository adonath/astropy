# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A set of standard astronomical equivalencies.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..constants import si as _si
from . import si
from . import cgs
from . import astrophys

__all__ = ['parallax', 'spectral', 'spectral_density', 'doppler_radio',
           'doppler_optical', 'doppler_relativistic', 'mass_energy', 'brightness_temperature']


def parallax():
    """
    Returns a list of equivalence pairs that handle the conversion
    between parallax angle and distance.
    """
    return [
        (si.arcsecond, astrophys.parsec, lambda x: 1. / x)
    ]


def spectral():
    """
    Returns a list of equivalence pairs that handle spectral
    wavelength, frequency, and energy equivalences.

    Allows conversions between wavelength units, frequency units and
    energy units as they relate to light.
    """

    return [
        (si.m, si.Hz, lambda x: _si.c.value / x),
        (si.m, si.J, lambda x: (_si.c.value * _si.h.value) / x),
        (si.Hz, si.J, lambda x: _si.h.value * x)
    ]


def spectral_density(sunit, sfactor):
    """
    Returns a list of equivalence pairs that handle spectral density
    with regard to wavelength and frequency.
    """
    c_Aps = _si.c.value * 10 ** 10

    fla = cgs.erg / si.angstrom / si.cm ** 2 / si.s
    fnu = cgs.erg / si.Hz / si.cm ** 2 / si.s
    nufnu = cgs.erg / si.cm ** 2 / si.s
    lafla = nufnu

    def converter(x):
        return x * (sunit.to(si.AA, sfactor, spectral()) ** 2 / c_Aps)

    def iconverter(x):
        return x / (sunit.to(si.AA, sfactor, spectral()) ** 2 / c_Aps)

    def converter_fnu_nufnu(x):
        return x * sunit.to(si.Hz, sfactor, spectral())

    def iconverter_fnu_nufnu(x):
        return x / sunit.to(si.Hz, sfactor, spectral())

    def converter_fla_lafla(x):
        return x * sunit.to(si.AA, sfactor, spectral())

    def iconverter_fla_lafla(x):
        return x / sunit.to(si.AA, sfactor, spectral())

    return [
        (si.AA, fnu, converter, iconverter),
        (fla, fnu, converter, iconverter),
        (si.AA, si.Hz, converter, iconverter),
        (fla, si.Hz, converter, iconverter),
        (fnu, nufnu, converter_fnu_nufnu, iconverter_fnu_nufnu),
        (fla, lafla, converter_fla_lafla, iconverter_fla_lafla),
    ]


def doppler_radio(rest):
    r"""
    Return the equivalency pairs for the radio convention for velocity.

    The radio convention for the relation between velocity and frequency is:
    
    :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    References
    ----------
    `NRAO site defining the conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> radio_CO_equiv = u.doppler_radio(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> radio_velocity = measured_freq.to(u.km/u.s, equivalencies=radio_CO_equiv)
    >>> radio_velocity
    <Quantity -31.2090920889... km / s>
    """

    ckms = _si.c.to('km/s').value

    def to_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        return (restfreq-x) / (restfreq) * ckms

    def from_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        voverc = x/ckms
        return restfreq * (1-voverc)


    def to_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return (x-restwav) / (x) * ckms

    def from_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return restwav * ckms / (ckms-x)


    def to_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        return (resten-x) / (resten) * ckms

    def from_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        voverc = x/ckms
        return resten * (1-voverc)

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]


def doppler_optical(rest):
    r"""
    Return the equivalency pairs for the optical convention for velocity.

    The optical convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + V/c )^{-1}`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    References
    ----------
    `NRAO site defining the conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> optical_CO_equiv = u.doppler_optical(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> optical_velocity = measured_freq.to(u.km/u.s, equivalencies=optical_CO_equiv)
    >>> optical_velocity
    <Quantity -31.205843488... km / s>
    """

    ckms = _si.c.to('km/s').value

    def to_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        return ckms * (restfreq-x) / x

    def from_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        voverc = x/ckms
        return restfreq / (1+voverc)


    def to_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return ckms * (x/restwav-1)

    def from_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        voverc = x/ckms
        return restwav * (1+voverc)


    def to_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        return ckms * (resten-x) / x

    def from_vel_en(x):
        resten = rest.to(si.eV, equivalencies=spectral()).value
        voverc = x/ckms
        return resten / (1+voverc)

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]


def doppler_relativistic(rest):
    r"""
    Return the equivalency pairs for the relativistic convention for velocity.

    The full relativistic convention for the relation between velocity and frequency is:

    :math:`V = c \frac{f_0^2 - f^2}{f_0^2 + f^2} ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

    Parameters
    ----------
    rest : Quantity
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency)

    References
    ----------
    `NRAO site defining the conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_

    Examples
    --------
    >>> import astropy.units as u
    >>> CO_restfreq = 115.27120*u.GHz  # rest frequency of 12 CO 1-0 in GHz
    >>> relativistic_CO_equiv = u.doppler_relativistic(CO_restfreq)
    >>> measured_freq = 115.2832*u.GHz
    >>> relativistic_velocity = measured_freq.to(u.km/u.s, equivalencies=relativistic_CO_equiv)
    >>> relativistic_velocity
    <Quantity -31.207467619...
    >>> measured_velocity = 1250 * u.km/u.s
    >>> relativistic_frequency = measured_velocity.to(u.GHz, equivalencies=relativistic_CO_equiv)
    >>> relativistic_frequency
    <Quantity 114.7915686...
    >>> relativistic_wavelength = measured_velocity.to(u.mm, equivalencies=relativistic_CO_equiv)
    >>> relativistic_wavelength
    <Quantity 2.6116243681...
    """

    ckms = _si.c.to('km/s').value

    def to_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        return (restfreq**2-x**2) / (restfreq**2+x**2) * ckms

    def from_vel_freq(x):
        restfreq = rest.to(si.Hz, equivalencies=spectral()).value
        voverc = x/ckms
        return restfreq * ((1-voverc) / (1+(voverc)))**0.5


    def to_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        return (x**2-restwav**2) / (restwav**2+x**2) * ckms

    def from_vel_wav(x):
        restwav = rest.to(si.AA, spectral()).value
        voverc = x/ckms
        return restwav * ((1+voverc) / (1-voverc))**0.5


    def to_vel_en(x):
        resten = rest.to(si.eV, spectral()).value
        return (resten**2-x**2) / (resten**2+x**2) * ckms

    def from_vel_en(x):
        resten = rest.to(si.eV, spectral()).value
        voverc = x/ckms
        return resten * ((1-voverc) / (1+(voverc)))**0.5

    return [(si.Hz, si.km/si.s, to_vel_freq, from_vel_freq),
            (si.AA, si.km/si.s, to_vel_wav, from_vel_wav),
            (si.eV, si.km/si.s, to_vel_en, from_vel_en),
            ]


def mass_energy():
    """
    Returns a list of equivalence pairs that handle the conversion
    between mass and energy.
    """

    return [(si.kg, si.J, lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.m ** 2, si.J / si.m ** 2 , 
             lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.m ** 3, si.J / si.m ** 3 , 
             lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
            (si.kg / si.s, si.J / si.s , lambda x: x * _si.c.value ** 2,
             lambda x: x / _si.c.value ** 2),
    ]

def brightness_temperature(beam_area, disp):
    """
    "Antenna Gain" or "sensitivity" equivalency: Defines the conversion between
    Jy/beam and "brightness temperature", :math:`T_B`, in Kelvins.  This is a
    unit very commonly used in radio astronomy.  Typically, the gain refers to
    the conversion between corrected antenna temperature :math:`T_A^*` and flux
    density.  See, e.g., "Tools of Radio Astronomy" (Wilson 2009) eqn 8.16 and
    eqn 8.19 (these pages are available on `google books
    <http://books.google.com/books?id=9KHw6R8rQEMC&pg=PA179&source=gbs_toc_r&cad=4#v=onepage&q&f=false>`__).

    :math:`T_B \equiv S_\\nu / \left(2 k \\nu^2 / c^2 \\right)`

    However, the beam area is essential for this computation: the brighntess
    temperature is inversely proportional to the beam area

    Parameters
    ----------
    beam_area : Beam Area equivalent
        Beam area in angular units, i.e. steradian equivalent
    disp : `Quantity` with spectral units
        The observed `spectral` equivalent `Unit` (e.g., frequency or
        wavelength)

    Examples
    --------
    Arecibo C-band beam gain ~ 7 K/Jy::

        >>> import numpy as np
        >>> from astropy import units as u
        >>> beam_area = np.pi*(50*u.arcsec)**2
        >>> freq = 5*u.GHz
        >>> u.Jy.to(u.K, equivalencies=u.brightness_temperature(beam_area,freq))
        7.052588858...
        >>> (1*u.Jy).to(u.K, equivalencies=u.brightness_temperature(beam_area,freq))
        <Quantity 7.05258...

    VLA synthetic beam::

        >>> beam_area = np.pi*(15*u.arcsec)**2
        >>> freq = 5*u.GHz
        >>> u.Jy.to(u.K, equivalencies=u.brightness_temperature(beam_area,freq))
        78.36209843...
    """
    beam = beam_area.to(si.sr).value
    nu = disp.to(si.GHz, spectral())

    def convert_Jy_to_K(x_jybm):
        factor = (2 * _si.k_B * si.K * nu**2 / _si.c**2).to(astrophys.Jy).value
        return (x_jybm / beam / factor)

    def convert_K_to_Jy(x_K):
        factor = (astrophys.Jy / (2 * _si.k_B * nu**2 / _si.c**2)).to(si.K).value
        return (x_K * beam / factor)

    return [(astrophys.Jy, si.K, convert_Jy_to_K, convert_K_to_Jy)]
