# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h" nogil:
    bint npy_isnan(double x)

cimport cython


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interpolate1d(np.ndarray[DTYPE_t, ndim=1] f,
                                np.ndarray[DTYPE_t, ndim=1] g):
    cdef int nx = f.shape[0]
    cdef int nkx = g.shape[0]
    cdef int wkx = nkx // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=1] fixed = f.copy()

    cdef unsigned int i, ii

    cdef int iimin, iimax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:
        # Need a first pass to replace NaN values with value convolved from
        # neighboring values
        for i in range(nx):
            if npy_isnan(f[i]) and i >= wkx and i < nx - wkx:
                top = 0.
                bot = 0.
                for ii in range(i - wkx, i + wkx + 1):
                    val = f[ii]
                    if not npy_isnan(val):
                        ker = g[<unsigned int>(wkx + ii - i)]
                        top += val * ker
                        bot += ker
                if bot != 0.:
                    fixed[i] = top / bot
                else:
                    fixed[i] = f[i]
            else:
                fixed[i] = f[i]
<<<<<<< HEAD:astropy/convolution/boundary_none.pyx

        # Now run the proper convolution
        for i in range(wkx, nx - wkx):
            if not npy_isnan(fixed[i]):
                top = 0.
                bot = 0.
                for ii in range(i - wkx, i + wkx + 1):
                    val = fixed[ii]
                    ker = g[<unsigned int>(wkx + ii - i)]
                    if not npy_isnan(val):
                        top += val * ker
                        bot += ker
                if bot != 0:
                    conv[i] = top / bot
                else:
                    conv[i] = fixed[i]
            else:
                conv[i] = fixed[i]
    # GIL acquired again here
=======
        else:
            fixed[i] = f[i]
    return fixed


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve1d(np.ndarray[DTYPE_t, ndim=1] f,
                             np.ndarray[DTYPE_t, ndim=1] g):
    cdef int nx = f.shape[0]
    cdef int nkx = g.shape[0]
    cdef int wkx = nkx // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=1] conv = np.zeros([nx], dtype=DTYPE)

    cdef unsigned int i, ii

    cdef int iimin, iimax

    cdef DTYPE_t top, bot, ker, val

    # Now run the proper convolution
    for i in range(wkx, nx - wkx):
        if not npy_isnan(f[i]):
            top = 0.
            bot = 0.
            for ii in range(i - wkx, i + wkx + 1):
                val = f[ii]
                ker = g[<unsigned int>(wkx + ii - i)]
                if not npy_isnan(val):
                    top += val * ker
                    bot += ker
            if bot != 0:
                conv[i] = top / bot
            else:
                conv[i] = f[i]
        else:
            conv[i] = f[i]
>>>>>>> Use np.pad for boundary handling:astropy/convolution/convolve_cython.pyx
    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interpolate2d(np.ndarray[DTYPE_t, ndim=2] f,
                                np.ndarray[DTYPE_t, ndim=2] g):
    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=2] fixed = f.copy()

    cdef unsigned int i, j, ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:
        # Need a first pass to replace NaN values with value convolved from
        # neighboring values
        for i in range(nx):
            for j in range(ny):
                if npy_isnan(f[i, j]) and i >= wkx and i < nx - wkx \
                and j >= wky and j < ny - wky:
                    top = 0.
                    bot = 0.
                    for ii in range(i - wkx, i + wkx + 1):
                        for jj in range(j - wky, j + wky + 1):
                            val = f[ii, jj]
                            if not npy_isnan(val):
                                ker = g[<unsigned int>(wkx + ii - i),
                                        <unsigned int>(wky + jj - j)]
                                top += val * ker
                                bot += ker
                    if bot != 0.:
                        fixed[i, j] = top / bot
                    else:
                        fixed[i, j] = f[i, j]
                else:
                    fixed[i, j] = f[i, j]
<<<<<<< HEAD

        # Now run the proper convolution
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                if not npy_isnan(fixed[i, j]):
                    top = 0.
                    bot = 0.
                    for ii in range(i - wkx, i + wkx + 1):
                        for jj in range(j - wky, j + wky + 1):
                            val = fixed[ii, jj]
                            ker = g[<unsigned int>(wkx + ii - i),
                                    <unsigned int>(wky + jj - j)]
                            if not npy_isnan(val):
                                top += val * ker
                                bot += ker
                    if bot != 0:
                        conv[i, j] = top / bot
                    else:
                        conv[i, j] = fixed[i, j]
                else:
                    conv[i, j] = fixed[i, j]
    # GIL acquired again here
=======
            else:
                continue
    return fixed


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve2d(np.ndarray[DTYPE_t, ndim=2] f,
                             np.ndarray[DTYPE_t, ndim=2] g):
    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.zeros([nx, ny], dtype=DTYPE)

    cdef unsigned int i, j, ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    # Now run the proper convolution
    for i in range(wkx, nx - wkx):
        for j in range(wky, ny - wky):
            if not npy_isnan(f[i, j]):
                top = 0.
                bot = 0.
                for ii in range(i - wkx, i + wkx + 1):
                    for jj in range(j - wky, j + wky + 1):
                        val = f[ii, jj]
                        ker = g[<unsigned int> (wkx + ii - i),
                                <unsigned int> (wky + jj - j)]
                        if not npy_isnan(val):
                            top += val * ker
                            bot += ker
                if bot != 0:
                    conv[i, j] = top / bot
                else:
                    conv[i, j] = f[i, j]
            else:
                conv[i, j] = f[i, j]

    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve2d_symmetric(np.ndarray[DTYPE_t, ndim=2] f,
                                       np.ndarray[DTYPE_t, ndim=2] g):
    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.zeros([nx, ny], dtype=DTYPE)

    cdef unsigned int i, j, ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    if wkx == 0:
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                val = f[i, j] * g[0, wky]
                for jj in range(1, wky + 1):
                    val += (f[i, j + jj] +
                            f[i, j - jj]) * g[0, jj + wky]
                conv[i, j] = val
    elif wky == 0:
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                val = f[i, j] * g[wkx, 0]
                for ii in range(1, wkx + 1):
                    val += (f[i + ii, j] +
                            f[i - ii, j]) * g[ii + wkx, 0]
                conv[i, j] = val
    else:
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                val = f[i, j] * g[wkx, wky]
                for ii in range(1, wkx + 1):
                    val += (f[i + ii, j] +
                            f[i - ii, j] +
                            f[i, j - ii] +
                            f[i, j + ii]) * g[ii + wkx, wkx]
                for ii in range(1, wkx + 1):
                    for jj in range(1, wky + 1):
                        val += (f[i + ii, j + jj] +
                                f[i - ii, j + jj] +
                                f[i + ii, j - jj] +
                                f[i - ii, j - jj]) * g[ii + wkx, jj + wky]
                conv[i, j] = val
    return conv


def convolve2d_separable():
    raise NotImplementedError


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve2d_symmetric_cache(np.ndarray[DTYPE_t, ndim=2] f,
                               np.ndarray[DTYPE_t, ndim=2] g):
    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.zeros([nx, ny], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] tmp = np.zeros([wkx + 1, ny], dtype=DTYPE)

    cdef unsigned int i, j, ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    for i in range(wkx, nx - wkx):
        # Setup tmp array
        for j in range(ny):
            for ii in range(1, wkx + 1):
                tmp[ii, j] = f[i + ii, j] + f[i - ii, j]
        for j in range(wky, ny - wky):
            val = f[i, j] * g[wkx, wky]
            for ii in range(1, wkx + 1):
                val += (tmp[ii, j] + f[i, j - ii] + f[i, j + ii]) * g[ii + wkx, wkx]
            for ii in range(1, wkx + 1):
                for jj in range(1, wky + 1):
                    val += (tmp[ii, j + jj] + tmp[ii, j - jj]) * g[ii + wkx, jj + wky]
            conv[i, j] = val
    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
<<<<<<< HEAD:astropy/convolution/boundary_none.pyx
def convolve2d_boundary_none_symmetric(np.ndarray[DTYPE_t, ndim=2] f,
                             np.ndarray[DTYPE_t, ndim=2] g):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.zeros([nx, ny], dtype=DTYPE)

    cdef unsigned int i, j, ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val
    
    if wkx == 0:
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                val = f[i, j] * g[0, wky]
                for jj in range(1, wky + 1):
                    val += (f[i, j + jj] +
                            f[i, j - jj]) * g[0, jj + wky]
                conv[i, j] = val
    elif wky == 0:
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                val = f[i, j] * g[wkx, 0]
                for ii in range(1, wkx + 1):
                    val += (f[i + ii, j] +
                            f[i - ii, j]) * g[ii + wkx, 0]
                conv[i, j] = val
    else:
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                val = f[i, j] * g[wkx, wky]
                for ii in range(1, wkx + 1):
                    for jj in range(0, wky + 1):
                        val += (f[i + ii, j + jj] +
                                f[i - jj, j + ii] +
                                f[i + jj, j - ii] +
                                f[i - ii, j - jj]) * g[ii + wkx, jj + wky]
                conv[i, j] = val
>>>>>>> Added separable convolution and code cleanup
    return conv



@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve3d_boundary_none(np.ndarray[DTYPE_t, ndim=3] f,
=======
def convolve3d(np.ndarray[DTYPE_t, ndim=3] f,
>>>>>>> Use np.pad for boundary handling:astropy/convolution/convolve_cython.pyx
                             np.ndarray[DTYPE_t, ndim=3] g):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1 or g.shape[2] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nz = f.shape[2]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int nkz = g.shape[2]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2
    cdef int wkz = nkz // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=3] fixed = np.zeros([nx, ny, nz], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=3] conv = np.zeros([nx, ny, nz], dtype=DTYPE)

    cdef unsigned int i, j, k, ii, jj, kk

    cdef int iimin, iimax, jjmin, jjmax, kkmin, kkmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:
        # Need a first pass to replace NaN values with value convolved from
        # neighboring values
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if npy_isnan(f[i, j, k]) and i >= wkx and i < nx - wkx \
                    and j >= wky and j < ny - wky and k >= wkz and k <= nz - wkz:
                        top = 0.
                        bot = 0.
                        for ii in range(i - wkx, i + wkx + 1):
                            for jj in range(j - wky, j + wky + 1):
                                for kk in range(k - wkz, k + wkz + 1):
                                    val = f[ii, jj, kk]
                                    if not npy_isnan(val):
                                        ker = g[<unsigned int>(wkx + ii - i),
                                                <unsigned int>(wky + jj - j),
                                                <unsigned int>(wkz + kk - k)]
                                        top += val * ker
                                        bot += ker
                        if bot != 0.:
                            fixed[i, j, k] = top / bot
                        else:
                            fixed[i, j, k] = f[i, j, k]
                    else:
                        fixed[i, j, k] = f[i, j, k]

        # Now run the proper convolution
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                for k in range(wkz, nz - wkz):
                    if not npy_isnan(fixed[i, j, k]):
                        top = 0.
                        bot = 0.
                        for ii in range(i - wkx, i + wkx + 1):
                            for jj in range(j - wky, j + wky + 1):
                                for kk in range(k - wkz, k + wkz + 1):
                                    val = fixed[ii, jj, kk]
                                    ker = g[<unsigned int>(wkx + ii - i),
                                            <unsigned int>(wky + jj - j),
                                            <unsigned int>(wkz + kk - k)]
                                    if not npy_isnan(val):
                                        top += val * ker
                                        bot += ker
                        if bot != 0:
                            conv[i, j, k] = top / bot
                        else:
                            conv[i, j, k] = fixed[i, j, k]
                    else:
                        conv[i, j, k] = fixed[i, j, k]
    # GIL acquired again here
    return conv


def interpolate3d():
    raise NotImplementedError


ndim_1 = {'x': {True: convolve1d,
                False: convolve1d},
          'none': {True: convolve1d,
                   False: convolve1d}}

ndim_2 = {'radial': {True: convolve2d_symmetric_cache,
                     False: convolve2d_symmetric_cache},
          'xy': {True: convolve2d_symmetric_cache,
                 False: convolve2d_symmetric_cache},
          'none': {True: convolve2d,
                   False: convolve2d}}

ndim_3 = {'none': {True: convolve3d,
                   False: convolve3d}}

convolve_functions = {1: ndim_1, 2: ndim_2, 3: ndim_3}

interpolate_functions = {1: interpolate1d,
                         2: interpolate2d,
                         3: interpolate3d}
