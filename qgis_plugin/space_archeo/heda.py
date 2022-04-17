# -*- coding: utf-8 -*-
"""

    ---------------------
    Date                 : 01.2021
    Copyright            : (C) 2021 by Rad.Data Spectral Analytics UG
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************

The original spectral EDGE Filtering techique is described here:
Bakker et al. 2002
https://doi.org/10.1016/S0924-2716(02)00060-6

Its based sensor PSF adoption is described here:
Mielke et al. 2016:
https://doi.org/10.3390/rs8020127
"""

import copy

import numpy
from scipy import ndimage
from scipy import signal


class MovingWindow:
    def __init__(self, image, winsize):
        if winsize % 2 == 0:
            winsize = winsize - 1
        self.image = image
        self.windowsize = (winsize // 2, winsize // 2)

    def untermat(self, i, j):
        self.fensterchen = self.image[:, i - self.windowsize[1]:i + self.windowsize[1] + 1,
                           j - self.windowsize[1]:j + self.windowsize[1] + 1]
        self.mitte = self.fensterchen[:, self.windowsize[0], self.windowsize[0]]
        self.fsh = self.fensterchen.shape


class MovingWindow2d:
    def __init__(self, image, winsize):
        if winsize % 2 == 0:
            winsize = winsize - 1
        self.image = image
        self.windowsize = (winsize // 2, winsize // 2)

    def untermat(self, i, j):
        self.fensterchen = self.image[i - self.windowsize[1]:i + self.windowsize[1] + 1,
                           j - self.windowsize[1]:j + self.windowsize[1] + 1]
        self.mitte = self.fensterchen[self.windowsize[0], self.windowsize[0]]
        self.fsh = self.fensterchen.shape


class Filters:
    @staticmethod
    def medianize_im(gscimage: numpy.ndarray, kernel: int = 3) -> numpy.ndarray:
        """
        Median Filtering of the result matrices
        :param gscimage: 2d numpy array
        :param kernel: filter kernel
        :return: filtered result
        """
        raus = copy.deepcopy(gscimage)
        result = signal.medfilt2d(raus, kernel_size=kernel)
        return result

    @staticmethod
    def mediannize_all(im: numpy.ndarray) -> numpy.ndarray:
        """
        Median Filtering of the result matrices
        :param gscimage: 2d numpy array
        :return: filtered result
        """
        lsh = im.shape
        L = []
        for i in numpy.arange(0, lsh[0], 1):
            ll = Filters.medianize_im(im[i, :, :])
            L.append(ll)
        return numpy.asarray(L)

    @staticmethod
    def define_LOG(size: int, sensor: str, sigma: float = 0.55248) -> numpy.ndarray:
        """
        Define a LoG Filter for hyperspectral Edge Filtering.
        :param size: Size of the filter; must be an odd integer
        :param sensor: name the sensor (prisma,hymap,hyperion) that you'd like to use to construct the PSF
        :param sigma: define your own sigma to construct the PSF
        :return: Kernel for edge filtering
        """
        if size % 2 == 0:
            print("Error Expected odd Integer as Size Input!")
        Kernel = numpy.zeros([size, size])
        if sensor == 'hyperion':
            sigma = 0.76433
        elif sensor == 'prisma':
            sigma = 0.52
        elif sensor == 'hymap':
            sigma = 0.5
        for i in numpy.arange((0 - size // 2), (0 + size // 2) + 1, 1):
            for j in numpy.arange((0 - size // 2), (0 + size // 2) + 1, 1):
                w = ((numpy.power(i, 2) + numpy.power(j, 2)) / (2 * numpy.power(sigma, 2)))
                Kernel[i, j] = -(1 / (numpy.pi * numpy.power(sigma, 4))) * (1 - w) * numpy.exp(-w)
        Kernel = numpy.roll(Kernel, -size // 2, axis=0)
        Kernel = numpy.roll(Kernel, -size // 2, axis=1)
        w = numpy.sum(Kernel)
        Kernel /= w
        return Kernel

    @staticmethod
    def define_LOG_generic(size, sigma):
        """
        Define a generic LoG Filter for hyperspectral Edge Filtering.
        :param size: Size of the filter; must be an odd integer
        :param sigma: define your own sigma to construct the PSF
        :return: Kernel for edge filtering
        """
        if size % 2 == 0:
            print("Error Expected odd Integer as Size Input!")
        Kernel = numpy.zeros([size, size])
        for i in numpy.arange((0 - size // 2), (0 + size // 2) + 1, 1):
            for j in numpy.arange((0 - size // 2), (0 + size // 2) + 1, 1):
                w = ((numpy.power(i, 2) + numpy.power(j, 2)) / (2 * numpy.power(sigma, 2)))
                Kernel[i, j] = -(1 / (numpy.pi * numpy.power(sigma, 4))) * (1 - w) * numpy.exp(-w)
        Kernel = numpy.roll(Kernel, -size // 2, axis=0)
        Kernel = numpy.roll(Kernel, -size // 2, axis=1)
        w = numpy.sum(Kernel)
        Kernel /= w
        return Kernel

    @staticmethod
    def vcorrcoef(M: numpy.ndarray, n: numpy.ndarray) -> numpy.ndarray:
        """
        Method for calculation of the corrcoeffs of the centerpoint to its neighbours
        :param M: whole matrix around center point
        :param n: centerpoint
        :return: Corrmatrix
        """
        Mm = numpy.reshape(numpy.mean(M, axis=1), (M.shape[0], 1))
        nm = numpy.mean(n)
        r_zahl = numpy.sum((M - Mm) * (n - nm), axis=1)
        r_nen = numpy.sqrt(numpy.sum((M - Mm) ** 2, axis=1) * numpy.sum((n - nm) ** 2))
        r = r_zahl / r_nen
        r[~numpy.isfinite(r)] = 0
        return r

    @staticmethod
    def bakker_filter_log_MSAM(image: int, sensor: str, winsize: int = 3):
        """
        Bakker Detector that calculates the "Correlation Edge"
        :param image: 3d-Data Cube
        :param sensor: sensor (hyperion,enmap,hymap) defines the PSF
        :param winsize: odd integer number for defining the size of the correlation window.
        :return: Output Result
        """
        if winsize % 2 == 0:
            winsize = winsize - 1
        lsh = image.shape
        slidwin = MovingWindow(image, winsize)
        if sensor != None:
            Kantenfilter = Filters.define_LOG(winsize, sensor)
        else:
            Kantenfilter = Filters.define_LOG(winsize, sensor='hyperion')
        outimage = numpy.zeros([lsh[1], lsh[2]])
        for i in numpy.arange(0 + winsize // 2, lsh[1] - winsize // 2, 1):
            for j in numpy.arange(0 + winsize // 2, lsh[2] - winsize // 2, 1):
                slidwin.untermat(i, j)
                fensterchen = slidwin.fensterchen
                alles = fensterchen.reshape(slidwin.fsh[0], slidwin.fsh[1] * slidwin.fsh[2])
                samall = Filters.vcorrcoef(alles.T, slidwin.mitte)
                samall = samall.reshape(slidwin.fsh[1], slidwin.fsh[2])
                mit_zentrum = numpy.dot(Kantenfilter, samall)
                ohne_zentrum = copy.deepcopy(mit_zentrum)
                ohne_zentrum[winsize // 2, winsize // 2] = 0
                outimage[i, j] = numpy.sum(ohne_zentrum)
        return outimage


class EdgeFiltering:
    """
    Class with functions to carry out the Edge Filtering
    """

    @staticmethod
    def combined_ops_bakker_MSAM_cubo_finale(cubo: numpy.ndarray, sensor: str, winsize: int = 5):
        """

        :param cubo: Numpy nd-array of the data that needs to be filtered
        :param sensor: (hymap,enmap,hyperion)
        :param winsize: odd integer that defines the corrwindowsize
        :return: LoG filtered cube
        """
        hyper_edge = Filters.bakker_filter_log_MSAM(cubo, sensor, winsize)
        return hyper_edge

    @staticmethod
    def heda_archstyle(cubo, sensor, winsize=3):
        """
        :param cubo: Numpy nd-array of the data that needs to be filtered
        :param sensor: (hymap,enmap,hyperion)
        :param winsize: odd integer that defines the corrwindowsize
        :return: LoG filtered cube
        """
        gradient = EdgeFiltering.combined_ops_bakker_MSAM_cubo_finale(cubo, sensor, winsize)
        abs_grad = numpy.abs(gradient)
        absgradsmooth = ndimage.filters.median_filter(abs_grad, size=winsize + 2)
        gradsmooth = ndimage.filters.median_filter(gradient, size=winsize + 2)
        return gradient, abs_grad, absgradsmooth, gradsmooth
