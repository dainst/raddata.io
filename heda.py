# -*- coding: utf-8 -*-
"""
Created on Sun Feb 15 12:22:01 2015

Still experimental need to know if this is of interest or not for final implementation
"""

import numpy, copy
from scipy import ndimage
from scipy import signal



class moving_window(object):
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


class moving_window2d(object):
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


# Medianfilterung Bild
def medianize_im(gscimage,
                 kernel=3):
    raus = copy.deepcopy(gscimage)
    result = signal.medfilt2d(raus, kernel_size=kernel)
    return result


def mediannize_all(im, kernel=3):
    lsh = im.shape
    L = []
    for i in numpy.arange(0, lsh[0], 1):
        ll = medianize_im(im[i, :, :])
        L.append(ll)
    return numpy.asarray(L)

def define_LOG(size, sensor, sigma=0.55248):
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


def define_LOG_generic(size, sigma):
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


# Fuer den MSAM?!?->Vectorized Correlate
def vcorrcoef(M, n):
    Mm = numpy.reshape(numpy.mean(M, axis=1), (M.shape[0], 1))
    nm = numpy.mean(n)
    r_zahl = numpy.sum((M - Mm) * (n - nm), axis=1)
    r_nen = numpy.sqrt(numpy.sum((M - Mm) ** 2, axis=1) * numpy.sum((n - nm) ** 2))
    r = r_zahl / r_nen
    r[~numpy.isfinite(r)] = 0
    return r

# Für 2d und nicht 3d Bildchen
class moving_window2d(object):
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


#Bakker Detector
def bakker_filter_log_MSAM(image, sensor, winsize=3):
    if winsize % 2 == 0:
        winsize = winsize - 1
    lsh = image.shape
    slidwin = moving_window(image, winsize)
    if sensor != None:
        Kantenfilter = define_LOG(winsize, sensor)
    else:
        Kantenfilter = define_LOG(winsize, sensor='hyperion')
    outimage = numpy.zeros([lsh[1], lsh[2]])
    for i in numpy.arange(0 + winsize // 2, lsh[1] - winsize // 2, 1):
        for j in numpy.arange(0 + winsize // 2, lsh[2] - winsize // 2, 1):
            slidwin.untermat(i, j)
            fensterchen = slidwin.fensterchen
            alles = fensterchen.reshape(slidwin.fsh[0], slidwin.fsh[1] * slidwin.fsh[2])
            samall = vcorrcoef(alles.T, slidwin.mitte)
            samall = samall.reshape(slidwin.fsh[1], slidwin.fsh[2])
            mit_zentrum = numpy.dot(Kantenfilter, samall)
            ohne_zentrum = copy.deepcopy(mit_zentrum)
            ohne_zentrum[winsize // 2, winsize // 2] = 0
            outimage[i, j] = numpy.sum(ohne_zentrum)
        # print i,'done'
    return outimage

##########
###########


def combined_ops_bakker_MSAM_cubo_finale(cubo, sensor, winsize=5):
    hyper_edge = bakker_filter_log_MSAM(cubo, sensor, winsize)
    return hyper_edge



def heda_archstyle(cubo, sensor, winsize=3):
    gradient= combined_ops_bakker_MSAM_cubo_finale(cubo, sensor, winsize)
    abs_grad = numpy.abs(gradient)
    absgradsmooth = ndimage.filters.median_filter(abs_grad, size=winsize + 2)
    gradsmooth = ndimage.filters.median_filter(gradient, size=winsize + 2)
    return gradient, abs_grad, absgradsmooth, gradsmooth
