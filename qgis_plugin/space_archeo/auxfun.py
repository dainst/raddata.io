from osgeo import gdal
import numpy


"""
    Date                 : Juli 2019
    Copyright            : (C) 2019 by Christian Mielke
    Email                : christian.mielke@gfz-potsdam.de
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

def schreibeBSQ(numpyarray,out):
    format = "ENVI"
    driver = gdal.GetDriverByName( format )
    h=numpy.shape(numpyarray)
    dst_ds = driver.Create(out, h[2], h[1], h[0], gdal.GDT_Float32 )
    ka=1
    while ka<=h[0]:
        hhh=ka-1
        dst_ds.GetRasterBand(ka).WriteArray(numpyarray[hhh,:,:])
        ka+=1
    dst_ds=None

def schreibeBSQ_int(numpyarray,out):
    format = "ENVI"
    driver = gdal.GetDriverByName( format )
    h=numpy.shape(numpyarray)
    dst_ds = driver.Create(out, h[2], h[1], h[0], gdal.GDT_UInt16)
    ka=1
    while ka<=h[0]:
        hhh=ka-1
        dst_ds.GetRasterBand(ka).WriteArray(numpyarray[hhh,:,:])
        ka+=1
    dst_ds=None

def schreibeBSQsingle(numpyarray,out):
    format = "ENVI"
    driver = gdal.GetDriverByName(format)
    h=numpy.shape(numpyarray)
    dst_ds = driver.Create(out,h[1],h[0],1, gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(numpyarray)
    dst_ds=None


##### I need this to parse header inforation:
###########
def return_header(filename):
    h = open(filename, 'r').readlines()
    j = []
    g = ''
    for i in h:
        i = i.replace('\n', '')
        i = i.replace('\r\n', '')
        g = g+i
        if ('{' and '}' in g) or ('=' and '{' not in g):
            j.append(g)
            g = ''
    return j


def return_header_list(filename):
    h = open(filename, 'r').readlines()
    j = []
    g = ''
    for i in h:
        g = g+i
        if ('{' and '}' in g) or ('=' and '{' not in g):
            j.append(g)
            g = ''
    return j


def splitter_hd(liste):
    hd_liste = {}
    prp = []
    val = []
    for i in liste:
        if '=' in i:
            i = i.replace('{', '')
            i = i.replace('}', '')
            i = i.replace('\r', '')
            i = i.replace('\r\n', '')
            i = i.replace('\n', '')
            i = i.strip()
            i = i.split('=')
            i[0] = i[0].strip()
            i[1] = i[1].strip()
            prp.append(i[0])
            val.append(i[1])
            hd_liste.update({i[0]: i[1]})
    w = hd_liste.keys()
    return hd_liste, w, prp, val

def splitter_hd2(liste):
    hd_liste = {}
    prp = []
    val = []
    for i in liste:
        if '=' in i:
            i = i.replace('{', '')
            i = i.replace('}', '')
            i = i.replace('\r', '')
            i = i.replace('\r\n', '')
            i = i.replace('\n', '')
            i = i.strip()
            i = i.split('=')
            i[0] = i[0].strip()
            i[1] = i[1].strip()
            prp.append(i[0])
            try:
                val.append(float(i[1]))
            except(ValueError):
                val.append(i[1])
            hd_liste.update({i[0]: i[1]})
    w = hd_liste.keys()
    return hd_liste, w, prp, val


class read_hdr(object):
    def __init__(self, filename=None):
        if type(filename) == str:
            self.filename = filename
            self.r_h = return_header(self.filename)
            self.obj = splitter_hd(self.r_h)
            self.li = self.obj[1]
            self.di = self.obj[0]
            self.prp = self.obj[2]
            self.val = self.obj[3]
            for i in self.li:
                setattr(self, i.replace(' ', '_'), self.di.get(i))

    def make_liste_w(self):
        pass


def str2npar(string):
    ah = string.split(',')
    data = numpy.asarray(ah)
    data = data.astype(float)
    return data


class read_hdr_flt(object):  #Read in Header and give its values back as floats
    def __init__(self, filename=None):
        if type(filename) == str:
            self.filename = filename
            self.r_h = return_header(self.filename)
            self.obj = splitter_hd(self.r_h)
            self.li = self.obj[1]
            self.di = self.obj[0]
            self.prp = self.obj[2]
            self.val = self.obj[3]
            for i in self.li:
                try:
                    setattr(self, i.replace(' ', '_'), float(self.di.get(i)))
                except(ValueError):
                    setattr(self, i.replace(' ', '_'), self.di.get(i))
            try:
                self.wavelength
                if type(self.wavelength) == str:
                    self.wavelength = str2npar(self.wavelength)
            except(AttributeError):
                pass
            try:
                self.Wavelength
                if type(self.Wavelength) == str:
                    self.Wavelength = str2npar(self.Wavelength)
            except(AttributeError):
                pass
            try:
                self.fwhm
                if type(self.fwhm) == str:
                    self.fwhm = str2npar(self.fwhm)
            except(AttributeError):
                pass
            try:
                self.bbl
                if type(self.bbl) == str:
                    self.bbl = str2npar(self.bbl)
            except(AttributeError):
                pass

    def make_liste_w(self):
        pass
########################
############################
#Header Parsing FInished!
###############################

#########I need these functions to prepare the data
################

