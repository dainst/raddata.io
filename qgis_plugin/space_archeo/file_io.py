import numpy
from osgeo import gdal

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


class FileIo:

    @staticmethod
    def schreibeBSQ(numpyarray: numpy.ndarray, format: str, out: str, dtype) -> None:
        """
        Function for reading and writing raster data.
        :param numpyarray:
        :param format: Output Format: "ENVI", "GTIFF", It is recomennded to avoid "GTIFF" and use "ENVI" in order to offer greater flexibility with metadata
        :param out: Outputfile
        :param dtype: gdal.GDT_Float32, gdal.GDT_UInt16, gdal.GDT_Float64, etc.
        :return: None
        """
        driver = gdal.GetDriverByName(format)
        h = numpy.shape(numpyarray)
        if len(h > 3):
            raise IndexError('Too many Dimensions funtions are for 2D and 3D arrays _only_!')
        if len(h) == 2:
            dst_ds = driver.Create(out, h[1], h[0], 1, dtype)
            dst_ds.GetRasterBand(1).WriteArray(numpyarray)
            dst_ds = None
        if len(h) == 3:
            dst_ds = driver.Create(out, h[2], h[1], h[0], dtype)
            ka = 1
            while ka <= h[0]:
                hhh = ka - 1
                dst_ds.GetRasterBand(ka).WriteArray(numpyarray[hhh, :, :])
                ka += 1
            dst_ds = None

    @staticmethod
    def update_CRS(xres, yres, nodata):
        pass


class MetaData:
    def __init__(self, filename:str):
        """
        Loads Metadata from ENVI HDR Files to NUMPY nd-arrays
        :param filename: Filename
        """
        if type(filename) == str:
            self.filename = filename
            self.r_h = MetaData.return_header(self.filename)
            self.obj = MetaData.splitter_hd(self.r_h)
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
                    self.wavelength = self.str2npar(self.wavelength)
            except(AttributeError):
                print("No wavelength Keyword in HDR")
            try:
                self.Wavelength
                if type(self.Wavelength) == str:
                    self.Wavelength = self.str2npar(self.Wavelength)
            except(AttributeError):
                print("No Wavelength Keyword in HDR")
            try:
                self.fwhm
                if type(self.fwhm) == str:
                    self.fwhm = self.str2npar(self.fwhm)
            except(AttributeError):
                print("No fwhm Keyword in HDR")
            try:
                self.bbl
                if type(self.bbl) == str:
                    self.bbl = self.str2npar(self.bbl)
            except(AttributeError):
                print("No bbl Keyword in HDR")

    @staticmethod
    def return_header(filename):
        h = open(filename, 'r').readlines()
        j = []
        g = ''
        for i in h:
            i = i.replace('\n', '')
            i = i.replace('\r\n', '')
            g = g + i
            if ('{' and '}' in g) or ('=' and '{' not in g):
                j.append(g)
                g = ''
        return j

    @staticmethod
    def return_header_list(filename):
        h = open(filename, 'r').readlines()
        j = []
        g = ''
        for i in h:
            g = g + i
            if ('{' and '}' in g) or ('=' and '{' not in g):
                j.append(g)
                g = ''
        return j

    @staticmethod
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

    @staticmethod
    def str2npar(string):
        ah = string.split(',')
        data = numpy.asarray(ah)
        data = data.astype(float)
        return data

    @staticmethod
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