import numpy

from osgeo import gdal, gdalnumeric
import auxfun as axf

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
"""

# wavelength={443, 482, 562, 655, 865, 1610, 2200}
# data[6,:]=band[5,:,:]
# band[4,:,:]=band[4,:,:]
# band[3,:,:]=band[3,:,:]
# band[2,:,:]=band[2,:,:]
# band[1,:,:]=band[1,:,:]
# band[0,:,:]=band[0,:,:]
# fwhm={20, 65, 75, 50, 40, 100, 200}
"""
def landsat7(filename):
    wave = [479, 560, 661, 834, 1650, 2209]
    fwhm = [70, 80, 60, 140, 200, 270]"""


def geoinfo(datahdr, targethdr):
    hdd = axf.read_hdr(datahdr)
    mf = hdd.map_info
    crs = hdd.coordinate_system_string
    mf = "map info ={" + mf + "}\n"
    crs = "coordinate system string = {" + crs + "}\n"
    L = []
    L.append(mf)
    L.append(crs)
    open(targethdr, 'a').writelines(L)
    return None


def writefile(origfile, band, suffix):
    axf.schreibeBSQsingle(numpy.nan_to_num(band), origfile + '_' + suffix)
    geoinfo(origfile + '.hdr', origfile + '_' + suffix + '.hdr')
    return None


def ls7_georat_single(file):
    band = gdalnumeric.LoadFile((file))
    band = band.astype(numpy.float32)
    ferric_iron1 = (band[5, :, :] / band[3, :, :]) + (band[1, :, :] / band[2, :, :])
    writefile(file, ferric_iron1, "ferric_iron1")
    del ferric_iron1
    ferrous_sillicates = (band[5, :, :] / band[4, :, :])
    writefile(file, ferrous_sillicates, "ferrous_sillicates")
    del ferrous_sillicates
    return None


def ls7_vegrat1_single(file):
    band = gdalnumeric.LoadFile((file))
    band = band.astype(numpy.float32)
    siwsi = (band[3, :, :] - band[4, :, :]) / (band[3, :, :] + band[4, :, :])
    writefile(file, siwsi, 'siwsi')
    del siwsi
    ndvi = (band[3, :, :] - band[2, :, :]) / (band[3, :, :] + band[2, :, :])
    writefile(file, ndvi, "ndvi")
    sbl = band[3, :, :] - 2.4 * band[2, :, :]
    writefile(file, sbl, 'sbl')
    del sbl
    rnr = band[2, :, :] / band[3, :, :]
    writefile(file, rnr, 'rnr')
    del rnr
    sw21 = band[5, :, :] / band[4, :, :]
    writefile(file, sw21, 'sw21')
    del sw21
    arvi = ((band[3, :, :] - band[2, :, :]) - (band[2, :, :] - band[0, :, :])) / (
                (band[3, :, :] + band[2, :, :]) - (band[2, :, :] - band[0, :, :]))
    writefile(file, arvi, "arvi")
    del arvi
    gvmi = ((band[3, :, :] + 0.1 * numpy.ones_like(ndvi)) - (band[4, :, :] + 0.02 * numpy.ones_like(ndvi))) / (
                (band[3, :, :] + 0.1 * numpy.ones_like(ndvi)) + (band[4, :, :] + 0.02 * numpy.ones_like(ndvi)))
    writefile(file, gvmi, "gvmi")
    del gvmi
    gli = (2 * band[1, :, :] - band[2, :, :] - band[0, :, :]) / (2 * band[1, :, :] + band[2, :, :] + band[0, :, :])
    writefile(file, gli, "gli")
    del gli
    ngi = band[1, :, :] / (band[3, :, :] + band[1, :, :] + band[2, :, :])
    writefile(file, ngi, 'ngi')
    del ngi
    nrot = band[2, :, :] / (band[3, :, :] + band[1, :, :] + band[2, :, :])
    writefile(file, nrot, 'nrot')
    del nrot
    pvr = (band[1, :, :] - band[2, :, :]) / (band[1, :, :] + band[2, :, :])
    writefile(file, pvr, 'pvr')
    del pvr
    ngrdi = (band[1, :, :] - band[2, :, :]) / (band[1, :, :] + band[2, :, :])
    writefile(file, ngrdi, 'ngrdi')
    del ngrdi
    sw2ndvi = (band[3, :, :] - band[5, :, :]) / (band[3, :, :] + band[5, :, :])
    writefile(file, sw2ndvi, 'sw2ndvi')
    del sw2ndvi
    sw1ndvi = (band[3, :, :] - band[4, :, :]) / (band[3, :, :] + band[4, :, :])
    writefile(file, sw1ndvi, 'sw1ndvi')
    del sw1ndvi
    sw21 = band[5, :, :] / band[4, :, :]
    writefile(file, sw21, 'sw21')
    del sw21
    wdvi = band[3, :, :] - 1.22 * band[2, :, :]
    writefile(file, wdvi, 'wdvi')
    del wdvi
    ndsi = (band[4, :, :] - band[5, :, :]) / (band[4, :, :] + band[5, :, :])
    writefile(file, ndsi, 'ndsi')
    del ndsi
    gvi = -0.2848 * band[0, :, :] - 0.2435 * band[1, :, :] - 0.5436 * band[2, :, :] + 0.7243 * band[3, :,
                                                                                               :] + 0.084 * band[4, :,
                                                                                                            :] - 0.18 * band[
                                                                                                                        5,
                                                                                                                        :,
                                                                                                                        :]
    writefile(file, gvi, 'gvi')
    del gvi
    wet = 0.1509 * band[0, :, :] + 0.1973 * band[1, :, :] + 0.3279 * band[2, :, :] + 0.3406 * band[3, :,
                                                                                              :] - 0.7112 * band[4, :,
                                                                                                            :] - 0.4572 * band[
                                                                                                                          5,
                                                                                                                          :,
                                                                                                                          :]
    writefile(file, wet, 'wet')
    del wet
    varigreen = (band[1, :, :] - band[2, :, :]) / (band[1, :, :] + band[2, :, :] - band[0, :, :])
    writefile(file, varigreen, 'varigreen')
    del varigreen
    return None
