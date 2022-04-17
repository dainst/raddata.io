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
# fwhm={20, 65, 75, 50, 40, 100, 200}

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


def writefile(origfile, data, suffix):
    axf.schreibeBSQsingle(numpy.nan_to_num(data), origfile + '_' + suffix)
    geoinfo(origfile + '.hdr', origfile + '_' + suffix + '.hdr')
    return None


def ls8_georat_single(file):
    data = gdalnumeric.LoadFile((file))
    data = data.astype(numpy.float32)
    ferric_iron1 = (data[6, :, :] / data[4, :, :]) + (data[2, :, :] / data[3, :, :])
    writefile(file, ferric_iron1, "ferric_iron1")
    del ferric_iron1
    ferrous_sillicates = (data[6, :, :] / data[5, :, :])
    writefile(file, ferrous_sillicates, "ferrous_sillicates")
    del ferrous_sillicates
    return None


def ls8_vegrat1_single(file):
    data = gdalnumeric.LoadFile((file))
    data = data.astype(numpy.float32)
    arvi = ((data[4, :, :] - data[3, :, :]) - (data[3, :, :] - data[1, :, :])) / (
                (data[4, :, :] + data[3, :, :]) - (data[3, :, :] - data[1, :, :]))
    writefile(file, arvi, "arvi")
    del arvi
    arvi2 = -0.18 * numpy.ones_like(data[4, :, :]) + 1.17 * (
                (data[4, :, :] - data[3, :, :]) / (data[4, :, :] + data[3, :, :]))
    writefile(file, arvi2, "arvi2")
    del arvi2
    gli = (2 * data[2, :, :] - data[3, :, :] - data[1, :, :]) / (2 * data[2, :, :] + data[3, :, :] + data[1, :, :])
    writefile(file, gli, "gli")
    del gli
    ng = data[2, :, :] / (data[4, :, :] + data[2, :, :] + data[3, :, :])
    writefile(file, ng, 'ng')
    del ng
    ndvi = (data[4, :, :] - data[3, :, :]) / (data[4, :, :] + data[3, :, :])
    writefile(file, ndvi, "ndvi")
    pvr = (data[2, :, :] - data[3, :, :]) / (data[2, :, :] + data[3, :, :])
    writefile(file, pvr, 'pvr')
    del pvr
    siwsi = (data[4, :, :] - data[5, :, :]) / (data[4, :, :] + data[5, :, :])
    writefile(file, siwsi, 'siwsi')
    del siwsi
    ngrdi = (data[2, :, :] - data[3, :, :]) / (data[2, :, :] + data[3, :, :])
    writefile(file, ngrdi, 'ngrdi')
    del ngrdi
    sw2ndvi = (data[4, :, :] - data[6, :, :]) / (data[4, :, :] + data[6, :, :])
    writefile(file, sw2ndvi, 'sw2ndvi')
    del sw2ndvi
    sw21 = data[6, :, :] / data[5, :, :]
    writefile(file, sw21, 'sw21')
    del sw21
    rnr = data[3, :, :] / data[4, :, :]
    writefile(file, rnr, 'rnr')
    del rnr
    sbl = data[4, :, :] - 2.4 * data[3, :, :]
    writefile(file, sbl, 'sbl')
    del sbl
    sw2ndvi = (data[4, :, :] - data[6, :, :]) / (data[4, :, :] + data[6, :, :])
    writefile(file, sw2ndvi, 'sw2ndvi')
    del sw2ndvi
    wet = 0.1509 * data[1, :, :] + 0.1973 * data[2, :, :] + 0.3279 * data[3, :, :] + 0.3406 * data[4, :,
                                                                                              :] - 0.7112 * data[5, :,
                                                                                                            :] - 0.4572 * data[
                                                                                                                          6,
                                                                                                                          :,
                                                                                                                          :]
    writefile(file, wet, 'wet')
    del wet
    gvmi = ((data[4, :, :] + 0.1 * numpy.ones_like(ndvi)) - (data[5, :, :] + 0.02 * numpy.ones_like(ndvi))) / (
                (data[4, :, :] + 0.1 * numpy.ones_like(ndvi)) + (data[5, :, :] + 0.02 * numpy.ones_like(ndvi)))
    writefile(file, gvmi, "gvmi")
    del gvmi
    varigreen = (data[2, :, :] - data[3, :, :]) / (data[2, :, :] + data[3, :, :] - data[1, :, :])
    writefile(file, varigreen, 'varigreen')
    del varigreen
    wdvi = data[4, :, :] - 1.22 * data[3, :, :]
    writefile(file, wdvi, 'wdvi')
    del wdvi
    return None
