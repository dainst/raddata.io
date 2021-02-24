import gdalnumeric,numpy

from osgeo import gdal
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
#LS8
#wavelength={443, 482, 562, 655, 865, 1610, 2200}
#fwhm={20, 65, 75, 50, 40, 100, 200}
"""
def sent2_20m(filename):
    wave = [490, 560, 665, 705, 740, 783, 865, 1610, 2190]
    fwhm = [65, 35, 30, 15, 15, 20, 20, 90, 180]
    dta = recalc_spectra_file(filename, wave, fwhm)
    return dta, wave, fwhm


def sent2_10m(filename):
    wave = [490, 560, 665, 842]
    fwhm = [65, 35, 30, 115]
    dta = recalc_spectra_file(filename, wave, fwhm)
    return dta, wave, fwhm

"""


def geoinfo(datahdr,targethdr):
    hdd=axf.read_hdr(datahdr)
    mf=hdd.map_info
    crs=hdd.coordinate_system_string
    mf="map info ={"+mf+"}\n"
    crs="coordinate system string = {"+crs+"}\n"
    L=[]
    L.append(mf)
    L.append(crs)
    open(targethdr,'a').writelines(L)
    return None

def writefile(origfile,data,suffix):
    axf.schreibeBSQsingle(numpy.nan_to_num(data),origfile+'_'+suffix)
    geoinfo(origfile+'.hdr',origfile+'_'+suffix+'.hdr')
    return None

"""
wavelength={490, 560, 665, 705, 740, 783, 865, 1610, 2190}
2, 3, 4, 5 6 7 8A, 11, 12
0, 1, 2, 3,4,5,6  ,7, 8"""

def S2_20mgeorat_single(file):
    data = gdalnumeric.LoadFile((file))
    data = data * 1.0
    ferric_iron1 = (data[8, :, :] / data[5, :, :]) + (data[1, :, :] / data[2, :, :])
    writefile(file,ferric_iron1,"ferric_iron1")
    ferric_oxides = (data[7, :, :] / data[6, :, :])
    writefile(file, ferric_oxides, "ferric_oxides")
    ferrous_iron = (data[8, :, :] / data[5, :, :]) + (data[1, :, :] / data[2, :, :])
    writefile(file, ferrous_iron, "ferrous_iron")
    ferrous_sillicates = (data[8, :, :] / data[7, :, :])
    writefile(file, ferrous_sillicates,"ferrous_sillicates")
    return None



def S2_10m_vegrat1_single(file):
    data=gdalnumeric.LoadFile((file))
    data=data*1.0
    avi=2.0*data[3,:,:]-data[2,:,:]
    writefile(file, avi, "avi")
    del avi
    ndvi=(data[3,:,:]-data[2,:,:])/(data[3,:,:]+data[2,:,:])
    writefile(file, ndvi, "ndvi")
    colorationindex=(data[2,:,:]-data[0,:,:])/data[2,:,:]
    writefile(file, colorationindex, "colorationindex")
    del colorationindex
    gli=(2*data[1, :, :]-data[2, :, :]-data[0, :, :])/(2*data[1, :, :]+data[2, :, :]+data[0, :, :])
    writefile(file, gli, "gli")
    del gli
    intensity=1/30.5*(data[0,:,:]+data[1,:,:]+data[2,:,:])
    writefile(file, intensity, "intensity")
    del intensity
    msavi=(2*(data[3,:,:])+numpy.ones_like(data[3,:,:])-numpy.sqrt((2*data[3,:,:]+numpy.ones_like(data[2,:,:]))**2-8*(data[3,:,:])-data[2,:,:]))/2
    writefile(file,msavi,'msavi')
    del msavi
    ngrdi = (data[1, :, :] - data[2, :, :]) / (data[1, :, :] + data[2, :, :])
    writefile(file,ngrdi,'ngrdi')
    del ngrdi
    ng=data[1,:,:]/(data[3,:,:]+data[1,:,:]+data[2,:,:])
    writefile(file,ng,'ng')
    del ng
    IF=(2*data[2,:,:]-data[1,:,:]-data[0,:,:])/(data[1,:,:]-data[0,:,:])
    writefile(file,IF,'IF')
    del IF
    nrot = data[2, :, :] / (data[3, :, :] + data[1, :, :] + data[2, :, :])
    writefile(file,nrot,'nrot')
    del nrot
    pvr=(data[1,:,:]-data[2,:,:])/(data[1,:,:]+data[2,:,:])
    writefile(file,pvr,'pvr')
    del pvr
    srgr=data[1,:,:]/data[2,:,:]
    writefile(file,srgr,'srgr')
    del srgr
    ri = (data[2, :, :] - data[1, :, :]) / (data[2, :, :] + data[1, :, :])
    writefile(file,ri,'ri')
    del ri
    bgi=data[0,:,:]/data[1,:,:]
    writefile(file,bgi,'bgi')
    del bgi
    rnr=data[2,:,:]/data[3,:,:]
    writefile(file,rnr,'rnr')
    ngrdi = (data[1, :, :] - data[2, :, :]) / (data[1, :, :] + data[2, :, :])
    writefile(file,ngrdi,'ngrdi')
    del rnr
    varigreen=(data[1,:,:]-data[2,:,:])/(data[0,:,:]+data[2,:,:]-data[0,:,:])
    writefile(file,varigreen,'varigreen')
    del varigreen
    return None

def S2_20m_vegrat1_single(file):
    data=gdalnumeric.LoadFile((file))
    data=data*1.0
    afri1600=data[6,:,:]-0.66*(data[7,:,:])/(data[6,:,:]+0.66*data[7,:,:])
    writefile(file, afri1600, "afri1600")
    del afri1600
    avi=2.0*data[6,:,:]-data[2,:,:]
    writefile(file, avi, "avi")
    ndvi=(data[6,:,:]-data[2,:,:])/(data[6,:,:]+data[2,:,:])
    writefile(file, ndvi, "ndvi")
    intensity=1/30.5*(data[0,:,:]+data[1,:,:]+data[2,:,:])
    writefile(file, intensity, "intensity")
    del intensity
    msavi=(2*(data[6,:,:])+numpy.ones_like(data[6,:,:])-numpy.sqrt((2*data[6,:,:]+numpy.ones_like(data[2,:,:]))**2-8*(data[6,:,:])-data[2,:,:]))/2
    writefile(file,msavi,'msavi')
    del msavi
    ng=data[1,:,:]/(data[6,:,:]+data[1,:,:]+data[2,:,:])
    writefile(file,ng,'ng')
    del ng
    nrot = data[2, :, :] / (data[6, :, :] + data[1, :, :] + data[2, :, :])
    writefile(file,nrot,'nrot')
    del nrot
    pvr=(data[1,:,:]-data[2,:,:])/(data[1,:,:]+data[2,:,:])
    writefile(file,pvr,'pvr')
    del pvr
    ngrdi = (data[1, :, :] - data[2, :, :]) / (data[1, :, :] + data[2, :, :])
    writefile(file,ngrdi,'ngrdi')
    del ngrdi
    bgi=data[0,:,:]/data[1,:,:]
    writefile(file,bgi,'bgi')
    del bgi
    srgr=data[1,:,:]/data[2,:,:]
    writefile(file,srgr,'srgr')
    del srgr
    sw21=data[8,:,:]/data[6,:,:]
    writefile(file,sw21,'sw21')
    del sw21
    nirg=data[6,:,:]/data[1,:,:]
    writefile(file,nirg,'nirg')
    del nirg
    rnr=data[2,:,:]/data[6,:,:]
    writefile(file,rnr,'rnr')
    del rnr
    s1nir=data[7,:,:]/data[6,:,:]
    writefile(file,s1nir,'s1nir')
    del s1nir
    sw2nir=data[8,:,:]/data[6,:,:]
    writefile(file,sw2nir,'sw2nir')
    del sw2nir
    sci=(data[7,:,:]-data[6,:,:])/(data[7,:,:]+data[6,:,:])
    writefile(file,sci,'sci')
    del sci         
    sbi=0.3037*data[0,:,:]+0.2793*data[1,:,:]+0.4743*data[2,:,:]+0.5585*data[6,:,:]+0.5082*data[7,:,:]+0.1863*data[8,:,:]
    writefile(file,sbi,'sbi')
    del sbi
    return None
