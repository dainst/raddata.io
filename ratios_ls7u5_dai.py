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

#wavelength={443, 482, 562, 655, 865, 1610, 2200}
#data[6,:]=band[5,:,:]
#band[4,:,:]=band[4,:,:]
#band[3,:,:]=band[3,:,:]
#band[2,:,:]=band[2,:,:]
#band[1,:,:]=band[1,:,:]
#band[0,:,:]=band[0,:,:]
#fwhm={20, 65, 75, 50, 40, 100, 200}
"""
def landsat7(filename):
    wave = [479, 560, 661, 834, 1650, 2209]
    fwhm = [70, 80, 60, 140, 200, 270]"""
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

def writefile(origfile,band,suffix):
    axf.schreibeBSQsingle(numpy.nan_to_num(band),origfile+'_'+suffix)
    geoinfo(origfile+'.hdr',origfile+'_'+suffix+'.hdr')
    return None


def ls7_georat_single(file):
    band = gdalnumeric.LoadFile((file))
    band = band * 1.0
    ferric_iron1 = (band[5, :, :] / band[3, :, :]) + (band[1, :, :] / band[2, :, :])
    writefile(file,ferric_iron1,"ferric_iron1")
    ferric_iron2 = (band[2, :, :] / band[1, :, :])
    writefile(file,ferric_iron2,"ferric_iron2")
    ferric_oxides = (band[4, :, :] / band[3, :, :])
    writefile(file, ferric_oxides, "ferric_oxides")
    ferrous_iron = (band[5, :, :] / band[3, :, :]) + (band[2, :, :] / band[1, :, :])
    writefile(file, ferrous_iron, "ferrous_iron")
    ferrous_sillicates = (band[5, :, :] / band[4, :, :])
    writefile(file, ferrous_sillicates,"ferrous_sillicates")
    gossan = band[4, :, :] / band[2, :, :]
    writefile(file, gossan, "gossan")
    laterite_alt = band[4, :, :] / band[5, :, :]
    writefile(file, laterite_alt, "laterite_alteration")
    return None


def ls7_vegrat1_single(file):
    band=gdalnumeric.LoadFile((file))
    band=band*1.0
    afri1600=band[3,:,:]-0.66*(band[4,:,:])/(band[3,:,:]+0.66*band[4,:,:])
    writefile(file, afri1600, "afri1600")
    del afri1600
    avi=2.0*band[3,:,:]-band[2,:,:]
    writefile(file, avi, "avi")
    del avi
    arvi=((band[3,:,:]-band[2,:,:])-(band[2,:,:]-band[0,:,:]))/((band[3,:,:]+band[2,:,:])-(band[2,:,:]-band[0,:,:]))
    writefile(file, arvi, "arvi")
    del arvi
    arvi2=-0.18*numpy.ones_like(band[3,:,:])+1.17*((band[3,:,:]-band[2,:,:])/(band[3,:,:]+band[2,:,:]))
    writefile(file, arvi2, "arvi2")
    del arvi2
    bwdrvi=(0.1*(band[3,:,:]-band[0,:,:])/0.1*(band[3,:,:]+band[0,:,:]))
    writefile(file, bwdrvi, "bwdrvi")
    del bwdrvi
    clgreen=band[3,:,:]/band[1,:,:]-numpy.ones_like(band[1,:,:])
    writefile(file, clgreen, "clgreen")
    del clgreen
    cvi=(band[2,:,:]/band[1,:,:]**2)*band[3,:,:]
    writefile(file, cvi, "cvi")
    del cvi
    colorationindex=(band[2,:,:]-band[0,:,:])/band[2,:,:]
    writefile(file, colorationindex, "colorationindex")
    del colorationindex
    ndvi=(band[3,:,:]-band[2,:,:])/(band[3,:,:]+band[2,:,:])
    writefile(file, ndvi, "ndvi")
    ctvi=((ndvi+numpy.ones_like(ndvi))/numpy.abs(ndvi+numpy.ones_like(ndvi)))*numpy.sqrt(numpy.abs(ndvi+0.5*numpy.ones_like(ndvi)))
    writefile(file, ctvi, "ctvi")
    del ctvi
    cri550=1/band[0,:,:]-1/band[1,:,:]
    writefile(file, cri550, "cri550")
    del cri550
    gdvi = band[3, :, :] - band[1, :, :]
    writefile(file, gdvi, "gdvi")
    del gdvi
    dvimss = 2.4*band[3, :, :] - band[2, :, :]
    writefile(file, dvimss, "dvimss")
    del dvimss
    evi=2.5*((band[3,:,:]-band[2,:,:])/((band[3,:,:]+6*band[2,:,:]-7.5*band[0,:,:])+numpy.ones_like(ndvi)))
    writefile(file, evi, "evi")
    del evi
    evi2=2.4*((band[3,:,:]-band[2,:,:])/(band[3,:,:]+band[2,:,:]+numpy.ones_like(ndvi)))
    writefile(file, evi2, "evi2")
    del evi2
    evi22=2.5*((band[3,:,:]-band[2,:,:])/(band[3,:,:]+2.4*band[2,:,:]+numpy.ones_like(ndvi)))
    writefile(file, evi22, "evi22")
    del evi22
    gvmi=((band[3, :, :]+0.1*numpy.ones_like(ndvi))-(band[4, :, :]+0.02*numpy.ones_like(ndvi)))/((band[3, :, :]+0.1*numpy.ones_like(ndvi))+(band[4, :, :]+0.02*numpy.ones_like(ndvi)))
    writefile(file, gvmi, "gvmi")
    del gvmi
    gari=(band[3, :, :]-(band[1, :, :]-(band[0, :, :]-band[2, :, :])))/(band[3, :, :]-(band[1, :, :]+(band[0, :, :]-band[2, :, :])))
    writefile(file, gari, "gari")
    del gari
    gli=(2*band[1, :, :]-band[2, :, :]-band[0, :, :])/(2*band[1, :, :]+band[2, :, :]+band[0, :, :])
    writefile(file, gli, "gli")
    del gli
    gndvi =(band[3, :, :] - band[1, :, :]) / (band[3, :, :] + band[1, :, :])
    writefile(file, gndvi, "gndvi")
    del gndvi
    gbndvi=(band[3, :, :]-(band[1, :, :]+band[0, :, :]))/(band[3, :, :]+(band[1, :, :]+band[0, :, :]))
    writefile(file, gbndvi, "gbndvi")
    del gbndvi
    grndvi = (band[3, :, :] - (band[1, :, :] + band[2, :, :])) / (band[3, :, :] + (band[1, :, :] + band[2, :, :]))
    writefile(file, grndvi, "grndvi")
    del grndvi
    hue=numpy.arctan(((2*band[2,:,:]-band[1,:,:]-band[0,:,:])/30.5)*(band[1,:,:]-band[0,:,:]))
    writefile(file, hue, "hue")
    del hue
    ivi=(band[3, :, :]-0.03*numpy.ones_like(band[3,:,:]))/(1.22*band[2,:,:])
    writefile(file, ivi, "ivi")
    del ivi
    ipvi=((ndvi+numpy.ones_like(ndvi))*((band[3,:,:])))/(((band[3,:,:]+band[2,:,:]))/(2*numpy.ones_like(band[3,:,:])))
    writefile(file, ipvi, "ipvi")
    del ipvi
    intensity=1/30.5*(band[0,:,:]+band[1,:,:]+band[2,:,:])
    writefile(file, intensity, "intensity")
    del intensity
    logratio=numpy.log(band[3,:,:]/band[2,:,:])
    writefile(file, logratio, "logratio")
    del logratio
    msnirred=((band[3, :, :] / band[2, :, :])-numpy.ones_like(band[3,:,:]))/(numpy.sqrt((band[3, :, :] / band[2, :, :])+numpy.ones_like(band[3,:,:])))
    writefile(file, msnirred, "msnirred")
    del msnirred
    mcrig=(1/band[0,:,:]-1/band[1,:,:])*band[3,:,:]
    writefile(file, mcrig, "mcrig")
    msavi=(2*(band[3,:,:])+numpy.ones_like(band[3,:,:])-numpy.sqrt((2*band[3,:,:]+numpy.ones_like(band[2,:,:]))**2-8*(band[3,:,:])-band[2,:,:]))/2
    writefile(file,msavi,'msavi')
    del msavi
    nli=(band[3,:,:]**2-band[2,:,:])/(band[3,:,:]**2+band[2,:,:])
    writefile(file,nli,'nli')
    del nli
    ngi=band[1,:,:]/(band[3,:,:]+band[1,:,:]+band[2,:,:])
    writefile(file,ngi,'ngi')
    del ngi
    nnir=band[3,:,:]/(band[3,:,:]+band[1,:,:]+band[2,:,:])
    writefile(file,nnir,'nnir')
    del nnir
    nrot = band[2, :, :] / (band[3, :, :] + band[1, :, :] + band[2, :, :])
    writefile(file,nrot,'nrot')
    del nrot
    ppr=(band[1,:,:]-band[0,:,:])/(band[1,:,:]+band[0,:,:])
    writefile(file,ppr,'ppr')
    del ppr
    pvr=(band[1,:,:]-band[2,:,:])/(band[1,:,:]+band[2,:,:])
    writefile(file,pvr,'pvr')
    del pvr
    siwsi=(band[3,:,:]-band[4,:,:])/(band[3,:,:]+band[4,:,:])
    writefile(file,siwsi,'siwsi')
    del siwsi
    ngrdi = (band[1, :, :] - band[2, :, :]) / (band[1, :, :] + band[2, :, :])
    writefile(file,ngrdi,'ngrdi')
    del ngrdi
    bndvi= (band[3, :, :] - band[0, :, :]) / (band[3, :, :] + band[0, :, :])
    writefile(file,bndvi,'bndvi')
    del bndvi
    gndvi = (band[3, :, :] - band[1, :, :]) / (band[3, :, :] + band[1, :, :])
    writefile(file,gndvi,'gndvi')
    del gndvi
    sw2ndvi = (band[3, :, :] - band[5, :, :]) / (band[3, :, :] + band[5, :, :])
    writefile(file,sw2ndvi,'sw2ndvi')
    del sw2ndvi
    sw1ndvi = (band[3, :, :] - band[4, :, :]) / (band[3, :, :] + band[4, :, :])
    writefile(file,sw1ndvi,'sw1ndvi')
    del sw1ndvi
    ri = (band[2, :, :] - band[1, :, :]) / (band[2, :, :] + band[1, :, :])
    writefile(file,ri,'ri')
    del ri
    ndsi= (band[4, :, :] - band[5, :, :]) / (band[4, :, :] + band[5, :, :])
    writefile(file,ndsi,'ndsi')
    del ndsi
    panndvi=(band[3,:,:]-(band[2,:,:]+band[1,:,:]+band[0,:,:]))/(band[3,:,:]+(band[2,:,:]+band[1,:,:]+band[0,:,:]))
    writefile(file,panndvi,'panndvi')
    del panndvi
    rbndvi=(band[3,:,:]-(band[2,:,:]+band[0,:,:]))/(band[3,:,:]-(band[2,:,:]+band[0,:,:]))
    writefile(file,rbndvi,'rbndvi')
    del rbndvi
    IF=(2*band[2,:,:]-band[1,:,:]-band[0,:,:])/(band[1,:,:]-band[0,:,:])
    writefile(file,IF,'IF')
    del IF
    s1s2=band[4,:,:]/band[5,:,:]
    writefile(file,s1s2,'s1s2')
    del s1s2
    srgr=band[1,:,:]/band[2,:,:]
    writefile(file,srgr,'srgr')
    del srgr
    srgnir=band[3,:,:]/band[1,:,:]
    writefile(file,srgnir,'srgnir')
    del srgnir
    sw2nir=band[5,:,:]/band[3,:,:]
    writefile(file,sw2nir,'sw2nir')
    del sw2nir
    sw2r=band[5,:,:]/band[2,:,:]
    writefile(file,sw2r,'sw2r')
    del sw2r
    sw21=band[5,:,:]/band[4,:,:]
    writefile(file,sw21,'sw21')
    del sw21
    nirg=band[3,:,:]/band[1,:,:]
    writefile(file,nirg,'nirg')
    del nirg
    nirsw2=band[3,:,:]/band[5,:,:]
    writefile(file,nirsw2,'nirsw2')
    del nirsw2
    nirr=band[3,:,:]/band[2,:,:]
    writefile(file,nirr,'nirr')
    del nirr
    rb=band[2,:,:]/band[0,:,:]
    writefile(file,rb,'rb')
    del rb
    rg=band[2,:,:]/band[1,:,:]
    writefile(file,rg,'rg')
    del rg
    rnr=band[2,:,:]/band[3,:,:]
    writefile(file,rnr,'rnr')
    del rnr
    s1nir=band[4,:,:]/band[3,:,:]
    writefile(file,s1nir,'s1nir')
    del s1nir
    sarvi2=2.5*((band[3,:,:]-band[2,:,:])/(numpy.ones_like(band[0,:,:])+6*band[2,:,:]-7.5*band[0,:,:]))
    writefile(file,sarvi2,'sarvi2')
    del sarvi2
    sbl=band[3,:,:]-2.4*band[2,:,:]
    writefile(file,sbl,'sbl')
    del sbl
    sci=(band[4,:,:]-band[3,:,:])/(band[4,:,:]+band[3,:,:])
    writefile(file,sci,'sci')
    del sci
    savi2=band[3,:,:]/(band[2,:,:]+(1.22*0.03)*numpy.ones_like(band[0,:,:]))
    writefile(file,savi2,'savi2')
    del savi2
    slavi=band[3,:,:]/(band[2,:,:]+band[5,:,:])
    writefile(file,slavi,'slavi')
    del slavi
    sqrtirr=numpy.sqrt(band[3,:,:]/band[2,:,:])
    writefile(file,sqrtirr,'sqrtirr')
    del sqrtirr
    sbi=0.3037*band[0,:,:]+0.2793*band[1,:,:]+0.4743*band[2,:,:]+0.5585*band[3,:,:]+0.5082*band[4,:,:]+0.1863*band[5,:,:]
    writefile(file,sbi,'sbi')
    del sbi
    gvi=-0.2848*band[0,:,:]-0.2435*band[1,:,:]-0.5436*band[2,:,:]+0.7243*band[3,:,:]+0.084*band[4,:,:]-0.18*band[5,:,:]
    writefile(file,gvi,'gvi')
    del gvi
    wet = 0.1509 * band[0, :, :] +0.1973 * band[1, :, :] +0.3279 * band[2, :, :] + 0.3406 * band[3, :,:] - 0.7112 * band[4, :,:] - 0.4572 * band[5,:,:]
    writefile(file,wet,'wet')
    del wet
    tvi=numpy.sqrt(ndvi+0.5*numpy.ones_like(band[0,:,:]))
    writefile(file,tvi,'tvi')
    del tvi
    varigreen=(band[1,:,:]-band[2,:,:])/(band[1,:,:]+band[2,:,:]-band[0,:,:])
    writefile(file,varigreen,'varigreen')
    del varigreen
    wdvi=band[3,:,:]-1.22*band[2,:,:]
    writefile(file,wdvi,'wdvi')
    del wdvi
    wdrvi=(0.1*band[3,:,:]-band[3,:,:])/(0.1*band[3,:,:]+band[3,:,:])
    writefile(file,wdrvi,"wdrvi")
    del wdrvi
    return None

