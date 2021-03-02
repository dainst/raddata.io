import gdalnumeric,numpy
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
#fwhm={20, 65, 75, 50, 40, 100, 200}

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


def ls8_georat_single(file):
    data = gdalnumeric.LoadFile((file))
    data = data * 1.0
    ferric_iron1 = (data[6, :, :] / data[4, :, :]) + (data[2, :, :] / data[3, :, :])
    writefile(file,ferric_iron1,"ferric_iron1")
    ferric_iron2 = (data[3, :, :] / data[2, :, :])
    writefile(file,ferric_iron2,"ferric_iron2")
    ferric_oxides = (data[5, :, :] / data[4, :, :])
    writefile(file, ferric_oxides, "ferric_oxides")
    ferrous_iron = (data[6, :, :] / data[4, :, :]) + (data[3, :, :] / data[2, :, :])
    writefile(file, ferrous_iron, "ferrous_iron")
    ferrous_sillicates = (data[6, :, :] / data[5, :, :])
    writefile(file, ferrous_sillicates,"ferrous_sillicates")
    gossan = data[5, :, :] / data[3, :, :]
    writefile(file, gossan, "gossan")
    laterite_alt = data[5, :, :] / data[6, :, :]
    writefile(file, laterite_alt, "laterite_alteration")
    return None


def ls8_vegrat1_single(file):
    data=gdalnumeric.LoadFile((file))
    data=data*1.0
    atsavi=1.22*((data[4,:,:]-1.22*data[3,:,:]-0.03*numpy.ones_like(data[0,:,:]))/(1.22*data[4,:,:]+data[3,:,:]-1.22*0.03*numpy.ones_like(data[0,:,:]+(0.08+0.08*1.22**2)*numpy.ones_like(data[0,:,:]))))
    writefile(file, atsavi, "atsavi")
    del atsavi
    afri1600=data[4,:,:]-0.66*(data[5,:,:])/(data[4,:,:]+0.66*data[5,:,:])
    writefile(file, afri1600, "afri1600")
    del afri1600
    avi=2.0*data[4,:,:]-data[3,:,:]
    writefile(file, avi, "avi")
    del avi
    arvi=((data[4,:,:]-data[3,:,:])-(data[3,:,:]-data[1,:,:]))/((data[4,:,:]+data[3,:,:])-(data[3,:,:]-data[1,:,:]))
    writefile(file, arvi, "arvi")
    del arvi
    arvi2=-0.18*numpy.ones_like(data[4,:,:])+1.17*((data[4,:,:]-data[3,:,:])/(data[4,:,:]+data[3,:,:]))
    writefile(file, arvi2, "arvi2")
    del arvi2
    bwdrvi=(0.1*(data[4,:,:]-data[1,:,:])/0.1*(data[4,:,:]+data[1,:,:]))
    writefile(file, bwdrvi, "bwdrvi")
    del bwdrvi
    clgreen=data[4,:,:]/data[2,:,:]-numpy.ones_like(data[2,:,:])
    writefile(file, clgreen, "clgreen")
    del clgreen
    cvi=(data[3,:,:]/data[2,:,:]**2)*data[4,:,:]
    writefile(file, cvi, "cvi")
    del cvi
    colorationindex=(data[3,:,:]-data[1,:,:])/data[3,:,:]
    writefile(file, colorationindex, "colorationindex")
    del colorationindex
    ndvi=(data[4,:,:]-data[3,:,:])/(data[4,:,:]+data[3,:,:])
    writefile(file, ndvi, "ndvi")
    ctvi=((ndvi+numpy.ones_like(ndvi))/numpy.abs(ndvi+numpy.ones_like(ndvi)))*numpy.sqrt(numpy.abs(ndvi+0.5*numpy.ones_like(ndvi)))
    writefile(file, ctvi, "ctvi")
    del ctvi
    cri550=1/data[1,:,:]-1/data[2,:,:]
    writefile(file, cri550, "cri550")
    del cri550
    gdvi = data[4, :, :] - data[2, :, :]
    writefile(file, gdvi, "gdvi")
    del gdvi
    dvimss = 2.4*data[4, :, :] - data[3, :, :]
    writefile(file, dvimss, "dvimss")
    del dvimss
    evi=2.5*((data[4,:,:]-data[3,:,:])/((data[4,:,:]+6*data[3,:,:]-7.5*data[1,:,:])+numpy.ones_like(ndvi)))
    writefile(file, evi, "evi")
    del evi
    evi2=2.4*((data[4,:,:]-data[3,:,:])/(data[4,:,:]+data[3,:,:]+numpy.ones_like(ndvi)))
    writefile(file, evi2, "evi2")
    del evi2
    evi22=2.5*((data[4,:,:]-data[3,:,:])/(data[4,:,:]+2.4*data[3,:,:]+numpy.ones_like(ndvi)))
    writefile(file, evi22, "evi22")
    del evi22
    gvmi=((data[4, :, :]+0.1*numpy.ones_like(ndvi))-(data[5, :, :]+0.02*numpy.ones_like(ndvi)))/((data[4, :, :]+0.1*numpy.ones_like(ndvi))+(data[5, :, :]+0.02*numpy.ones_like(ndvi)))
    writefile(file, gvmi, "gvmi")
    del gvmi
    gari=(data[4, :, :]-(data[2, :, :]-(data[1, :, :]-data[3, :, :])))/(data[4, :, :]-(data[2, :, :]+(data[1, :, :]-data[3, :, :])))
    writefile(file, gari, "gari")
    del gari
    gli=(2*data[2, :, :]-data[3, :, :]-data[1, :, :])/(2*data[2, :, :]+data[3, :, :]+data[1, :, :])
    writefile(file, gli, "gli")
    del gli
    gndvi =(data[4, :, :] - data[2, :, :]) / (data[4, :, :] + data[2, :, :])
    writefile(file, gndvi, "gndvi")
    del gndvi
    gbndvi=(data[4, :, :]-(data[2, :, :]+data[1, :, :]))/(data[4, :, :]+(data[2, :, :]+data[1, :, :]))
    writefile(file, gbndvi, "gbndvi")
    del gbndvi
    grndvi = (data[4, :, :] - (data[2, :, :] + data[3, :, :])) / (data[4, :, :] + (data[2, :, :] + data[3, :, :]))
    writefile(file, grndvi, "grndvi")
    del grndvi
    hue=numpy.arctan(((2*data[3,:,:]-data[2,:,:]-data[1,:,:])/30.5)*(data[2,:,:]-data[1,:,:]))
    writefile(file, hue, "hue")
    del hue
    ivi=(data[4, :, :]-0.03*numpy.ones_like(data[4,:,:]))/(1.22*data[3,:,:])
    writefile(file, ivi, "ivi")
    del ivi
    ipvi=((ndvi+numpy.ones_like(ndvi))*((data[4,:,:])))/(((data[4,:,:]+data[3,:,:]))/(2*numpy.ones_like(data[4,:,:])))
    writefile(file, ipvi, "ipvi")
    del ipvi
    intensity=1/30.5*(data[1,:,:]+data[2,:,:]+data[3,:,:])
    writefile(file, intensity, "intensity")
    del intensity
    logratio=numpy.log(data[4,:,:]/data[3,:,:])
    writefile(file, logratio, "logratio")
    del logratio
    msnirred=((data[4, :, :] / data[3, :, :])-numpy.ones_like(data[4,:,:]))/(numpy.sqrt((data[4, :, :] / data[3, :, :])+numpy.ones_like(data[4,:,:])))
    writefile(file, msnirred, "msnirred")
    del msnirred
    mcrig=(1/data[1,:,:]-1/data[2,:,:])*data[4,:,:]
    writefile(file, mcrig, "mcrig")
    msavi=(2*(data[4,:,:])+numpy.ones_like(data[4,:,:])-numpy.sqrt((2*data[4,:,:]+numpy.ones_like(data[3,:,:]))**2-8*(data[4,:,:])-data[3,:,:]))/2
    writefile(file,msavi,'msavi')
    del msavi
    nli=(data[4,:,:]**2-data[3,:,:])/(data[4,:,:]**2+data[3,:,:])
    writefile(file,nli,'nli')
    del nli
    ng=data[2,:,:]/(data[4,:,:]+data[2,:,:]+data[3,:,:])
    writefile(file,ng,'ng')
    del ng
    nnir=data[4,:,:]/(data[4,:,:]+data[2,:,:]+data[3,:,:])
    writefile(file,nnir,'nnir')
    del nnir
    nrot = data[3, :, :] / (data[4, :, :] + data[2, :, :] + data[3, :, :])
    writefile(file,nrot,'nrot')
    del nrot
    ppr=(data[2,:,:]-data[1,:,:])/(data[2,:,:]+data[1,:,:])
    writefile(file,ppr,'ppr')
    del ppr
    pvr=(data[2,:,:]-data[3,:,:])/(data[2,:,:]+data[3,:,:])
    writefile(file,pvr,'pvr')
    del pvr
    siwsi=(data[4,:,:]-data[5,:,:])/(data[4,:,:]+data[5,:,:])
    writefile(file,siwsi,'siwsi')
    del siwsi
    ngrdi = (data[2, :, :] - data[3, :, :]) / (data[2, :, :] + data[3, :, :])
    writefile(file,ngrdi,'ngrdi')
    del ngrdi
    bndvi= (data[4, :, :] - data[1, :, :]) / (data[4, :, :] + data[1, :, :])
    writefile(file,bndvi,'bndvi')
    del bndvi
    gndvi = (data[4, :, :] - data[2, :, :]) / (data[4, :, :] + data[2, :, :])
    writefile(file,gndvi,'gndvi')
    del gndvi
    sw2ndvi = (data[4, :, :] - data[6, :, :]) / (data[4, :, :] + data[6, :, :])
    writefile(file,sw2ndvi,'sw2ndvi')
    del sw2ndvi
    sw1ndvi = (data[4, :, :] - data[5, :, :]) / (data[4, :, :] + data[5, :, :])
    writefile(file,sw1ndvi,'sw1ndvi')
    del sw1ndvi
    ri = (data[3, :, :] - data[2, :, :]) / (data[3, :, :] + data[2, :, :])
    writefile(file,ri,'ri')
    del ri
    ndsi= (data[5, :, :] - data[6, :, :]) / (data[5, :, :] + data[6, :, :])
    writefile(file,ndsi,'ndsi')
    del ndsi
    panndvi=(data[4,:,:]-(data[3,:,:]+data[2,:,:]+data[1,:,:]))/(data[4,:,:]+(data[3,:,:]+data[2,:,:]+data[1,:,:]))
    writefile(file,panndvi,'panndvi')
    del panndvi
    rbndvi=(data[4,:,:]-(data[3,:,:]+data[1,:,:]))/(data[4,:,:]-(data[3,:,:]+data[1,:,:]))
    writefile(file,rbndvi,'rbndvi')
    del rbndvi
    IF=(2*data[3,:,:]-data[2,:,:]-data[1,:,:])/(data[2,:,:]-data[1,:,:])
    writefile(file,IF,'IF')
    del IF
    s1s2=data[5,:,:]/data[6,:,:]
    writefile(file,s1s2,'s1s2')
    del s1s2
    bgi=data[0,:,:]/data[2,:,:]
    writefile(file,bgi,'bgi')
    del bgi
    srgr=data[2,:,:]/data[3,:,:]
    writefile(file,srgr,'srgr')
    del srgr
    srgnir=data[4,:,:]/data[2,:,:]
    writefile(file,srgnir,'srgnir')
    del srgnir
    sw2nir=data[6,:,:]/data[4,:,:]
    writefile(file,sw2nir,'sw2nir')
    del sw2nir
    sw2r=data[6,:,:]/data[3,:,:]
    writefile(file,sw2r,'sw2r')
    del sw2r
    sw21=data[6,:,:]/data[5,:,:]
    writefile(file,sw21,'sw21')
    del sw21
    nirg=data[4,:,:]/data[2,:,:]
    writefile(file,nirg,'nirg')
    del nirg
    nirsw2=data[4,:,:]/data[6,:,:]
    writefile(file,nirsw2,'nirsw2')
    del nirsw2
    nirr=data[4,:,:]/data[3,:,:]
    writefile(file,nirr,'nirr')
    del nirr
    rb=data[3,:,:]/data[1,:,:]
    writefile(file,rb,'rb')
    del rb
    rg=data[3,:,:]/data[2,:,:]
    writefile(file,rg,'rg')
    del rg
    rnr=data[3,:,:]/data[4,:,:]
    writefile(file,rnr,'rnr')
    del rnr
    s1nir=data[5,:,:]/data[4,:,:]
    writefile(file,s1nir,'s1nir')
    del s1nir
    sarvi2=2.5*((data[4,:,:]-data[3,:,:])/(numpy.ones_like(data[1,:,:])+6*data[3,:,:]-7.5*data[1,:,:]))
    writefile(file,sarvi2,'sarvi2')
    del sarvi2
    sbl=data[4,:,:]-2.4*data[3,:,:]
    writefile(file,sbl,'sbl')
    del sbl
    sci=(data[5,:,:]-data[4,:,:])/(data[5,:,:]+data[4,:,:])
    writefile(file,sci,'sci')
    del sci
    savi2=data[4,:,:]/(data[3,:,:]+(1.22*0.03)*numpy.ones_like(data[1,:,:]))
    writefile(file,savi2,'savi2')
    del savi2
    slavi=data[4,:,:]/(data[3,:,:]+data[6,:,:])
    writefile(file,slavi,'slavi')
    del slavi
    sqrtirr=numpy.sqrt(data[4,:,:]/data[3,:,:])
    writefile(file,sqrtirr,'sqrtirr')
    del sqrtirr
    sbi=0.3037*data[1,:,:]+0.2793*data[2,:,:]+0.4743*data[3,:,:]+0.5585*data[4,:,:]+0.5082*data[5,:,:]+0.1863*data[6,:,:]
    writefile(file,sbi,'sbi')
    del sbi
    gvi=-0.2848*data[1,:,:]-0.2435*data[2,:,:]-0.5436*data[3,:,:]+0.7243*data[4,:,:]+0.084*data[5,:,:]-0.18*data[6,:,:]
    writefile(file,gvi,'gvi')
    del gvi
    wet = 0.1509 * data[1, :, :] +0.1973 * data[2, :, :] +0.3279 * data[3, :, :] + 0.3406 * data[4, :,:] - 0.7112 * data[5, :,:] - 0.4572 * data[6,:,:]
    writefile(file,wet,'wet')
    del wet
    tvi=numpy.sqrt(ndvi+0.5*numpy.ones_like(data[1,:,:]))
    writefile(file,tvi,'tvi')
    del tvi
    varigreen=(data[2,:,:]-data[3,:,:])/(data[2,:,:]+data[3,:,:]-data[1,:,:])
    writefile(file,varigreen,'varigreen')
    del varigreen
    wdvi=data[4,:,:]-1.22*data[3,:,:]
    writefile(file,wdvi,'wdvi')
    del wdvi
    wdrvi=(0.1*data[4,:,:]-data[4,:,:])/(0.1*data[4,:,:]+data[4,:,:])
    writefile(file,wdrvi,'wdrvi')
    del wdrvi
    return None
