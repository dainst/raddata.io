

import h5py,numpy,gdal,copy,utm
from scipy.ndimage import gaussian_filter1d
import cv2, numpy, gdalnumeric
from scipy.signal import bartlett,triang,gaussian


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
#Applies only for fil IO below (please see auxfun.py)
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


def schreibeBSQ(numpyarray,out): ##mal bitte nur fuer Hyspex nehmen!!!
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

def schreibeBSQsingle(numpyarray,out):
    format = "ENVI"
    driver = gdal.GetDriverByName(format)
    h=numpy.shape(numpyarray)
    dst_ds = driver.Create(out,h[1],h[0],1, gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(numpyarray)
    dst_ds=None

def schreibetif(numpyarray,out): ##mal bitte nur fuer Hyspex nehmen!!!
    format = "GTiff"
    driver = gdal.GetDriverByName( format )
    h=numpy.shape(numpyarray)
    dst_ds = driver.Create(out, h[2], h[1], h[0], gdal.GDT_Byte )
    ka=1
    while ka<=h[0]:
        hhh=ka-1
        dst_ds.GetRasterBand(ka).WriteArray(numpyarray[hhh,:,:])
        ka+=1
    dst_ds=None

def schreibetifsingle(numpyarray,out):
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    h=numpy.shape(numpyarray)
    dst_ds = driver.Create(out,h[1],h[0],1, gdal.GDT_Byte )
    dst_ds.GetRasterBand(1).WriteArray(numpyarray)
########## FILE IO End
#################
######################
  """  ---------------------
    Date                 : 01.2021
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""
def histeq(im,bins=256):
    im=numpy.nan_to_num(im)
    numpy.place(im,im<-99,0)
    numpy.place(im,im>9999999,0)
    im=numpy.nan_to_num(im)
    idx=numpy.where(im!=0)
    datah=numpy.percentile(im[idx].flatten(),[5,95],interpolation='nearest')
    image_histo,bins=numpy.histogram(im[idx].flatten(),bins,range=datah)
    imasch=numpy.interp(im.flatten(),bins[:-1],numpy.arange(0,256,1))
    imasch=imasch.reshape(im.shape)
    imasch=numpy.round(imasch)
    imasch=imasch.astype('uint8',casting='unsafe')
    return imasch

#so funzt es:
def recode(im,recode_classes=256):
    im = numpy.nan_to_num(im)
    numpy.place(im, im < -99, 0)
    #numpy.place(im, im > 9999999999, 0)
    im = numpy.nan_to_num(im)
    imshape = im.shape
    a=numpy.max(im)
    b=numpy.min(im)
    div=(a-b)/256
    im=im/div
    im=numpy.nan_to_num(im)
    return numpy.round(im)

def do_sectorbilder(im):
    nn=recode(im)
    unq=numpy.unique(nn)
    output=numpy.zeros_like(im)
    l=[]
    for j in enumerate(unq):
        idx=numpy.where(nn==j[1])
        output[idx]=recode(im[idx])
        w=numpy.zeros_like(output)
        w[idx]=output[idx]
        l.append(w)
    return output,l


#Routine zum Histeq scaling auf uint8 bilder
def d3d(imasch):
    shp=imasch.shape
    if len(shp)==2:
        L=histeq(imasch)#[0]
    if len(shp)==3:
        L=[]
        for j in numpy.arange(0,shp[0],1):
            a=histeq(imasch[j,:,:])
            L.append(a)
        L=numpy.asarray(L)
    return L


def shrp_res(rgb_bild,panbild):
    imasch=cv2.imread(rgb_bild)
    panimasch = cv2.imread(panbild)
    imsh=imasch.shape
    psh=panimasch.shape
    #print(psh,imsh)
    scl=cv2.resize(imasch,(psh[1],psh[0]),interpolation = cv2.INTER_CUBIC)#Dimensionen 1 und 0 vertauscht deswegen hatte es NICHT funktioniert
    data=cv2.cvtColor(scl,cv2.COLOR_BGR2HSV_FULL)
    print(data.shape)
    print(imasch.shape)
    pangray=cv2.cvtColor(panimasch,cv2.COLOR_BGR2HSV)
    print(pangray.shape)
    data[:,:,2]=panimasch[:,:,2]
    out=cv2.cvtColor(data,cv2.COLOR_HSV2BGR)
    jpnam=rgb_bild.split('.')[0]+'_shrp_result.jpg'
    print(jpnam)
    cv2.imwrite(jpnam,out)
    return None


    dst_ds=None

def do_headerkram(keyword,arrayss):
    arrayss=str(list(arrayss))
    arrayss=arrayss.replace("[","{")
    arrayss = arrayss.replace("]", "}")
    arrayss=keyword+'='+arrayss+"\n"
    return arrayss

class open_prisma(object):# Nur für L2D daten!!!!!!!!!!!!
    def __init__(self,fname):
        self.file=h5py.File(fname,'r')
        if '.' in fname:
            self.bname=fname.split('.')[0]
        else:
            self.bname=fname
        self.attrbute=list(self.file.attrs.keys())
        self.vswir=list(self.file["HDFEOS"]["SWATHS"]['PRS_L2D_HCO']["Data Fields"])
        self.pan=list(self.file["HDFEOS"]["SWATHS"]['PRS_L2D_PCO']["Data Fields"])
        for j in enumerate (self.attrbute):
            #print(j[1])
            setattr(self,j[1].replace(' ','_'), self.file.attrs.get(j[1])) #strings als attributnamen und dann dem Attribut einen Wert zuweisen
        for j in enumerate(self.vswir):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L2D_HCO']["Data Fields"].get(j[1]))
        for j in enumerate(self.pan):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L2D_PCO']["Data Fields"].get(j[1]))
        self.cubeswir=numpy.swapaxes(self.SWIR_Cube,0,1)[:-5,:,:]#30m, ersten 2 Kanäle immer Grütze # Achrung Hardcoded könnte man anhand der flags auch machen
        self.cubevnir=numpy.swapaxes(self.VNIR_Cube,0,1)[5:,:,:]#30m
        cubevnir_err=numpy.swapaxes(self.VNIR_PIXEL_L2_ERR_MATRIX,0,1)[5:,:,:]
        cubeswir_err=numpy.swapaxes(self.SWIR_PIXEL_L2_ERR_MATRIX,0,1)[:-5,:,:]
        self.cubepan=self.Cube#5m
        self.vnir_wv=self.List_Cw_Vnir[5:]#letzten beiden sind immer schrott
        self.swir_wv = self.List_Cw_Swir[:-5]
        self.sortvnir=numpy.argsort(self.vnir_wv)
        self.sortswir=numpy.argsort(self.swir_wv)
        self.vnir_fwhm = self.List_Fwhm_Vnir[5:]
        self.swir_fwhm = self.List_Fwhm_Swir[:-5]
        self.vnir_wv=self.vnir_wv[self.sortvnir]
        self.swir_wv = self.swir_wv[self.sortswir]
        self.vnir_fwhm = self.vnir_fwhm[self.sortvnir]
        self.swir_fwhm = self.swir_fwhm[self.sortswir]
        self.fwhmall=numpy.concatenate([self.vnir_fwhm,self.swir_fwhm])
        self.wvall = numpy.concatenate([self.vnir_wv, self.swir_wv])
        self.cubeswir=self.cubeswir[self.sortswir,:,:]
        self.cubevnir =self.cubevnir[self.sortvnir,:,:]
        div=numpy.nan_to_num(self.cubevnir[-1,:,:]/self.cubeswir[0,:,:])#Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall=numpy.concatenate([self.cubevnir,(self.cubeswir*div)],axis=0)
        self.errall=numpy.concatenate([cubevnir_err,cubeswir_err],axis=0)
        self.east=self.Product_ULcorner_easting
        self.north = self.Product_ULcorner_northing
        self.zone= self.Projection_Id
        self.ellips=self.Reference_Ellipsoid
        corlat=self.Product_LLcorner_lat
        if corlat>0:
            self.zzone='North'
        else:
            self.zzone="South"
        self.crstringhysp="map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",30,30,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n" #UTM Zonendefinition checken!!!!
        self.crstringpan = "map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",5,5,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n"
    def write_pan(self):
        schreibeBSQsingle(numpy.asarray(self.cubepan),self.bname+'_pan_cube')
        open(self.bname +'_pan_cube' + '.hdr', 'a').writelines(self.crstringpan)
    def write_all(self,data):
        schreibeBSQ(self.cubeall,self.bname+'_ref_cube')#
        schreibeBSQ(data,self.bname+'_ref_cube_polished')
        eins=do_headerkram('wavelength',self.wvall)
        zwo=do_headerkram('fwhm',self.fwhmall)
        g=open(self.bname+'_ref_cube'+'.hdr','a')
        h = open(self.bname +'_ref_cube_polished' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
        h.writelines(eins)
        h.writelines(zwo)
        h.writelines(self.crstringhysp)
    def pn_shrp(self,indices,bild,panbild):
        if len(indices)!=3:
            print("error")
        a= bild[indices[0],:,:]
        b = bild[indices[1], :, :]
        c= bild[indices[2], :, :]
        lbb=numpy.asarray([a,b,c])
        strechedrgb=d3d(lbb)
        strechedpan=d3d(panbild)
        clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(24,24))
        strechedrgb[0,:,:]=clahe.apply(strechedrgb[0,:,:])
        strechedrgb[1, :, :] = clahe.apply(strechedrgb[1, :, :])
        strechedrgb[2, :, :] = clahe.apply(strechedrgb[2, :, :])
        strechedpan=clahe.apply(strechedpan)
        print(strechedrgb.shape,strechedpan.shape)
        schreibetifsingle(strechedpan,self.bname+'_panbild')
        schreibetif(strechedrgb,self.bname+'_rgbbild')
        shrp_res(self.bname+'_rgbbild',self.bname+'_panbild')
    def spectral_polish(self):
        cpdat=copy.deepcopy(self.cubeall)
        dsh=cpdat.shape
        dneo=cpdat.reshape(dsh[0],dsh[1]*dsh[2])
        #print(dneo.shape)
        self.spec_pol=numpy.zeros_like(dneo)
        wvred=numpy.concatenate([self.wvall[0:56],self.wvall[62:74],self.wvall[82:]],axis=0)# mal 54 bis 68 probieren auszuschließen!!!!!
        for j in numpy.arange(0,dsh[1]*dsh[2],1):
            if numpy.sum(dneo[:,j])==0:
                continue
            spc=dneo[:,j]
            #print(spc.shape)
            sp1=numpy.concatenate([spc[0:56],spc[62:74],spc[82:]],axis=0)
            #print(self.wvall.shape,wvred.shape,sp1.shape)
            restor1=numpy.interp(self.wvall,wvred,sp1)
            d1=gaussian_filter1d(restor1,sigma=0.66)
            d3 = gaussian_filter1d(restor1, sigma=0.4)
            d2=gaussian_filter1d(restor1, sigma=1.66)
            self.spec_pol[:30,j]=d1[:30]
            self.spec_pol[30:85, j]=d2[30:85]
            self.spec_pol[85:180, j]=d1[85:180]
            self.spec_pol[180:, j] = d3[180:]
        self.spec_pol=self.spec_pol.reshape(dsh)
    def doit_swir(self):
        pass
    def doit_vnir(self):
        pass

def labelerrors2(chan,bild):
    l=numpy.unique(chan)
    dl=numpy.where(chan!=0)
    bcp=copy.deepcopy(bild)
    for j in enumerate(dl[0]):
        kaest=chan[dl[1][j[0]]-3:dl[1][j[0]]+3,j[1]-3:j[1]+3]
        bkaest=bild[dl[1][j[0]]-3:dl[1][j[0]]+3,j[1]-3:j[1]+3]
        zeugs=numpy.where(kaest!=0)
        if len(zeugs)==1:
            continue
        nn=numpy.nanmedian(bkaest)
        bcp[dl[0][j[0]],dl[1][j[0]]]=nn
    return bcp

#Mit flag noch zusaetzliches entstreiden anknippsen!
############Staenz et al. 2008 Destriping!
###################
##############################
def multall(dat,gains):
    g=0
    for i in gains:
        if i==numpy.nan:
            i=1
        dat[:,g]=dat[:,g]*i
        g+=1
    return dat

def gained(cubi,gains):
    dim=cubi.shape
    for o in range(1,dim[0],1):
        cubi_gained=multall(cubi[o,:,:],gains[o,:])
        cubi[o,:,:]=cubi_gained
    return cubi

def triang_normed(N):  #Dreicksfilter mit der Flaeche 1 und edge != 0
    if N%2 == 0:
        print("you need to input an odd number")
        exit(-1)
    filt=triang(N)
    filt=filt/((N+1)/2)
    return filt

def average_columns(data):
    dim=data.shape
    data=data/1.0
    daten=numpy.ndarray((dim[0],dim[2]),dtype='float')
    for i in range(0,dim[0],1):
        dat=data[i,:,:]
        dat=numpy.nanmedian(dat,axis=0)#median statt average um extremwerte nicht mit einzubeziehen
        daten[i,:]=dat
    return daten

def Im_cheating(data,filter_normed):
    dat=average_columns(data)
    sata=dat
    dim=dat.shape
    cache=copy.deepcopy(dat)
    for i in range(0,dim[0],1):
        dat=cache[i,:]
        dat=numpy.abs(numpy.convolve(dat,filter_normed,'same'))
        sata[i,:]=dat
    return sata
#Das hier ist das Destriping!
def stanze_hyperion(data):
    fido=triang_normed(33)#55 statt 11
    avg_dat=average_columns(data)
    avg_dat_filt=Im_cheating(data,fido)
    gains=abs(avg_dat_filt/avg_dat)
    numpy.place(gains,gains==numpy.nan,1)
    destriped=gained(data,gains)
    return destriped#,gains


#L2C Daten sind in Sensor-Geomterie, d.h. sie müssen noch georeferenziert werden -> Hier kann man destriping machen!!!!

class open_prismal2c(object):
    def __init__(self,fname):
        self.file=h5py.File(fname,'r')
        if '.' in fname:
            self.bname=fname.split('.')[0]
        else:
            self.bname=fname
        self.attrbute=list(self.file.attrs.keys())
        self.vswir=list(self.file["HDFEOS"]["SWATHS"]['PRS_L2C_HCO']["Data Fields"])
        self.pan=list(self.file["HDFEOS"]["SWATHS"]['PRS_L2C_PCO']["Data Fields"])
        for j in enumerate (self.attrbute):
            #print(j[1])
            setattr(self,j[1].replace(' ','_'), self.file.attrs.get(j[1])) #strings als attributnamen und dann dem Attribut einen Wert zuweisen
        for j in enumerate(self.vswir):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L2C_HCO']["Data Fields"].get(j[1]))
        for j in enumerate(self.pan):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L2C_PCO']["Data Fields"].get(j[1]))
        self.cubeswir=numpy.swapaxes(self.SWIR_Cube,0,1)[:-5,:,:]#30m, ersten 2 Kanäle immer Grütze # Achrung Hardcoded könnte man anhand der flags auch machen
        self.cubevnir=numpy.swapaxes(self.VNIR_Cube,0,1)[5:,:,:]#30m
        cubevnir_err=numpy.swapaxes(self.VNIR_PIXEL_L2_ERR_MATRIX,0,1)[5:,:,:]
        cubeswir_err=numpy.swapaxes(self.SWIR_PIXEL_L2_ERR_MATRIX,0,1)[:-5,:,:]
        self.cubepan=self.Cube#5m
        self.vnir_wv=self.List_Cw_Vnir[5:]#letzten beiden sind immer schrott
        self.swir_wv = self.List_Cw_Swir[:-5]
        self.sortvnir=numpy.argsort(self.vnir_wv)
        self.sortswir=numpy.argsort(self.swir_wv)
        self.vnir_fwhm = self.List_Fwhm_Vnir[5:]
        self.swir_fwhm = self.List_Fwhm_Swir[:-5]
        self.vnir_wv=self.vnir_wv[self.sortvnir]
        self.swir_wv = self.swir_wv[self.sortswir]
        self.vnir_fwhm = self.vnir_fwhm[self.sortvnir]
        self.swir_fwhm = self.swir_fwhm[self.sortswir]
        self.fwhmall=numpy.concatenate([self.vnir_fwhm,self.swir_fwhm])
        self.wvall = numpy.concatenate([self.vnir_wv, self.swir_wv])
        self.cubeswir=self.cubeswir[self.sortswir,:,:]
        self.cubevnir =self.cubevnir[self.sortvnir,:,:]
        div=numpy.nan_to_num(self.cubevnir[-1,:,:]/self.cubeswir[0,:,:])#Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall=numpy.concatenate([self.cubevnir,(self.cubeswir*div)],axis=0)
        self.errall=numpy.concatenate([cubevnir_err,cubeswir_err],axis=0)
        #self.east=self.Product_ULcorner_easting
        #self.north = self.Product_ULcorner_northing
        lon=self.Product_ULcorner_long
        lat=self.Product_ULcorner_lat
        self.zone=utm.latlon_to_zone_number(lat,lon)
        self.east,self.north,self.zone,self.letter=utm.from_latlon(lat,lon)
        #self.zone= self.Projection_Id
        self.ellips='WGS84'#self.Reference_Ellipsoid
        corlat=self.Product_LLcorner_lat
        if corlat>0:
            self.zzone='North'
        else:
            self.zzone="South"
        self.crstringhysp="map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",30,30,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n" #UTM Zonendefinition checken!!!!
        self.crstringpan = "map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",5,5,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n"
    def write_pan(self):
        schreibeBSQsingle(numpy.asarray(self.cubepan),self.bname+'_pan_cube')
        open(self.bname +'_pan_cube' + '.hdr', 'a').writelines(self.crstringpan)
    def write_all(self,data):
        schreibeBSQ(self.cubeall,self.bname+'_ref_cube')#
        schreibeBSQ(data,self.bname+'_ref_cube_polished')
        eins=do_headerkram('wavelength',self.wvall)
        zwo=do_headerkram('fwhm',self.fwhmall)
        g=open(self.bname+'_ref_cube'+'.hdr','a')
        h = open(self.bname +'_ref_cube_polished' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
        h.writelines(eins)
        h.writelines(zwo)
        h.writelines(self.crstringhysp)
    def pn_shrp(self,indices,bild,panbild):
        if len(indices)!=3:
            print("error")
        a= bild[indices[0],:,:]
        b = bild[indices[1], :, :]
        c= bild[indices[2], :, :]
        lbb=numpy.asarray([a,b,c])
        strechedrgb=d3d(lbb)
        strechedpan=d3d(panbild)
        print(strechedrgb.shape,strechedpan.shape)
        schreibetifsingle(strechedpan,self.bname+'_panbild')
        schreibetif(strechedrgb,self.bname+'_rgbbild')
        shrp_res(self.bname+'_rgbbild',self.bname+'_panbild')
    def spectral_polish(self):
        cpdat=copy.deepcopy(self.cubeall)
        dsh=cpdat.shape
        dneo=cpdat.reshape(dsh[0],dsh[1]*dsh[2])
        #print(dneo.shape)
        self.spec_pol=numpy.zeros_like(dneo)
        wvred=numpy.concatenate([self.wvall[0:56],self.wvall[62:74],self.wvall[82:]],axis=0)# mal 54 bis 68 probieren auszuschließen!!!!!
        for j in numpy.arange(0,dsh[1]*dsh[2],1):
            if numpy.sum(dneo[:,j])==0:
                continue
            spc=dneo[:,j]
            #print(spc.shape)
            sp1=numpy.concatenate([spc[0:56],spc[62:74],spc[82:]],axis=0)
            #print(self.wvall.shape,wvred.shape,sp1.shape)
            restor1=numpy.interp(self.wvall,wvred,sp1)
            d1=gaussian_filter1d(restor1,sigma=0.66)
            d3 = gaussian_filter1d(restor1, sigma=0.4)
            d2=gaussian_filter1d(restor1, sigma=1.66)
            self.spec_pol[:30,j]=d1[:30]
            self.spec_pol[30:85, j]=d2[30:85]
            self.spec_pol[85:180, j]=d1[85:180]
            self.spec_pol[180:, j] = d3[180:]
        self.spec_pol=self.spec_pol.reshape(dsh)
    def doit_swir(self):
        pass
    def doit_vnir(self):
        pass

#Staenz destriping implementiert, georef, fehlt!
class open_prismal2b(object):
    def __init__(self,fname):
        self.file=h5py.File(fname,'r')
        if '.' in fname:
            self.bname=fname.split('.')[0]
        else:
            self.bname=fname
        self.attrbute=list(self.file.attrs.keys())
        self.vswir=list(self.file["HDFEOS"]["SWATHS"]['PRS_L2B_HCO']["Data Fields"])
        self.pan=list(self.file["HDFEOS"]["SWATHS"]['PRS_L2B_PCO']["Data Fields"])
        for j in enumerate (self.attrbute):
            #print(j[1])
            setattr(self,j[1].replace(' ','_'), self.file.attrs.get(j[1])) #strings als attributnamen und dann dem Attribut einen Wert zuweisen
        for j in enumerate(self.vswir):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L2B_HCO']["Data Fields"].get(j[1]))
        for j in enumerate(self.pan):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L2B_PCO']["Data Fields"].get(j[1]))
        self.cubeswir=numpy.swapaxes(self.SWIR_Cube,0,1)[:-5,:,:]#30m, ersten 2 Kanäle immer Grütze # Achrung Hardcoded könnte man anhand der flags auch machen
        self.cubevnir=numpy.swapaxes(self.VNIR_Cube,0,1)[5:,:,:]#30m
        cubevnir_err=numpy.swapaxes(self.VNIR_PIXEL_L2_ERR_MATRIX,0,1)[5:,:,:]
        cubeswir_err=numpy.swapaxes(self.SWIR_PIXEL_L2_ERR_MATRIX,0,1)[:-5,:,:]
        self.cubepan=self.Cube#5m
        self.vnir_wv=self.List_Cw_Vnir[5:]#letzten beiden sind immer schrott
        self.swir_wv = self.List_Cw_Swir[:-5]
        self.sortvnir=numpy.argsort(self.vnir_wv)
        self.sortswir=numpy.argsort(self.swir_wv)
        self.vnir_fwhm = self.List_Fwhm_Vnir[5:]
        self.swir_fwhm = self.List_Fwhm_Swir[:-5]
        self.vnir_wv=self.vnir_wv[self.sortvnir]
        self.swir_wv = self.swir_wv[self.sortswir]
        self.vnir_fwhm = self.vnir_fwhm[self.sortvnir]
        self.swir_fwhm = self.swir_fwhm[self.sortswir]
        self.fwhmall=numpy.concatenate([self.vnir_fwhm,self.swir_fwhm])
        self.wvall = numpy.concatenate([self.vnir_wv, self.swir_wv])
        self.cubeswir=self.cubeswir[self.sortswir,:,:]
        self.cubevnir =self.cubevnir[self.sortvnir,:,:]
        div=numpy.nan_to_num(self.cubevnir[-1,:,:]/self.cubeswir[0,:,:])#Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall=numpy.concatenate([self.cubevnir,(self.cubeswir*div)],axis=0)
        self.errall=numpy.concatenate([cubevnir_err,cubeswir_err],axis=0)
        #self.east=self.Product_ULcorner_easting
        #self.north = self.Product_ULcorner_northing
        lon=self.Product_ULcorner_long
        lat=self.Product_ULcorner_lat
        self.zone=utm.latlon_to_zone_number(lat,lon)
        self.east,self.north,self.zone,self.letter=utm.from_latlon(lat,lon)
        #self.zone= self.Projection_Id
        self.ellips='WGS84'#self.Reference_Ellipsoid
        corlat=self.Product_LLcorner_lat
        if corlat>0:
            self.zzone='North'
        else:
            self.zzone="South"
        self.crstringhysp="map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",30,30,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n" #UTM Zonendefinition checken!!!!
        self.crstringpan = "map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",5,5,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n"
    def write_pan(self):
        schreibeBSQsingle(numpy.asarray(self.cubepan),self.bname+'_pan_cube')
        open(self.bname +'_pan_cube' + '.hdr', 'a').writelines(self.crstringpan)
    def write_all(self,data):
        schreibeBSQ(self.cubeall,self.bname+'_ref_cube')#
        schreibeBSQ(data,self.bname+'_ref_cube_polished')
        eins=do_headerkram('wavelength',self.wvall)
        zwo=do_headerkram('fwhm',self.fwhmall)
        g=open(self.bname+'_ref_cube'+'.hdr','a')
        h = open(self.bname +'_ref_cube_polished' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
        h.writelines(eins)
        h.writelines(zwo)
        h.writelines(self.crstringhysp)
    def pn_shrp(self,bild,panbild,indices=[190,59,20]):
        if len(indices)!=3:
            print("error")
        a= bild[indices[0],:,:]
        b = bild[indices[1], :, :]
        c= bild[indices[2], :, :]
        lbb=numpy.asarray([a,b,c])
        strechedrgb=d3d(lbb)
        strechedpan=d3d(panbild)
        print(strechedrgb.shape,strechedpan.shape)
        schreibetifsingle(strechedpan,self.bname+'_panbild')
        schreibetif(strechedrgb,self.bname+'_rgbbild')
        shrp_res(self.bname+'_rgbbild',self.bname+'_panbild')
    def spectral_polish(self):
        cpdat=copy.deepcopy(self.cubeall)
        dsh=cpdat.shape
        dneo=cpdat.reshape(dsh[0],dsh[1]*dsh[2])
        #print(dneo.shape)
        self.spec_pol=numpy.zeros_like(dneo)
        wvred=numpy.concatenate([self.wvall[0:56],self.wvall[62:74],self.wvall[82:]],axis=0)# mal 54 bis 68 probieren auszuschließen!!!!!
        for j in numpy.arange(0,dsh[1]*dsh[2],1):
            if numpy.sum(dneo[:,j])==0:
                continue
            spc=dneo[:,j]
            #print(spc.shape)
            sp1=numpy.concatenate([spc[0:56],spc[62:74],spc[82:]],axis=0)
            #print(self.wvall.shape,wvred.shape,sp1.shape)
            restor1=numpy.interp(self.wvall,wvred,sp1)
            d1=gaussian_filter1d(restor1,sigma=0.66)
            d3 = gaussian_filter1d(restor1, sigma=0.4)
            d2=gaussian_filter1d(restor1, sigma=1.66)
            self.spec_pol[:30,j]=d1[:30]
            self.spec_pol[30:85, j]=d2[30:85]
            self.spec_pol[85:180, j]=d1[85:180]
            self.spec_pol[180:, j] = d3[180:]
        self.spec_pol=self.spec_pol.reshape(dsh)
    def vnir_destripe(self):
        self.cubevnir_destr=self.cubevnir
        self.cubevnir_destr=numpy.rot90(self.cubevnir,k=1,axes=(1,2))
        self.cubevnir_destr=stanze_hyperion(self.cubevnir_destr)
        self.cubevnir_destr=numpy.rot90(self.cubevnir_destr,k=-1,axes=(1,2))
    def swir_destripe(self):
        self.cubeswir_destr = self.cubeswir
        self.cubeswir_destr = numpy.rot90(self.cubeswir, k=1, axes=(1, 2))
        self.cubeswir_destr = stanze_hyperion(self.cubeswir_destr)
        self.cubeswir_destr = numpy.rot90(self.cubeswir_destr, k=-1, axes=(1, 2))
    def do_destr_first(self):
        self.vnir_destripe()
        self.swir_destripe()
        div = numpy.nan_to_num(self.cubevnir_destr[-1, :, :] / self.cubeswir_destr[0, :, :])  # Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall_destripe = numpy.concatenate([self.cubevnir_destr,(self.cubeswir_destr*div)], axis=0)
        schreibeBSQ(self.cubeall_destripe,self.bname + '_ref_cube_destriped')  #
        eins = do_headerkram('wavelength', self.wvall)
        zwo = do_headerkram('fwhm', self.fwhmall)
        g = open(self.bname + '_ref_cube_destriped' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
#Siehe L2B jedoch nun für L2A
class open_prismal1(object):
    def __init__(self,fname):
        self.file=h5py.File(fname,'r')
        if '.' in fname:
            self.bname=fname.split('.')[0]
        else:
            self.bname=fname
        self.attrbute=list(self.file.attrs.keys())
        self.vswir=list(self.file["HDFEOS"]["SWATHS"]['PRS_L1_HCO']["Data Fields"])
        self.pan=list(self.file["HDFEOS"]["SWATHS"]['PRS_L1_PCO']["Data Fields"])
        for j in enumerate (self.attrbute):
            #print(j[1])
            setattr(self,j[1].replace(' ','_'), self.file.attrs.get(j[1])) #strings als attributnamen und dann dem Attribut einen Wert zuweisen
        for j in enumerate(self.vswir):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L1_HCO']["Data Fields"].get(j[1]))
        for j in enumerate(self.pan):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L1_PCO']["Data Fields"].get(j[1]))
        self.cubeswir=numpy.swapaxes(self.SWIR_Cube,0,1)[:-5,:,:]#30m, ersten 2 Kanäle immer Grütze # Achrung Hardcoded könnte man anhand der flags auch machen
        self.cubevnir=numpy.swapaxes(self.VNIR_Cube,0,1)[5:,:,:]#30m
        cubevnir_err=numpy.swapaxes(self.VNIR_PIXEL_SAT_ERR_MATRIX,0,1)[5:,:,:]
        cubeswir_err=numpy.swapaxes(self.SWIR_PIXEL_SAT_ERR_MATRIX,0,1)[:-5,:,:]
        self.cubepan=self.Cube#5m
        self.vnir_wv=self.List_Cw_Vnir[5:]#letzten beiden sind immer schrott
        self.swir_wv = self.List_Cw_Swir[:-5]
        self.sortvnir=numpy.argsort(self.vnir_wv)
        self.sortswir=numpy.argsort(self.swir_wv)
        self.vnir_fwhm = self.List_Fwhm_Vnir[5:]
        self.swir_fwhm = self.List_Fwhm_Swir[:-5]
        self.vnir_wv=self.vnir_wv[self.sortvnir]
        self.swir_wv = self.swir_wv[self.sortswir]
        self.vnir_fwhm = self.vnir_fwhm[self.sortvnir]
        self.swir_fwhm = self.swir_fwhm[self.sortswir]
        self.fwhmall=numpy.concatenate([self.vnir_fwhm,self.swir_fwhm])
        self.wvall = numpy.concatenate([self.vnir_wv, self.swir_wv])
        self.cubeswir=self.cubeswir[self.sortswir,:,:]
        self.cubevnir =self.cubevnir[self.sortvnir,:,:]
        div=numpy.nan_to_num(self.cubevnir[-1,:,:]/self.cubeswir[0,:,:])#Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall=numpy.concatenate([self.cubevnir,(self.cubeswir*div)],axis=0)
        self.errall=numpy.concatenate([cubevnir_err,cubeswir_err],axis=0)
        #self.east=self.Product_ULcorner_easting
        #self.north = self.Product_ULcorner_northing
        lon=" "#self.Product_ULcorner_long
        lat=" "#self.Product_ULcorner_lat
        self.zone=' '#+utm.latlon_to_zone_number(lat,lon)
        self.east,self.north,self.zone,self.letter=' ',' ',' ',' '#utm.from_latlon(lat,lon)
        #self.zone= self.Projection_Id
        self.ellips='WGS84'#self.Reference_Ellipsoid
        corlat='a'#self.Product_LLcorner_lat
        """if corlat>0:
            self.zzone='North'
        else:
            self.zzone="South"""
        self.zzone = "South"
        self.crstringhysp="\n"#"map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",30,30,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n" #UTM Zonendefinition checken!!!!
        self.crstringpan ="\n"# "map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",5,5,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n"
    def write_pan(self):
        schreibeBSQsingle(numpy.asarray(self.cubepan),self.bname+'_pan_cube')
        open(self.bname +'_pan_cube' + '.hdr', 'a').writelines(self.crstringpan)
    def write_all(self,data):
        schreibeBSQ(self.cubeall,self.bname+'_ref_cube')#Error Matrizen des Sensors sollte man noch Anwenden sonst sieht es scheiße aus!
        schreibeBSQ(data,self.bname+'_ref_cube_polished')
        eins=do_headerkram('wavelength',self.wvall)
        zwo=do_headerkram('fwhm',self.fwhmall)
        g=open(self.bname+'_ref_cube'+'.hdr','a')
        h = open(self.bname +'_ref_cube_polished' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
        h.writelines(eins)
        h.writelines(zwo)
        h.writelines(self.crstringhysp)
    def pn_shrp(self,bild,panbild,indices=[190,59,20]):
        if len(indices)!=3:
            print("error")
        a= bild[indices[0],:,:]
        b = bild[indices[1], :, :]
        c= bild[indices[2], :, :]
        lbb=numpy.asarray([a,b,c])
        strechedrgb=d3d(lbb)
        strechedpan=d3d(panbild)
        print(strechedrgb.shape,strechedpan.shape)
        schreibetifsingle(strechedpan,self.bname+'_panbild')
        schreibetif(strechedrgb,self.bname+'_rgbbild')
        shrp_res(self.bname+'_rgbbild',self.bname+'_panbild')
    def spectral_polish(self):
        cpdat=copy.deepcopy(self.cubeall)
        dsh=cpdat.shape
        dneo=cpdat.reshape(dsh[0],dsh[1]*dsh[2])
        #print(dneo.shape)
        self.spec_pol=numpy.zeros_like(dneo)
        wvred=numpy.concatenate([self.wvall[0:56],self.wvall[62:74],self.wvall[82:]],axis=0)# mal 54 bis 68 probieren auszuschließen!!!!!
        for j in numpy.arange(0,dsh[1]*dsh[2],1):
            if numpy.sum(dneo[:,j])==0:
                continue
            spc=dneo[:,j]
            #print(spc.shape)
            sp1=numpy.concatenate([spc[0:56],spc[62:74],spc[82:]],axis=0)
            #print(self.wvall.shape,wvred.shape,sp1.shape)
            restor1=numpy.interp(self.wvall,wvred,sp1)
            d1=gaussian_filter1d(restor1,sigma=0.66)
            d3 = gaussian_filter1d(restor1, sigma=0.4)
            d2=gaussian_filter1d(restor1, sigma=1.66)
            self.spec_pol[:30,j]=d1[:30]
            self.spec_pol[30:85, j]=d2[30:85]
            self.spec_pol[85:180, j]=d1[85:180]
            self.spec_pol[180:, j] = d3[180:]
        self.spec_pol=self.spec_pol.reshape(dsh)
    def vnir_destripe(self):
        self.cubevnir_destr=self.cubevnir
        self.cubevnir_destr=numpy.rot90(self.cubevnir,k=1,axes=(1,2))
        self.cubevnir_destr=stanze_hyperion(self.cubevnir_destr)
        self.cubevnir_destr=numpy.rot90(self.cubevnir_destr,k=-1,axes=(1,2))
    def swir_destripe(self):
        self.cubeswir_destr = self.cubeswir
        self.cubeswir_destr = numpy.rot90(self.cubeswir, k=1, axes=(1, 2))
        self.cubeswir_destr = stanze_hyperion(self.cubeswir_destr)
        self.cubeswir_destr = numpy.rot90(self.cubeswir_destr, k=-1, axes=(1, 2))
    def do_destr_first(self):
        self.vnir_destripe()
        self.swir_destripe()
        div = numpy.nan_to_num(self.cubevnir_destr[-1, :, :] / self.cubeswir_destr[0, :, :])  # Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall_destripe = numpy.concatenate([self.cubevnir_destr,(self.cubeswir_destr*div)], axis=0)
        schreibeBSQ(self.cubeall_destripe,self.bname + '_ref_cube_destriped')  # Error Matrizen des Sensors sollte man noch Anwenden sonst sieht es scheiße aus!
        eins = do_headerkram('wavelength', self.wvall)
        zwo = do_headerkram('fwhm', self.fwhmall)
        g = open(self.bname + '_ref_cube_destriped' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
#Das sind die Verstreiften Daten; HCO Daten sind unverstreift!!!
#Siehe L2B jedoch nun für L2A
class open_prismal1_hrc(object):
    def __init__(self,fname):
        self.file=h5py.File(fname,'r')
        if '.' in fname:
            self.bname=fname.split('.')[0]
        else:
            self.bname=fname
        self.attrbute=list(self.file.attrs.keys())
        self.vswir=list(self.file["HDFEOS"]["SWATHS"]['PRS_L1_HRC']["Data Fields"])
        self.pan=list(self.file["HDFEOS"]["SWATHS"]['PRS_L1_PRC']["Data Fields"])
        for j in enumerate (self.attrbute):
            #print(j[1])
            setattr(self,j[1].replace(' ','_'), self.file.attrs.get(j[1])) #strings als attributnamen und dann dem Attribut einen Wert zuweisen
        for j in enumerate(self.vswir):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L1_HRC']["Data Fields"].get(j[1]))
        for j in enumerate(self.pan):
            setattr(self, j[1].replace(' ', '_'), self.file["HDFEOS"]["SWATHS"]['PRS_L1_PRC']["Data Fields"].get(j[1]))
        self.cubeswir=numpy.swapaxes(self.SWIR_Cube,0,1)[:-5,:,:]#30m, ersten 2 Kanäle immer Grütze # Achrung Hardcoded könnte man anhand der flags auch machen
        self.cubevnir=numpy.swapaxes(self.VNIR_Cube,0,1)[5:,:,:]#30m
        cubevnir_err=numpy.swapaxes(self.VNIR_PIXEL_SAT_ERR_MATRIX,0,1)[5:,:,:]
        cubeswir_err=numpy.swapaxes(self.SWIR_PIXEL_SAT_ERR_MATRIX,0,1)[:-5,:,:]
        self.cubepan=self.Cube#5m
        self.vnir_wv=self.List_Cw_Vnir[5:]#letzten beiden sind immer schrott
        self.swir_wv = self.List_Cw_Swir[:-5]
        self.sortvnir=numpy.argsort(self.vnir_wv)
        self.sortswir=numpy.argsort(self.swir_wv)
        self.vnir_fwhm = self.List_Fwhm_Vnir[5:]
        self.swir_fwhm = self.List_Fwhm_Swir[:-5]
        self.vnir_wv=self.vnir_wv[self.sortvnir]
        self.swir_wv = self.swir_wv[self.sortswir]
        self.vnir_fwhm = self.vnir_fwhm[self.sortvnir]
        self.swir_fwhm = self.swir_fwhm[self.sortswir]
        self.fwhmall=numpy.concatenate([self.vnir_fwhm,self.swir_fwhm])
        self.wvall = numpy.concatenate([self.vnir_wv, self.swir_wv])
        self.cubeswir=self.cubeswir[self.sortswir,:,:]
        self.cubevnir =self.cubevnir[self.sortvnir,:,:]
        div=numpy.nan_to_num(self.cubevnir[-1,:,:]/self.cubeswir[0,:,:])#Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall=numpy.concatenate([self.cubevnir,(self.cubeswir*div)],axis=0)
        self.errall=numpy.concatenate([cubevnir_err,cubeswir_err],axis=0)
        #self.east=self.Product_ULcorner_easting
        #self.north = self.Product_ULcorner_northing
        lon=" "#self.Product_ULcorner_long
        lat=" "#self.Product_ULcorner_lat
        self.zone=' '#+utm.latlon_to_zone_number(lat,lon)
        self.east,self.north,self.zone,self.letter=' ',' ',' ',' '#utm.from_latlon(lat,lon)
        #self.zone= self.Projection_Id
        self.ellips='WGS84'#self.Reference_Ellipsoid
        corlat='a'#self.Product_LLcorner_lat
        """if corlat>0:
            self.zzone='North'
        else:
            self.zzone="South"""
        self.zzone = "South"
        self.crstringhysp="\n"#"map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",30,30,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n" #UTM Zonendefinition checken!!!!
        self.crstringpan ="\n"# "map info= {UTM,1,1,"+str(self.east)+','+str(self.north)+",5,5,"+str(int(self.zone))+","+self.zzone+", WGS84}"+"\n"
    def write_pan(self):
        schreibeBSQsingle(numpy.asarray(self.cubepan),self.bname+'_pan_cube')
        open(self.bname +'_pan_cube' + '.hdr', 'a').writelines(self.crstringpan)
    def write_all(self,data):
        schreibeBSQ(self.cubeall,self.bname+'_ref_cube')#Error Matrizen des Sensors sollte man noch Anwenden sonst sieht es scheiße aus!
        schreibeBSQ(data,self.bname+'_ref_cube_polished')
        eins=do_headerkram('wavelength',self.wvall)
        zwo=do_headerkram('fwhm',self.fwhmall)
        g=open(self.bname+'_ref_cube'+'.hdr','a')
        h = open(self.bname +'_ref_cube_polished' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
        h.writelines(eins)
        h.writelines(zwo)
        h.writelines(self.crstringhysp)
    def pn_shrp(self,bild,panbild,indices=[190,59,20]):
        if len(indices)!=3:
            print("error")
        a= bild[indices[0],:,:]
        b = bild[indices[1], :, :]
        c= bild[indices[2], :, :]
        lbb=numpy.asarray([a,b,c])
        strechedrgb=d3d(lbb)
        strechedpan=d3d(panbild)
        print(strechedrgb.shape,strechedpan.shape)
        schreibetifsingle(strechedpan,self.bname+'_panbild')
        schreibetif(strechedrgb,self.bname+'_rgbbild')
        shrp_res(self.bname+'_rgbbild',self.bname+'_panbild')
    def spectral_polish(self):
        cpdat=copy.deepcopy(self.cubeall)
        dsh=cpdat.shape
        dneo=cpdat.reshape(dsh[0],dsh[1]*dsh[2])
        #print(dneo.shape)
        self.spec_pol=numpy.zeros_like(dneo)
        wvred=numpy.concatenate([self.wvall[0:56],self.wvall[62:74],self.wvall[82:]],axis=0)# mal 54 bis 68 probieren auszuschließen!!!!!
        for j in numpy.arange(0,dsh[1]*dsh[2],1):
            if numpy.sum(dneo[:,j])==0:
                continue
            spc=dneo[:,j]
            #print(spc.shape)
            sp1=numpy.concatenate([spc[0:56],spc[62:74],spc[82:]],axis=0)
            #print(self.wvall.shape,wvred.shape,sp1.shape)
            restor1=numpy.interp(self.wvall,wvred,sp1)
            d1=gaussian_filter1d(restor1,sigma=0.66)
            d3 = gaussian_filter1d(restor1, sigma=0.4)
            d2=gaussian_filter1d(restor1, sigma=1.66)
            self.spec_pol[:30,j]=d1[:30]
            self.spec_pol[30:85, j]=d2[30:85]
            self.spec_pol[85:180, j]=d1[85:180]
            self.spec_pol[180:, j] = d3[180:]
        self.spec_pol=self.spec_pol.reshape(dsh)
    def vnir_destripe(self):
        self.cubevnir_destr=self.cubevnir
        self.cubevnir_destr=numpy.rot90(self.cubevnir,k=1,axes=(1,2))
        self.cubevnir_destr=stanze_hyperion(self.cubevnir_destr)
        self.cubevnir_destr=numpy.rot90(self.cubevnir_destr,k=-1,axes=(1,2))
    def swir_destripe(self):
        self.cubeswir_destr = self.cubeswir
        self.cubeswir_destr = numpy.rot90(self.cubeswir, k=1, axes=(1, 2))
        self.cubeswir_destr = stanze_hyperion(self.cubeswir_destr)
        self.cubeswir_destr = numpy.rot90(self.cubeswir_destr, k=-1, axes=(1, 2))
    def do_destr_first(self):
        self.vnir_destripe()
        self.swir_destripe()
        div = numpy.nan_to_num(self.cubevnir_destr[-1, :, :] / self.cubeswir_destr[0, :, :])  # Factor for SWIR on VNIR Niveau Adjustment
        self.cubeall_destripe = numpy.concatenate([self.cubevnir_destr,(self.cubeswir_destr*div)], axis=0)
        schreibeBSQ(self.cubeall_destripe,self.bname + '_ref_cube_destriped')  # 
        eins = do_headerkram('wavelength', self.wvall)
        zwo = do_headerkram('fwhm', self.fwhmall)
        g = open(self.bname + '_ref_cube_destriped' + '.hdr', 'a')
        g.writelines(eins)
        g.writelines(zwo)
        g.writelines(self.crstringhysp)
