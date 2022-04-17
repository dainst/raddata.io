import subprocess, auxfun
import os, sys, shlex, shutil, gdalnumeric, numpy
import glob
import auxfun as axf
import ratios_ls5_reduced
import ratios_ls7_reduced
import ratios_ls8_reduced
import sent2_ratios_reduced
import ratios_ls8_dai, ratios_ls7u5_dai, sent2_ratios_dai
import heda
import zipfile, tarfile

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
LANDSAT9PREFIX = 'LC09'
LANDSAT9FOLDERPREFIX = 'LS9_L02'
LANDSAT8PREFIX = 'LC08'
LANDSAT8FOLDERPREFIX = 'LS8_L02'
LANDSAT7PREFIX = 'LE07'
LANDSAT7FOLDERPREFIX = 'LS7_L02'
LANDSAT5PREFIX = 'LT05'
LANDSAT5FOLDERPREFIX = 'LS5_L02'
LANDSAT4PREFIX = 'LT04'
LANDSAT4FOLDERPREFIX = 'LS4_L02'
LSentinel2PREFIX = 'S2'
Sentinel2FOLDERPREFIX = 'S2_L02'
SENSORPREFIX = [LANDSAT4PREFIX, LANDSAT5PREFIX, LANDSAT7PREFIX, LANDSAT8PREFIX, LANDSAT9PREFIX, LSentinel2PREFIX]
L2FOLDERPREFIX = [LANDSAT4FOLDERPREFIX, LANDSAT5FOLDERPREFIX, LANDSAT7FOLDERPREFIX, LANDSAT8FOLDERPREFIX,
                  LANDSAT9FOLDERPREFIX, Sentinel2FOLDERPREFIX]

LANDSATREFLECTANCESUFFIX = "srenv_stack_file"
SENTINEL220MREFLECTANCESUFFIX = "_20m_stack_file"
SENTINEL210MREFLECTANCESUFFIX = "_10m_stack_file"
SENTINEL260MREFLECTANCESUFFIX = "_60m_stack_file"


# for python gt 3.5
# patternsearchrecursive

class SensorCalculationsStatic:
    @staticmethod
    def stack_all(folder):
        folderitems = FolderFunctions.subdirs(folder)
        sensordict = FolderFunctions.match_sensorpaths(folderitems)
        for sensor in sensordict.keys():
            if sensor == 'S2_L02':
                stack(sensordict[sensor])
            if sensor == 'LS4_L02':
                stackls4_7(sensordict[sensor])
            if sensor == 'LS7_L02':
                stackls4_7(sensordict[sensor])
            if sensor == 'LS5_L02':
                stackls4_7(sensordict[sensor])
            if sensor == 'LS8_L02':
                stackls8(sensordict[sensor])
            if sensor == 'LS9_L02':
                stackls8(sensordict[sensor])
        return sensordict

    @staticmethod
    def ratios(folder, flagg):
        folderitems = FolderFunctions.subdirs(folder)
        sensordict = FolderFunctions.match_sensorpaths(folderitems)
        for sensor in sensordict.keys():
            if sensor == 'S2_L02':
                stack_analyzeSen2_reduced(sensordict[sensor], flag=flagg, flagm=[20, 10])
                print('Finished S2 Ratios')
            if sensor == 'LS4_L02':
                stack_analyzeS5_reduced(sensordict[sensor], flag=flagg)
                print('Finished LS4 Ratios')
            if sensor == 'LS7_L02':
                stack_analyzeS7_reduced(sensordict[sensor], flag=flagg)
                print('Finished LS7 Ratios')
            if sensor == 'LS5_L02':
                stack_analyzeS5_reduced(sensordict[sensor], flag=flagg)
                print('Finished LS5 Ratios')
            if sensor == 'LS8_L02':
                stack_analyzeS8_reduced(sensordict[sensor], flag=flagg)
                print('Finished LS8 Ratios')
            if sensor == 'LS9_L02':
                stack_analyzeS8_reduced(sensordict[sensor], flag=flagg)
                print('Finished LS9 Ratios')
        return sensordict


class FolderFunctions:
    @staticmethod
    def sort_compressed(folder):
        folderlist = os.listdir(folder)
        for name in folderlist:
            if os.path.isfile(f'{folder}/{name}'):
                for sensoridx, sensor in enumerate(SENSORPREFIX):
                    try:
                        if name.startswith(sensor):
                            os.mkdir(f'{folder}/{L2FOLDERPREFIX[sensoridx]}')
                            os.rename(f'{folder}/{name}', f'{folder}/{L2FOLDERPREFIX[sensoridx]}/{name}')
                    except Exception:
                        try:
                            os.rename(f'{folder}/{name}', f'{folder}/{L2FOLDERPREFIX[sensoridx]}/{name}')
                        except Exception:
                            print(
                                f'Cant move file{name} somethings wrong with the destination or premissions and directory {L2FOLDERPREFIX[sensoridx]} already exists!')

    @staticmethod
    def sort_uncompressed(folder):
        folderlist = os.listdir(folder)
        for name in folderlist:
            if os.path.isdir(f'{folder}/{name}'):
                for sensoridx, sensor in enumerate(SENSORPREFIX):
                    try:
                        if name.startswith(sensor):
                            os.mkdir(f'{folder}/{L2FOLDERPREFIX[sensoridx]}')
                            os.renames(f'{folder}/{name}', f'{folder}/{L2FOLDERPREFIX[sensoridx]}/{name}')
                    except Exception:
                        try:
                            os.renames(f'{folder}/{name}', f'{folder}/{L2FOLDERPREFIX[sensoridx]}/{name}')
                        except Exception:
                            print(
                                f'Cant move file{name} somethings wrong with the destination or premissions and directory {L2FOLDERPREFIX[sensoridx]} already exists!')

    @staticmethod
    def subdirs(toplevelfolder: str) -> list:
        directory_content = os.listdir(toplevelfolder)
        directories = []
        for item in directory_content:
            if os.path.isdir(f'{toplevelfolder}/{item}'):
                dirr = (f'{toplevelfolder}/{item}')
                directories.append(os.path.normpath(dirr))
        return directories

    @staticmethod
    def match_sensorpaths(folderlist: list) -> dict:
        sensorfolder = []
        sensors = []
        for folder in folderlist:
            for item in L2FOLDERPREFIX:
                if item in folder:
                    sensorfolder.append(folder)
                    sensors.append(item)
        return dict(zip(sensors, sensorfolder))

    @staticmethod
    def unzip_python(feil, dir_extract_to):
        if os.path.isfile(feil):
            with zipfile.ZipFile(feil, 'r') as zipped:
                zipped.extractall(dir_extract_to)

    @staticmethod
    def untragz_python(feil, dir_extract_to):
        if feil.endswith(".tar") and os.path.isfile(feil):
            with tarfile.open(feil, "r:") as tarr:
                tarr.extractall(dir_extract_to)
        if "tar.gz" in feil and os.path.isfile(feil):
            with tarfile.open(feil, "r:gz") as tarr:
                tarr.extractall(dir_extract_to)


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


def psr(pattern, my_path):
    files = glob.glob(my_path + '/**/*' + pattern + "*", recursive=True)
    return files


def psrp(pattern, my_path):
    files = glob.glob(my_path + '/**/*' + pattern, recursive=True)
    return files


def psrp2(pattern1, pattern2, my_path):
    files = glob.glob(my_path + '/**/*' + pattern1 + '*' + pattern2, recursive=True)
    return files


def unzip(path, pattern="zip"):
    fliste = psr(pattern, path)
    for i in enumerate(fliste):
        a = i[1]
        # print(a)
        subprocess.call(["unzip", a, "-d", path])
    return None


# Für L2C_T1 Landsat!
def gunzip(path, pattern="tar"):
    fliste = psr(pattern, path)
    for i in enumerate(fliste):
        a = i[1]
        # print(a)
        dn = a.split('/')[-1]
        dn = dn.split('.')[0]
        os.mkdir(path + '/' + dn + "_srdata_")
        subprocess.call(["tar", "xvf", a, "-C", path + '/' + dn + "_srdata_"])
    return None


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


def atms(path, pattern="L1C", b="60"):
    fliste = psr(pattern, path)
    for i in enumerate(fliste):
        a = i[1]
        subprocess.call(
            ['/home/christian/neo_sen2cor/Sen2Cor-02.08.00-Linux64/bin/L2A_Process', "--resolution", b, a])  # Anpassen!
    return None


def atmsall(path, pattern="L1C"):
    fliste = psr(pattern, path)
    for i in enumerate(fliste):
        a = i[1]
        subprocess.call(['/home/christian/neo_sen2cor/Sen2Cor-02.08.00-Linux64/bin/L2A_Process', a])  # Anpassen!
    return None


def append_to_hdr(feilname, string1, string2):  # Funktioniert um sachen zu appenden!
    z = open(feilname, 'a')
    if '\n' not in string1:
        string1 = string1 + '\n'
    if '\n' not in string2:
        string2 = string2 + '\n'
    z.writelines(string1)
    z.writelines(string2)
    z.close()
    return 0


def locglob(path, pat):
    L = []
    # print(pat)
    # print(path)
    for j in enumerate(pat):
        a = glob.glob(path + "/" + j[1])
        a = a[0]
        L.append(a)
    return L


def stack(path):
    p60 = "60m"
    p20 = "20m"
    p10 = "10m"
    fliste10 = psrp(p10, path)
    fliste20 = psrp(p20, path)
    fliste60 = psrp(p60, path)
    L20 = []
    L10 = []
    L60 = []
    L20.append("wavelength={490, 560, 665, 705, 740, 783, 865, 1610, 2190}" + '\n')
    L20.append("fwhm={65, 35, 30, 15, 15, 20, 20, 90, 180}" + '\n')
    L10.append("wavelength={490, 560, 665, 842}" + '\n')
    L10.append("fwhm={65, 35, 30, 115}" + '\n')
    L60.append("wavelength={443, 490, 560, 665, 705, 740, 783, 865, 945, 1610, 2190}" + '\n')
    L60.append("fwhm={20, 65, 35, 30, 15, 15, 20, 20, 20,90, 180}" + '\n')
    try:
        for i in enumerate(fliste10):
            a = i[1]
            os.chdir(a)
            bnm = a.split("/")[-3]
            bnnm = a + "/" + bnm
            b10 = bnnm + "_10m_stack.vrt"
            pats = ["*B02*.jp2", "*B03*.jp2", "*B04*.jp2", "*B08*.jp2"]
            x = locglob(a, pats)
            print(type(b10))
            print(type(x[3]))
            subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", b10] + x)
            subprocess.call(["gdal_translate", "-of", "ENVI", bnnm + "_10m_stack.vrt", bnnm + "_10m_stack_file"])
            append_to_hdr(bnnm + "_10m_stack_file.hdr", L10[0], L10[1])
    except(Exception):
        print('Fehler!')
    try:
        for i in enumerate(fliste20):
            a = i[1]
            os.chdir(a)
            bnm = a.split("/")[-3]
            bnnm = a + "/" + bnm
            b10 = bnnm + "_20m_stack.vrt"
            pats = ["*B02*.jp2", "*B03*.jp2", "*B04*.jp2", "*B05*.jp2", "*B06*.jp2", "*B07*.jp2", "*B8A*.jp2",
                    "*B11*.jp2", "*B12*.jp2"]
            x = locglob(a, pats)
            # print(x)
            subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", bnnm + "_20m_stack.vrt", b10] + x)
            subprocess.call(["gdal_translate", "-of", "ENVI", bnnm + "_20m_stack.vrt", bnnm + "_20m_stack_file"])
            append_to_hdr(bnnm + "_20m_stack_file.hdr", L20[0], L20[1])
    except(Exception):
        print('Fehler!')
    for i in enumerate(fliste60):
        a = i[1]
        os.chdir(a)
        bnm = a.split("/")[-3]
        bnnm = a + "/" + bnm
        b10 = bnnm + "_60m_stack.vrt"
        pats = ["*B01*.jp2", "*B02*.jp2", "*B03*.jp2", "*B04*.jp2", "*B05*.jp2", "*B06*.jp2", "*B07*.jp2", "*B8A*.jp2",
                "*B09*.jp2", "*B11*.jp2", "*B12*.jp2"]
        x = locglob(a, pats)
        subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", bnnm + "_60m_stack.vrt", b10] + x)
        subprocess.call(["gdal_translate", "-of", "ENVI", bnnm + "_60m_stack.vrt", bnnm + "_60m_stack_file"])
        append_to_hdr(bnnm + "_60m_stack_file.hdr", L60[0], L60[1])
    return None


def deliver_bandnames(csvdata):
    Namen = []
    listus = open(csvdata).readlines()
    for j in enumerate(listus):
        nam = j[1].split(',')
        Namen.append(nam)
    Namen = str(Namen)
    Namen = Namen.replace(']', '}\n')
    Namen = Namen.replace('[', 'band names={')
    return Namen


# Nur in Verbindung mit Untar Skript, sonst ist Alles anzupassen.
def stackls8(path):
    p10 = "_srdata_"
    # p10="SC202002"
    fliste10 = psrp(p10, path)
    L10 = []
    L10.append("wavelength={443, 482, 562, 655, 865, 1610, 2200}" + '\n')
    L10.append("fwhm={20, 65, 75, 50, 40, 100, 200}" + '\n')
    print(fliste10)
    try:
        for i in enumerate(fliste10):
            a = i[1]
            print(a)
            os.chdir(a)
            bnm = a.split("/")[-1]
            bnnm = a + "/" + bnm
            b10 = bnnm + "_sr_stack.vrt"
            pats = ["*SR_B1*", "*SR_B2*", "*SR_B3*", "*SR_B4*", "*SR_B5*", "*SR_B6*", "*SR_B7*"]
            x = locglob(a, pats)
            print(type(b10))
            print(type(x[3]))
            subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", b10] + x)
            subprocess.call(["gdal_translate", "-of", "ENVI", b10, bnnm + "_srenv_stack_file"])
            append_to_hdr(bnnm + "_srenv_stack_file.hdr", L10[0], L10[1])
    except(Exception):
        print('Huch')
    return None


# Nur in Verbindung mit Untar Skript, sonst ist Alles anzupassen.
def stackls4_7(path):
    p10 = "_srdata_"
    # p10="SC202002"
    fliste10 = psrp(p10, path)
    L10 = []
    L10.append("wavelength={479, 560, 661, 834, 1650, 2209}" + '\n')
    L10.append("fwhm={70, 80, 60, 140, 200, 270}" + '\n')
    print(fliste10)
    try:
        for i in enumerate(fliste10):
            a = i[1]
            print(a)
            os.chdir(a)
            bnm = a.split("/")[-1]
            bnnm = a + "/" + bnm
            b10 = bnnm + "_sr_stack.vrt"
            pats = ["*SR_B1*", "*SR_B2*", "*SR_B3*", "*SR_B4*", "*SR_B5*", "*SR_B7*"]
            x = locglob(a, pats)
            print(type(b10))
            print(type(x[3]))
            subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", b10] + x)
            subprocess.call(["gdal_translate", "-of", "ENVI", b10, bnnm + "_srenv_stack_file"])
            append_to_hdr(bnnm + "_srenv_stack_file.hdr", L10[0], L10[1])
    except(Exception):
        print('Huch')
    return None


# Um die Koordinateninfos auf alle Header zu verteilen!!!

def listing(path):  # listing all files in this directory, which ends with .hdr
    return [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".hdr")]


def listing_ftiff(path):  # listing all files in this directory, which ends with _geotiff
    return [os.path.join(path, f) for f in os.listdir(path) if f.endswith("_geotiff")]


def rw_gtif(pfad):
    pff = listing_ftiff(pfad)
    for i in enumerate(pff):
        x = i[1].strip()  # entfernt führenden und nachfolgenden Leerzeichen
        print(x)
        x2 = x + '_env'
        # subprocess.Popen(['python','/usr/bin/gdal_merge.py','-o',basename+'VNIR_SWIR_TOAR_stack','-ps','30','30','-separate','-n','0','-of','ENVI',basename+'_VNIRTOAR',basename+'_SWIRTOAR'])
        pp = subprocess.Popen(['gdal_translate', '-of', 'ENVI', x, x2])
        pp.wait()
    return None


def rewrite_headers(
        inputbildhdr):  # read orginal header and read map info and crs, create a new list, append map info and crs into empty list (L[]), write L into all hdr files
    pf = inputbildhdr.rfind('/')
    pfad = inputbildhdr[:pf]
    liste = listing(pfad)
    inputs = axf.read_hdr_flt(inputbildhdr)
    L = []
    try:
        mapinfo = inputs.map_info
        mapinfo = 'map info={' + mapinfo + '}\n'
        L.append(mapinfo)
    except(AttributeError):
        print('Keine Koordinateninfo vorhanden!!!!')
        return -1
    try:
        crs = inputs.coordinate_system_string
        crs = 'coordinate system string={' + crs + '}\n'
        L.append(crs)
    except(AttributeError):
        print('Keine CRS Info vorhanden, halb so wild!')
    for j in enumerate(liste):
        if j[1] == inputbildhdr:
            continue
        open(j[1], 'a').writelines(L)
    return None


class SensorReflectanceAnalysis:
    def __init__(self, parentfolder):
        self.parentfolder = parentfolder

    def get_all_ref_files(self):
        self.subfolders = FolderFunctions.subdirs(self.parentfolder)
        self.sensordict = FolderFunctions.match_sensorpaths(self.subfolders)
        refdata_files = []
        sensor_keys = []
        for sensor in self.sensordict.keys():
            if sensor.startswith('L'):
                filelist = psrp("srenv_stack_file", self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append(sensor)
            if sensor.startswith('S'):
                filelist = psrp(SENTINEL210MREFLECTANCESUFFIX, self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append("S2_10m")
                filelist = psrp(SENTINEL220MREFLECTANCESUFFIX, self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append("S2_20m")
                filelist = psrp(SENTINEL260MREFLECTANCESUFFIX, self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append("S2_60m")
        self.sensor_ref_filedict = dict(zip(sensor_keys, refdata_files))

    def get_all_ref_files(self):
        self.subfolders = FolderFunctions.subdirs(self.parentfolder)
        self.sensordict = FolderFunctions.match_sensorpaths(self.subfolders)
        refdata_files = []
        sensor_keys = []
        for sensor in self.sensordict.keys():
            if sensor.startswith('L'):
                filelist = psrp("srenv_stack_file", self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append(sensor)
            if sensor.startswith('S'):
                filelist = psrp(SENTINEL210MREFLECTANCESUFFIX, self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append("S2_10m")
                filelist = psrp(SENTINEL220MREFLECTANCESUFFIX, self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append("S2_20m")
                filelist = psrp(SENTINEL260MREFLECTANCESUFFIX, self.sensordict[sensor])
                if len(filelist) > 0:
                    refdata_files.append(filelist)
                    sensor_keys.append("S2_60m")
        self.sensor_ref_filedict = dict(zip(sensor_keys, refdata_files))

    # remember to include the full ratiosuffs of S2 (e.g. 10m_stacked_file_ndvi_, instead of only _ndvi_)!
    def get_all_ref_files_ratio_anaysis(self, param_suffix_list):
        self.subfolders = FolderFunctions.subdirs(self.parentfolder)
        self.sensordict = FolderFunctions.match_sensorpaths(self.subfolders)
        refdata_files = []
        sensor_keys = []
        for ratiosuffix in param_suffix_list:
            for sensor in self.sensordict.keys():
                if sensor.startswith('L'):
                    filelist = psrp(ratiosuffix, self.sensordict[sensor])
                    if len(filelist) > 0:
                        stack_analyze(self.sensordict[sensor], ratiosuffix)
                        refdata_files.append(filelist)
                        sensor_keys.append(f'{sensor}_{ratiosuffix}')
                if sensor.startswith('S'):
                    filelist = psrp(ratiosuffix, self.sensordict[sensor])
                    if len(filelist) > 0 and ("60m" in filelist[0] or "10m" in filelist[0] or "20m" in filelist[0]):
                        stack_analyzeSen2(self.sensordict[sensor], ratiosuffix)
                        refdata_files.append(filelist)
                        sensor_keys.append(f'{sensor}_{ratiosuffix}')
        self.sensor_ref_filedict = dict(zip(sensor_keys, refdata_files))


# Bitte hier weiter batch ratio,batch ratio stacks

def stack_analyzeS5_reduced(path, flag=1):
    p60 = "srenv_stack_file"
    fliste60 = psrp(p60, path)
    L = []
    namliste = []
    for i in enumerate(fliste60):
        a = i[1]
        if '.hdr' in a:
            continue
        if flag == 1:
            ratios_ls5_reduced.ls7_georat_single(a)
            ratios_ls5_reduced.ls7_vegrat1_single(a)
        if flag == 0:
            ratios_ls7u5_dai.ls7_georat_single(a)
            ratios_ls7u5_dai.ls7_vegrat1_single(a)
    return None


def stack_analyzeS7_reduced(path, flag=1):
    p60 = "srenv_stack_file"
    fliste60 = psrp(p60, path)
    L = []
    namliste = []
    for i in enumerate(fliste60):
        a = i[1]
        if '.hdr' in a:
            continue
        if flag == 1:
            ratios_ls7_reduced.ls7_georat_single(a)
            ratios_ls7_reduced.ls7_vegrat1_single(a)
        if flag == 0:
            ratios_ls7u5_dai.ls7_georat_single(a)
            ratios_ls7u5_dai.ls7_vegrat1_single(a)
    return None


def stack_analyzeS8_reduced(path, flag=1):
    p60 = "srenv_stack_file"
    fliste60 = psrp(p60, path)
    L = []
    namliste = []
    for i in enumerate(fliste60):
        a = i[1]
        if '.hdr' in a:
            continue
        if flag == 1:
            ratios_ls8_reduced.ls8_georat_single(a)
            ratios_ls8_reduced.ls8_vegrat1_single(a)
        if flag == 0:
            ratios_ls8_dai.ls8_georat_single(a)
            ratios_ls8_dai.ls8_vegrat1_single(a)
    return None


def stack_analyzeSen2_reduced(path, flag=1, flagm=[20, 10]):
    if 20 in flagm:
        p60 = "_20m_stack_file"
        fliste60 = psrp(p60, path)
        L = []
        namliste = []
        for i in enumerate(fliste60):
            a = i[1]
            if '.hdr' in a:
                continue
            if flag == 1:
                sent2_ratios_reduced.S2_20mgeorat_single(a)
                sent2_ratios_reduced.S2_20m_vegrat1_single(a)
            if flag == 0:
                sent2_ratios_dai.S2_20mgeorat_single(a)
                sent2_ratios_dai.S2_20m_vegrat1_single(a)
    if 10 in flagm:
        p60 = "_10m_stack_file"
        fliste60 = psrp(p60, path)
        L = []
        namliste = []
        for i in enumerate(fliste60):
            a = i[1]
            if '.hdr' in a:
                continue
            if flag == 1:
                sent2_ratios_reduced.S2_10m_vegrat1_single(a)
            if flag == 0:
                sent2_ratios_dai.S210m_georat_single(a)
                sent2_ratios_dai.S2_10m_vegrat1_single(a)
        return None


# Time slice Stacks for single Data Files
# Can be used with LS5 File Structure maybe needs to be adapted to fit other sensors for correct data name
def stack_analyze(path, ratiosuff):
    p60 = ratiosuff
    fliste60 = psrp(p60, path)
    L = []
    namliste = []
    dataliste = []
    for i in enumerate(fliste60):
        a = i[1]
        if '.hdr' in a:
            continue
        namliste.append(a.split('/')[-1][17:26])  # get just the data names
        dataliste.append(a)
    dt = str(namliste).replace("[", "{").replace("]", "}")
    dt = "band names=" + dt + '\n'
    b10 = path + "/_" + ratiosuff + "_stack.vrt"
    pats = dataliste
    subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", b10] + dataliste)
    subprocess.call(["gdal_translate", "-of", "ENVI", b10, path + "/_" + ratiosuff + "_stacked_file"])
    open(path + "/_" + ratiosuff + "_stacked_file.hdr", 'a').write(dt)
    return None


# when using the ratiosuff please also supply the full resolution specifying string (e.g. 10m_stacked_file_ndvi_, instead of only _ndvi_)!
def stack_analyzeSen2(path, ratiosuff):
    p60 = ratiosuff
    fliste60 = psrp(p60, path)
    L = []
    namliste = []
    dataliste = []
    for i in enumerate(fliste60):
        a = i[1]
        if '.hdr' in a:
            continue
        namliste.append(a.split('/')[-1][19:35])  # get just the data names
        dataliste.append(a)
    dt = str(namliste).replace("[", "{").replace("]", "}")
    dt = "band names=" + dt + '\n'
    b10 = path + "/_" + ratiosuff + "_stack.vrt"
    pats = dataliste
    subprocess.call(["gdalbuildvrt", "-r", "cubic", "-separate", b10] + dataliste)
    subprocess.call(["gdal_translate", "-of", "ENVI", b10, path + "/_" + ratiosuff + "_stacked_file"])
    open(path + "/_" + ratiosuff + "_stacked_file.hdr", 'a').write(dt)
    return None


# Still experimental
# Hypertemproal edge detection:
# gradient,abs_grad,absgradsmooth,cannyedge,gradsmooth
def do_hypertemporal_edge(data, geoinfoflag=0):
    cubo = gdalnumeric.LoadFile(data)
    out = heda.heda_archstyle(cubo, 'Hyperion', winsize=3)
    out = numpy.asarray(out)
    auxfun.schreibeBSQ(out, data + '_hteda')
    nn = "band names= {gradient,abs_grad,absgradsmooth,gradsmooth}" + '\n'
    open(data + '_hteda.hdr', 'a').writelines(nn)
    if geoinfoflag == 0:
        try:
            geoinfo(data + '.hdr', data + '_hteda.hdr')
        except(Exception):
            print("no SRS string found")
    return None
