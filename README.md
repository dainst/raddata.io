# Scripts for Simple Ratio Transforms 

These python files should be used to calculate ratio transformations from NASA's Landast 5/7 and 8 fleet of instruments. The ratios are according to the Index Database Project (1). The prupose of the scripts is the scientific exploration of remote sensing data in the scientific field of archeoology. 


## Script Layout and Function

### auxfun.py
Offers the necessary file I/O components for the ratio scripts. It Includes the following functions/methods:

- schreibeBSQ
-> For writing 3d spectrosocpic image data cubes in float.
- schreibeBSQ_int 
-> For writing 3d spectrosocpic image data cubes in 8bit unsigned integer.
- schreibeBSQsingle
-> For writing 2d image data in float.

#### read_hdr and read_hdr_flt objects
-Read a header file associated with a data cube and either output the header as object, or as object with attached arrays for parameters such as wavelength and fwhm.

### batch_proc
Offers batch processing capabilites for the ratios_* scripts, which themselves are designed to run for a specific sensor. The usage of the batch_proc script is recommended as a necessary prerequisite for a successful batch processing. I.e. it is necessary to uncompress the data and put it into an according subdirectory (for each image acquisition) with the correct dataset name. In addition data needs to be stacked to a 3d image cube for diplay in any GIS or for further ratio processing all of this is done by the following functions:
- unzip and gunzip (For Landsat Data and Sentinel-2)
- atms and atms_all Sen2Cor Processing for Sentinel-2 using ESAs Sen2cor Application (2) on a local machine. 
- stack (For stacking of Sentinel-2, all available resolutions)

### ls8_ratios

Calculates the ratios from a single stacked LS8 dataset. The calculation process is split into two functions:
- ls8_georat_single calculates the soil/geologic ratios of a single stacked landsat8 file
- ls8_vegrat1_single calculates the vegetation ratios of a single stacked landsat8 file

### ratios_ls7u5

Calculates the ratios from a single stacked LS4,5 and 7 dataset. The calculation process is split into two functions:
- ls7_georat_single calculates the soil/geologic ratios of a single stacked landsat4,5,7 file
- ls7_vegrat1_single calculates the vegetation ratios of a single stacked landsat4,5,7 file

### prisma

Solely offers file I/O and transformation capabilities for ASIs hyperspectral Sensor PRISMA (3). The most interesting point is the provision of HSV Pan-Sharpening via HSV using the 5 m Panchromatic Band of PRISMA, which improves the spatial resolution for visual interpretation quite dramatically. Obviously this file still needs further improovment and testing of the methods. However, it provides data handling capability to a state of the art hyperspectral system with interesting capabilities for archeological research. 

## References

- 1.) Henrich, V., Krauss, G., Götze, C., Sandow, C. (2012): IDB - www.indexdatabase.de, Entwicklung einer Datenbank für Fernerkundungsindizes. AK Fernerkundung, Bochum, 4.-5. 10. 2012.
- (2) MuellerWilm U., Devignot O., Pessiot L., S2 MPC, Sen2Cor Configuration and User Manual Ref. S2-PDGS-MPC-L2A-SUM-V2.5.5 (http://step.esa.int/thirdparties/sen2cor/2.5.5/docs/S2-PDGS-MPC-L2A-SUM-V2.5.5_V2.pdf)
- (3) Loizzo, R.; Guarini, R.; Longo, F.; Scopa, T.; Formaro, R.; Facchinetti, C.; Varacalli, G. PRISMA: The Italian Hyperspectral Mission. In Proceedings of the International Geoscience and Remote Sensing Symposium on Observing, Understanding and Forecasting the Dynamics of our Planet (IGARSS), Valencia, Spain, 22–27 July 2018.