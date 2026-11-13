# Data
Description of the data used in this project.

The cells were imaged on an inverted motorized stand Zeiss AxioObserver7 equipped with CSU-W1 Confocal Scanner Unit (Yokogawa) with Dual T2, a 50um pinhole disk unit, 2 sCMOS cameras (Prime 95B, photometrics),a Visitron VS-Homogenizer ,an ASI MS2000 X,Y, ZPiezo drive (300um travel range) and controlled with VisiView (version 6.0.0.35) imaging software (Visitron Systems GmbH). Illumination was achieved with 561 Cobolt Jive (Cobolt, 200mW). For FRAP, illumination was done using a 473 nm FRAP laser.
A plan Apochromat 40X/1.3 oil objective was used, resulting in a pixel size of 0.275 μm in x,y. A 405/488/561 dichroic was used. The emission filter in the scan head was a 575 nm long pass.  Imaging was done in a humidified incubation chamber at 37° C supplied with 8% CO2. 
A ROI was defined on the field of view and few cells (~10) were selected for FRAP. A rectangular ROI was defined for each nuclei encompassing around 80% of the nuclei area at the focal plane. 2D timelapse images were acquired every 30 s for 250 frames, further referred to as movies. Images were taken at 30% max laser power (~ 3.528mW) to limit bleaching. The focus was kept using a Zeiss definiteFocus 2 system between every frame. After 4 frames, the FRAP sequence was initiated on the selected nuclei at 100% laser power. For FRAP, a laser of 22 um was focused on the focal plane and was shined on each pixel for 300ms. The focusing of the FRAP laser to the imaging focal plane and calibration of the x,y direction of the laser was done before every acquisition using a dedicated well not used for the imaging. The calibration was done by moving the beam to specific locations on the field of view and adjusting its focus to obtain the smallest beam at the imaging focal plane. 
The cells were imaged on a µ-Slide 8 Well glass bottom (Ibidi, 80807). Briefly, a day before imaging ~200.000 cells were split in Gibco™ FluoroBrite™ DMEM (Gibco™) onto the dish. On the day of the imaging 100nM of Halo-JF-549 ligand (Janelia Fluor® NHS ester Tocris lab, 6147) was added for 30 minutes to the imaging medium. Then cells were washed three times with 1X PBS followed by addition of fresh imaging media.

## Data structure

The raw data is structured as a .ome.tf2 file with axis (t,y,x). The actual shape of the image analyzed can be found in the .log file of each run.

## Image analysis

For the complete analysis steps see the processing_steps page.



