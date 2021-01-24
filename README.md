# Epidermal Z-stack Peeler: EZ-Peeler

First published use: Zoulias N, Brown J, Rowe J, Casson SA (2020) HY5 is not integral to light mediated stomatal development in Arabidopsis. PLOS ONE 15(1): e0222480. https://doi.org/10.1371/journal.pone.0222480

Intially started whilst working in Stuart Casson's lab (Sheffield University), I have continued to work on this as a learning exercised in Alexander Jones' lab (SLCU Cambrige). EZ-Peeler is a top down surface segmentation tool that can be used to extract a non flat surface from a confocal stack. Data is exported as a segmented surface, a surface projection, an areamap and a heightmap. The extraction of a surface layer allows users to apply their chosen 2D image segmentation tools and an areamap is produced to account for area differences from the flattening the image. Whilst the code isn't very clean (writing EZ-Peeler has been a learning exercise in image analysis) it is functional and useful.



**Introduction**
Understanding gas exchange in plants requires an in-depth knowledge of stomatal development and dynamics. Investigating stomatal development requires mechanistic understanding of how the leaf epidermis is patterned. Clearing leaves and imaging under a DIC microscope or taking epidermal impressions are often used to quantify this development.
This method is cheap and effective but doesn’t allow quantification of development or its TF regulators to be tracked in situ. Confocal microscopes allow high resolution imaging of the epidermal cells, with optical sectioning, enabling us to see cell patterning, morphology and protein localisation.  Unfortunately, the epidermis is often contoured (Figure 1), so thin optical sections can fail to capture the whole surface in a single optical plane (Figure 1) and a lot of expensive microscopy time may be spent searching for flat sections. Taking a series of these optical sections at different Z positions (Z-stacks) is the most common solution, but projecting, interpreting and quantifying these three dimensional datasets in an insightful way is difficult. Projecting z-stacks in two dimensions leads to inclusion of out of focus planes, the noise from which can make interpretation difficult (Figure 1) and often includes mesophyll which contain large autofluorescent chloroplasts (Figure 1).


![Figure 1: contoured surfaces present a problem](https://github.com/JimageJ/EZ-Peeler/blob/master/EZPF1.png)
*Figure 1	Optical sectioning presents a problem for contoured surfaces
Top: An XZ-slice from an Arabidopsis abaxial cotyledon surface. Purple: Chlorophyll autofluorescence, grey: Propidium iodide, Green: GFP-RGA
Bottom: Three XY sum projections of the same Arabidopsis abaxial cotyledon surface and their corresponding position in the Z stack. Projection A shows mostly epidermis, but does not include large sections below the stack depth. Projection B shows the missing epidermal sections, but a lot of mesophyll bleedthrough. Projection C shows the entire epidermis, with considerable amounts of mesophyll bleedthrough.*


**Software development**
The FIJI distribution of ImageJ was chosen as a development platform as it is already widely used in the academic community and its capabilities can be expanded by coding plugins in a number of languages (e.g. Java, Python, Groovy) or by recording macros, allowing complex tasks to be automated (Rueden et al., 2017; Schindelin et al., 2012). As well as being free, open source and cross platform, the vast collection of ImageJ analysis plugins allow new plugins to be integrated into existing workflows. EZ Peeler is written in Jython and available as a .py source code. EZ Peeler is dependent on CLIJ, which is used to process images directly on the GPU, allowing repetitive image processing tasks to be highly parallelised, improving processing speed.  In the following sections, the image processing workflow is outlined, as summarised in Figure 2.

![Figure 2: plugin workflow](https://github.com/JimageJ/EZ-Peeler/blob/master/EZPF2.png)
*Figure 2.	EZ Peeler plugin workflow
A) Pre-processing: User selects a timepoint and channel to be used for the initial segmentation, Apply optional XY(Z) smoothing filters apply optional XZ edge detection filter B) Segmentation / error correction: Raster through Z, X then Y, searching for the first voxel at every XY position that is brighter than a user defined threshold value. Replace missing data/segmentation errors with linear interpolations. Make a binary mask based of defined thickness on these points. Generate a heightmap of the surface. C) Apply segmentation to original image frame: Apply mask to original image.
Sum Z projection to create surface projection.*

**Pre-processing**
The channel to be used for segmentation is duplicated, allowing it to be filtered and manipulated for segmentation whilst leaving the original image unchanged. A pre-processing smoothing filter may be applied in XY or XYZ (2D Gaussian, 2D median or 3D median) to remove noise and improve segmentation. As the signal intensity of the surface may not be uniform, edge detection filters can be applied in the XZ plane (3X3 1D Sobel, 3X3 2D Sobel, 3 X 3 Laplace filter, 1 X 3 Gradient, 1 X 7 Gradient, 3 X 7 Gradient, Figure S1). Most of these filters convert the image into a differential of signal intensity in the Z direction, making segmentation less brightness-dependent, but also increasing noise. We recommended using a XY(Z) smoothing filter that reduces specular noise (e.g. 2D median) when using edge detection filters, to reduce segmentation errors.

**Segmentation**
To detect the epidermis at a particular X,Y coordinate, a column of voxels is scanned in the Z dimension. The first voxel from the top with greater than a threshold intensity is recorded as the ‘top’ of the epidermis (Figure 2 B). The columns are rastered through in X then Y directions. When no pixel in the column is above the specified threshold, linear interpolation along the X-axis is used to automatically fill in the missing data. The user can specify if they would like to apply smoothing interpolation to the segmentation results. This consists of applying a 2D 3-point moving average to the original detected points and using bilinear interpolation to space vertices at user defined interval. 

The results of the segmentation are presented as multipoint ROIs on the XZ slices of the processed stack, allowing fast visual inspection for segmentation errors. Temporary XY heightmap and divergence maps are also generated for segmentation error identification and removal. 




**Error checking, semiautomated correction and trichome removal**

Two forms of error can be corrected in a semiautomated way. Particulate matter floating above the epidermis may create floating spots in the surface, or if the top surface of the epidermis is not found and so the next threshold value is a deeper aberrant object, resulting in holes (Figure 3A, B, C). In both cases, the aberrant surface height will differ significantly from the surrounding surface heights. Generating a heightmap (Figure 3C, D) and applying a large sigma Gaussian blur, we get a smooth approximation of the surface. By subtracting the observed heightmap from the smoothed heightmap, a surface ‘divergence map’ is generated, where large positive or negative pixel values represent sudden changes in surface height, or hotspots (Figure 3E, F). A divergence threshold is set by the user to select these hotspots and successive gaussian blurs  (without blurring non selected pixels) are used to approximate the missing data (Figure 3 G, H and I). These hotspots normally represent errors in segmentation, but the same methods can be used to remove protruding leaf features such as large trichomes with moderate success.

![Figure 3: Error checking](https://github.com/JimageJ/EZ-Peeler/blob/master/EZPF3.png)
*Figure 3	Segmentation choice and error correction
A) XZ slice illustrating a segmentation error, as the epidermal surface is not recognised B and C) Segmentation errors on heightmaps, D and E) Errors on divergence maps, F and G) Heightmaps affter error removal I) Corrected segmentation J and K) Woodgrain effect removal by Gaussian filtering of binary mask.*
    
**Applying segmentation choices to final image**
The user can then define how thick a Z peel they would like (in voxels) and large an offset from the epidermal surface the peel should be (Figure 3J). Two methods are available for selecting the peel, linear Z offset and 3D erosion (Figure 3K, L). Linear Z offset is computationally much faster and provides good segmentation of surfaces such as leaf epidermis which are approximately flat XY planes. This is a similar approach to that of SurfCut (Erguvan et al., 2019). It also outputs a constant Z thickness, meaning that the surface projection is less variable in background intensity. This is especially useful when projecting/segmenting flat surfaces. 3D erosion based surface selection is broadly equivalent to that of MorphographX (de Reuille et al., 2015) and will give more accurate results for contoured surfaces.
A binary stack is generated of the epidermis, which is applied to the original image to generate a segmented surface stack (Figure 2C). Image scale properties from the original image are applied to the segmented surface stack. A Z-sum projection of this segmented surface stack is generated.
A ‘woodgrain’ artefact in final images was apparent and caused by a lack of antialiasing in the binary mask. This can be corrected by using an optional XY Gaussian filter of the binary mask (Figure 3M, N). 	

To account for curvature an ‘areamap’ is generated. Dividing the heightmap into triangular polygons, with each pixel a vertex, basic trigonometry is used to calculate the length of edges, then these lengths are used to calculate the surface area. The areamap uses the heightmap to calculate the surface area for each given pixel of the 2D surface projection. Steeper parts of the surface projection will have a surface area that is correspondingly larger.


**Time series data**
If the image provided is a time series, the user then has the option to apply the same segmentation settings to the whole time series to create a segmented time series (Figure 4). The user has the option to use an otsu threshold determined for each frame, to overcome the effects of bleaching and sample movement.

![Figure 4.	EZ Peeler timeseries workflow](https://github.com/JimageJ/EZ-Peeler/blob/master/EZPF6.png)

*Figure 4.	EZ Peeler timeseries workflow
Steps A) Pre-processing, B) Segmentation / error correction and C) Apply segmentation to original image frame are applied to a single timeframe with user input, as in Figure 2, generating images that can be checked by the user. If the user is happy with the segmentation, these segmentation setting are then used to segment every other timeframe in the image in a loop of steps D) Pre-processing, E) Segmentation / error correction and F) Apply segmentation to original time frame. The segmented time series is rendered as a surface projection timeseries and a series of heightmaps.*

**Results**

![Figure 5: EZ-Peeler surface segmentation of an Arabidopsis cotyledion epidermis](https://github.com/JimageJ/EZ-Peeler/blob/master/EZPF4.png)
*Figure 5.	Comparison of Arabidopsis proRGA::GFP-RGA EZ Peeler surface projection with corresponding unsegmented projections
    A) Sum Z projection of an unsegmented stack B) Sum Z surface projection of EZ Peeler processed stack. C) An XY slice from the unsegmented projection containing mesophyll chloroplasts (purple). D) Corresponding slice from EZ Peeler processed stack. E) XZ slice from the unprocessed stack. F)  XZ slice from EZ Peeled Z stack. G) Oblique angle projection of the unprocessed stack H) Oblique angle projection of the EZ Peeler processed stack. I) Heightmap generated by EZ Peeler J) Areamap generated by EZ Peeler*

3D and surface projections of RGApro::GFP-RGA cotyledons show successful segmentation of the epidermal cells, and multichannel projection Figure(X).  excluding the strongly fluorescent stained surface, the mesophyll cells containing autofluorescent chloroplasts and .


**Software testing**
EZ Peeler has been tested successfully on Arabidopsis thaliana cotyledon and Oryza sativa leaf Z stacks (Figure 5). Time series data has been tested with Arabidopsis thaliana true leaves and Nicotiana Benthamiana leaves (Figure 6). The majority of testing was performed on a 2013 Dell laptop (Windows Intel i5-3337U Dual Core 1.8GHz, 16GB RAM, Intel HD Graphics 4000/AMD Radeon HD 8730M (switchable graphics)), or a 2014 Gigabyte laptop (Ubuntu Intel i7-4710Q 2.5GHz Quad core, 16GB RAM, Nvidia GTX 860M 4gb) on which the software runs well. Considerable speed increases are present on more modern computing hardware. Dozens of Arabidopsis cotyledon z-stacks have been tested.

**Limitations**
EZ-Peeler is designed to segment and project surfaces in a similar orientation to the XY imaging plane. The greater the Z variation of the surface, the less accurate a 2D XY surface projection can be. The limitations of this ‘top down’ surface projection were explored well previously (Erguvan et al., 2019) and other options exist if 2.5D or true 3D segmentation are required (de Reuille et al., 2015).
As well as having an undulating surface, the epidermis is not of uniform thickness and cells of a large XY area also tend to be thicker in Z (Figure 2). Currently EZ Peeler only segments the top surface of the epidermal cells and creates a surface projection of uniform thickness meaning that if the user sets the peel thickness as too thin, important cellular details may be missing from the segmentation.  The basal epidermal periclinal wall of the epidermal cells is not always very bright in confocal images, so we reasoned it would be difficult to segment.


**Requirements and installation**
The FIJI (https://fiji.sc/) distribution of ImageJ is required to run EZ-Peeler. Installation can be achieved by placing the .py file in your Plugins folder. EZ Peeler is dependent on Clij, which is installed via an ImageJ update site.
The plugin has been extensively tested on Windows 7, 10 and Ubuntu 19.10, Fiji has been developed to be cross platform and EZ Peeler only uses plugins included in Fiji as default apart from CLIJ. Images must be converted from RGB before processing.

**Acknowledgements**
This work was funded by the Biotechnology and Biological Sciences Research Council (BBSRC) grants to SC (BB/N002393/1) and AMJ. Thanks to Albert Cardona and Robert Haase for their excellent ImageJ tutorials, resources and code examples.
It works well on Windows 7 and 10 and linux computers (64 bit is preferred as 32- bit processors limit image size in ImageJ). 
