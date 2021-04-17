# README DampingUI version 1#

The Damping tool, written for Mevislab 3.1.1, reads 4DPCA MRI data, from which velocity and flow curves can be extracted. This document will give an introduction on how to get the application running and how to use it.

### How do I get set up? ###

####Summary of set up####

Download the Damping repository. At least, all files named Colorbar, and Damping are needed. If you have an .exe installer for the Damping UI, you do not need anything besides the installer and a valid product key.

####Dependencies####

Mevislab (tested for v3.1.1)

Python tools:  
sklearn (v0.0)    
matplotlib (v1.4.3)  
scipy (v1.1.0)   
scikit-image(v0.15.0)

Python tools for Mevislab can be installed by adding the 'PythonPip' module to the network. Install the packages by typing their names, without version numbers, in the package name textbox, and subsequently pressing the install button. Note that MevisLab needs to be opened with admin rights for this to work.
PythonPip is not available for MevisLab versions lower than 3.1. There is, to my knowledge, no solution for this.

Fixed version of: CSOMarchingSquares.cpp. Copy this file to ...\MeVisLab3.1VS2017-64\Packages\MeVisLab\Standard\Sources\ML\MLCSO\CSOTools\CSOMarchingSquares.cpp. Hereafter, run MeVisLab Tool Runner-> File -> Add files from directory -> ...\MeVisLab3.1VS2017-64\Packages\MeVisLab\Standard\Sources\ML\MLCSOModules\ (Select Folder)-> MLCSOModules.pro (Run on Selection) ->  recompile the generated MLCSOModules.vcxproj in VisualStudio. 
####Deployment instructions####

Run the network by going to Scripting -> Start network script (Ctrl+R)

### How to use the Damping tool?###

####Import PCA Images####
Dicoms can be loaded in the Import screen. Drag a 4DPCA image into the GUI. The GUI will automatically find all other velocity encoded directions in the same folder, and load those as well. It will flag it as a 4DPCA image if the ImageType in the header contains the word 'PCA'. A subfolder is created 'MLoutput', in which the magnitude, modulus and phase images are stored separately.  
In addition, when loading the 4DPCA images, the tool will search the current image folder for a T1W TFE image, and will load this into the top right panel of the IMPORT PCA IMAGES tab. It is recommended to use the T1W image to draw centerlines on, since the FOV of the PCA images is often limited.
When the loading bar is at 100%, the Dicoms will be split and loaded. Note, that if the Dicoms are loaded a second time, the computational time is significantly less, since it will load directly from the previously created MLoutput folder.

*Options* 
 
- **Perform image corrections (Button)**: will correct an offset/gradient in the static tissue and will perform phase unwrapping. Pressing this button will overwrite the PPCA images in the MLoutput subfolder. Original images are untouched. 
 
- **Using  corrected scans (Checkbox)**: Is checked if the corrected images are used.  

- **Open registration (Button)**: Opens panel for registration of 4DPCA to optional T1W image. 

- **Save registration results (Button)**: Saves the registration results when registration is performed in the panel opened by the OPEN REGISTRATION button.

- **Inspect registration results (Button)**: Opens panel to check if the 4D PCA images are aligned with the T1W image.

- **Manual load image (Button)**: Option to load MFFE, PPCA and MPCA images manually.

- **Ignore existing ML output (Checkbox)**: If checked, loading dicoms will always resplit the dicoms and save them in the subfolder. If unchecked, the GUI will not resplit if the MLoutput folder already exists.

- **Override: use T1W(Checkbox)**: use T1W. If a T1W image is loaded in, but you would like to draw your vessel center lines based on the MFFE image, deselect this checkbox. If checked, the T1W image is shown in the following tabs.

- **Invert Phase (RL/FH/AP) (Checkboxes)**: Use this if the direction of  the phases is incorrectly shown. Checking this box will flip the phase for the given phase image. 

- **Help! (Button)**: Opens this readme file in an external pdf reader.

*Possible problems*  

- *Only the magnitude image is shown*: Look at the Debug output, if it gives a message 'non orthogonal/no velocity encoding directions found in header', your header information is incorrect, and the images will not load. 
 
- *Nothing is loaded*: Have you given a PCA image as input? It will not load anything if the image is not a PCA image. If your image is a PCA image, does the ImageType tag in the header contain the word 'PCA'?

####Center line tracking####
In this tab, the vessels can be selected from which the flow curves need to be calculated. The 3DViewer will show the vessels present in the loaded magnitude image. Image can be rotated and translated using mouse buttons.  

To select a certain vessel, hold the Alt key pressed. This will allow to select a point in the centerline of a vessel of interest. Select a second marker further downstream. If a path is detected between both points, a center line is automatically shown in green. Pressing the 'Compute Pulsatility' will compute the pulsatility and vessel contours for along the center line.  

Misclicked markers can be deleted by pressing the 'Delete all markers' button. Controls are also explained in the Controls panel.

*Options*  

- **Show vessel segmentation (Checkbox)**: Shows the raw segmentation of the vessel, after COMPUTE PULSATILITY has been pressed. This option is purely for visuals.

- **Show selection skeleton (Checkbox)**: will show the centerline surface generated when the Alt key is pressed while in the 3DViewer.

- **Skull strip (Checkbox)**: Basic skull stripping to make it easier to see the intracranial arteries.

- **Create manual vessel mask (Button)**: Opens window to draw and additional mask on a loaded T1W image. This option is only necessary if the intensity of the target vessel is too low, and the center line is not detected. The saved manualmask, can be loaded in using the load option on the same GUI tab.

- **LUT**: By changing the horizontal position of the second point (red), the rendering of the vessel changes. This does not change the detection of the vessels.

- **Experimental options (Buttons)**: Experimental. Will open a window with extra options:  
  
  * **Set non-vessel intensity threshold**: A higher value will delete more of the potential vessels, as it thresholds the output of the vesselness filtered image.

  * **Set FOV (mm)**: Will change the threshold used in the DTF method, to filter out things that are not vessel. If set too high, some smaller vessels, e.g. MCA, will be unselectable.

  * **Set DTF threshold (Slider)**: Will change the threshold used in the DTF method, to filter out things that are not vessel. If set too high, some smaller vessels, e.g. MCA, will be unselectable.
     
  * **FWHM calculation windows width**: Sets the number of MPR slices around the current slice, which are used to calculate the FWHM value to set the iso contour at. Only the current MPR slice is used to calculate the iso value, if the value of this option is set to 0. If this value is set larger than the number of slices in the MPR, a single isovalue is calculated, that is used for all MPR slices.

  * **Centerline smoothing (Slider)**: Will set the smoothing for the DTF method centerline.
  
  * **Use variable FWHM (Checkbox)**: Opens number edit if checked. If checked, the FWHM iso value is calculated per time point. The number edit can be used to set the level of the FWHM isovalue. A value of 0.5 sets the isovalue at the FWHM value. A value of 1 and 0would set the isovalue at the median foreground/background intensity, respectively. If unchecked the FWHM isovalue is determined on the first dynamic and used for all other time points. 
 

- **Load results (Browse)**: Can be used to load in already saved results using this GUI.

- **Load manual mask (Browse)**: Open a manual drawn vessel mask. This mask is added to the vessel segmentation performed by mevislab. 

####View Pulsatility####
In this tab the same 3Dviewer is present. This will show the measurement locations, given by the spheres, along the centerline, combined with a white outline, giving the orientation of the MPR plane. 

In the MPR panel, the magnitude and velocity maps are given in their respective tabs. The colorbar in the magnitude tab gives the velocity in cm/s. Scrolling through the MPR will make you go through all the slices. The left and right arrow buttons will allow for going through all dynamics.

In green a contour is given per slice. The velocity values in this contour are used for further computations. Holding the Ctrl button, and clicking using the left mouse button and scrolling will allow you to set your own isocontour, if the automatic method failed. 
Iso contour redrawing can only be done in the Magnitude tab. When scrolling through the slices, the 3DViewer will show the current slice. In the CURRENT SLICE tab, shown below the MPR panel, the velocity curve is given for the current slice. The ALL SLICES tab will show the max, mean and minimum velocity over all slices.

*Options*

- **Overlay velocity map (Checkbox)**: Will overlay the velocity map in the magnitude MPR image.

- **Sphere size as pulsatility (Checkbox)**:  Will scale the spheres in the 3Dviewer, with the measured pulsatility values.

- **Show pulsatility values (Checkbox)**: Will show the pulsatility values at each sphere.

- **Show MPR orientation (Checkbox)**: Will show the orientation of the MPR slice in both the MPR panel as well as the 3Dviewer.

- **Show vessel segmentation (Checkbox)**: Shows the raw segmentation of the vessel. This option is purely for visuals.

- **Anatomical location current slice (Radio buttons)**: Let's you select if you have performed the analysis on the left or right artery. In addition, you can select per slice to which anatomical location it belongs to (C1..C7). E.g. if for slice 1, C2 is selected, and for slice 15, C3 is selected, the AUTOFILL button will set the anatomical location of slices 1-14 to C1, and 15-end to C2.It is import to always set the anatomical location of the first slice correctly, as the autofill algorithm copies the location in the positive slice direction. If this is not done, the save box in the next tab will warn you that you have not selected an anatomical location for all slices. These variables are purely for analysis purposes outside the scope of this GUI.

- **Maximum, mean (Radio buttons)**: When maximum is used, the maximum velocity of the ROI is used in the Single and All tabs, if mean is used, the average velocity of the ROI is used in the Single and ALl tabs.

####Save####
Will let you save the results to a user defined location. The saved file is in a .npy and .csv format, which can be read by numpy or matlab / excel. The array contains:
- Centerline positions: centerline positions in world coordinates for all MPR slices.

- Centerline vectors: Orientation of all MPR slices.

- Leftright: if the analysis was performed on the left or right artery.

- Anatomical location: The C.. location per slice as integers, where 1=C1 and 7=C7.

- Velocity profiles: velocity profiles in cm/s for all slices (mean or max velocity depending on the used arithmetic)

- Arithmetic: If mean or max velocity was used for pulsatility calculation/ velocity profiles.

- Pulsatility: Pulsatility values per slice

- Path length: Length of the center line in mm.

- Area: Contour area per slice and timepoint in mm2.

- Slice world matrix: World matrix per slice.

- MPR world matrix: World matrix of the entire MPR stack.

- Variable FWHM: If a variable FWHM was used

- FWHM value: corresponding FWHM value. 0=median bg value, 0.5=FWHM, 1=median fg value.

### Who do I talk to? ###

You can contact me at k.m.vanhespen@umcutrecht.nl or by creating an issue on Github.
