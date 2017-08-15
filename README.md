# DCT-Photometry

Perform photometry on images from the Large Monolithic Imager at the Discovery Channel Telescope.

---

## [`LMI_Photometry.py`](LMI_Photometry.py)

This module contains three functions that produce useful photometry from raw `.FITS` images from LMI at DCT:

### Data_Reduction

*func* `LMI_Photometry.` **`Data_Reduction`** (*directory, filters, targets, save_to=None, dark_exp=1.0, subtract_dark=False*)

Creates and applies a master bias, flat, and dark (optional) frame to science images, saves in directory, and updates the `.FITS` header to make targets Simbad-compatible.

* **`Parameters:`**

   **directory** : str
    
   A directory containing raw `.FITS` images and calibration frames
    
   **filters** : dict
    
   Filters used and corresponding flat exposures
    
   *e.g. {'filter' : flat exposure}*
  
   **targets** : dict
    
   "SCITARG" name in `.FITS` header and corresponding name in Simbad
    
   *e.g. {'FITS target name' : 'Simbad target name'}*
            
   **save_to** : str, optional (default=`None`)
    
   Optional second directory to save calibrated frames to
    
   **dark_exp** : float, optional (default=`1.0`)
    
   Exposure time for dark frames
    
   **subtract_dark** : bool, optional (default=`False`)
    
   Set to `True` in order to subtract dark frame
    
   *Note: LMI has negligible dark current*
    
* **`Returns:`**

   None

### Aperture_Photometry

*func* `LMI_Photometry.` **`Aperture_Photometry`** (*directory, ap_radius, standards, show_figures=False*)

Measures raw electron counts for a target star and utilizes the .FITS header to calculate and save fluxes and instrumental magnitudes. Flags specified stars as standards to be used for standard magnitude transformations.

* **`Parameters:`**
   
   **directory** : str
    
   A directory containing reduced `.FITS` images
    
   **ap_radius** : int
    
   Radius of aperture used for photometry
  
   **standards** : dict
    
   Simbad-compatible name with list of standard star names in the field
    
   *e.g. {'Query Name' : ['Standard Query Name', 'Standard Query Name']}*
            
   **show_figures** : bool, optional (default=`False`)
    
   Display optional figures that are relevant
    
* **`Returns:`**

   None            
            
### Convert_Magnitudes

*func* `LMI_Photometry.` **`Convert_Magnitudes`** (*directory, filters, bin_size=10, show_figures=False*)

Reads magnitudes and airmass values saved in the `.FITS` headers of standard stars, calculates a magnitude transformation for each filter used, then applies the transformation to science images to convert their instrumental magnitudes to standard magnitudes. Saves measurements and uncertainties in a .txt table in `ascii` format.

* **`Parameters:`**

   **directory** : str
    
   A directory containing reduced `.FITS` images instrumental magnitudes appended to the `.FITS` header
    
   **filters** : list
    
   A list of filters used
   
   *e.g. ['filter 1', 'filter 2']*
   
   **bin_size** : int, optional (default=`10`)
   
   Number of epochs target is observed
   
   **show_figures** : bool, optional (default=`False`)
   
   Display optional figures that are relevant
    
* **`Returns:`**
   
   None            

## [`ViewImage.py`](ViewImage.py)

Displays a specified `.FITS` image. Useful for visual inspection for hot pixels, cosmic rays, saturation, etc.
