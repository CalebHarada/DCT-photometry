# DCT-Photometry

Perform photometry on images from the Large Monolithic Imager at the Discovery Channel Telescope.

---

## [`LMI_Photometry.py`](LMI_Photometry.py)

This module contains three functions that produce useful photometry from raw .FITS images from LMI at DCT:

#### Data_Reduction

*func* `LMI_Photometry.` **`Data_Reduction`** (*directory, filters, targets, save_to=None, dark_exp=1.0, subtract_dark=False*)

Creates and applies a master bias, flat, and dark (optional) frame to science images, and updates the .FITS header to make targets Simbad-compatible.

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
    
    **dark_exp** : float (default=`1.0`)
    
    Exposure time for dark frames
    
    **subtract_dark** : bool (default=`False`)
    
    Set to `True` in order to subtract dark frame
    
    *Note: LMI has negligible dark current*
            
            
            
            
`Aperture_Photometry` 
Measures raw electron counts for a target star and utilizes the .FITS header to calculate and save fluxes and instrumental magnitudes. It flags specified stars as standards to be used for standard magnitude transformations.

`Convert_Magnitudes` 
Reads the magnitudes and airmass values saved in the .FITS headers of standard stars, calculates a magnitude transformation for each filter used, then applies the transformation to science images to convert their instrumental magnitudes to standard magnitudes. It saves measurements and uncertainties in a .txt table in 'ascii' format.


## [`ViewImage.py`](ViewImage.py)

Displays a specified `.FITS` image. Useful for visual inspection for hot pixels, cosmic rays, saturation, etc.
