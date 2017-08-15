# DCT-Photometry

Perform photometry on images from the Large Monolithic Imager at the Discovery Channel Telescope.

[`LMI_Photometry`](LMI_Photometry.py)
---

This module contains three functions that produce useful photometry from raw .FITS images from LMI at DCT:

`Data_Reduction` 
Creates and applies a master bias, flat, and dark (optional) frame to science images, and updates the .FITS header to make targets Simbad-compatible.

    param directory: str
            A directory containing raw .FITS images and calibration frames
    param filters: dict
            Filters used and corresponding flat exposures
            {'filter' : flat exposure}
    param targets: dict
            "SCITARG" name in .FITS header and corresponding name in Simbad
            {'FITS target name' : 'Simbad target name'}
    param save_to: str, opt
            Optional second directory to save calibrated frames to
    param dark_exp: float
            Exposure time for dark frames, default is 1.0 sec
    param subtract_dark: bool
            Set to True to subtract dark frame, default is False
            Note: LMI has negligible dark current
    return: None

`Aperture_Photometry` 
Measures raw electron counts for a target star and utilizes the .FITS header to calculate and save fluxes and instrumental magnitudes. It flags specified stars as standards to be used for standard magnitude transformations.

`Convert_Magnitudes` 
Reads the magnitudes and airmass values saved in the .FITS headers of standard stars, calculates a magnitude transformation for each filter used, then applies the transformation to science images to convert their instrumental magnitudes to standard magnitudes. It saves measurements and uncertainties in a .txt table in 'ascii' format.


[`ViewImage`](ViewImage.py)
---

Displays a specified `.FITS` image. Useful for visual inspection for hot pixels, cosmic rays, saturation, etc.
