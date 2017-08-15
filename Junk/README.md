# Junk
This folder contains old pieces of code that were prototypes of [`LMI_Photometry.py`](https://github.com/CalebHarada/DCT-photometry/blob/master/LMI_Photometry.py). Much of it is probably buggy and not very useful.


[`master_bias.py`](master_bias.py)
---

Create a master bias frame for `data_reduction.py`.

[`master_dark.py`](master_dark.py)
---

Create a master dark frame for `data_reduction.py`.

[`master_flat.py`](master_flat.py)
---

Create a master flat frame for `data_reduction.py`.

[`data_reduction.py`](data_reduction.py)
---

Apply bias subtraction and flat field correction to raw `.fits` files.


[`aperture_photometry`](aperture_photometry.py)
---

Do aperture photometry on reduced `.fits` images.

[`PSF_Photometry`](PSF_Photometry.py)
---

Create PSF models from `.fits` images and fit to stars to calculate flux.

[`PSFModelExample.py`](PSFModelExample.py)
---

An example of fitting 2D Gaussian models to contrived data.
