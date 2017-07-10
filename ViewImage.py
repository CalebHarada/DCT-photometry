import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


# Choose .FITS file
filename = '...'
hdu = fits.open(filename)
data = hdu[0].data

# display histogram
plt.figure()
plt.hist(np.asarray(data).flatten(), 500)
plt.yscale('Log')

# display image
plt.figure()
plt.imshow(data, cmap='gray', vmin=np.median(data), vmax=np.median(data) + 500)
plt.colorbar()
plt.title(file)

plt.show()
