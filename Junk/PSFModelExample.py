from photutils.datasets import make_random_gaussians
from photutils.datasets import make_noise_image
from photutils.datasets import make_gaussian_sources
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
from photutils import find_peaks, DAOStarFinder
from photutils.background import MADStdBackgroundRMS, MMMBackground
from photutils import CircularAperture
from astropy.table import Table
from mpl_toolkits.mplot3d import Axes3D

#=======================================================================================================================

# generate fake data

num_sources = 50
min_flux = 100
max_flux = 500000
min_xmean = 5
max_xmean = 395
sigma_psf = 1.72
starlist = make_random_gaussians(num_sources, [min_flux, max_flux],
                                 [min_xmean, max_xmean],
                                 [min_xmean, max_xmean],
                                 [sigma_psf, sigma_psf],
                                 [sigma_psf, sigma_psf],
                                 random_state=123)

shape = (400, 400)
image = (make_gaussian_sources(shape, starlist) +
         make_noise_image(shape, type='poisson', mean=50., random_state=123) +
         make_noise_image(shape, type='gaussian', mean=0, stddev=2.0, random_state=123))

np.random.seed(0)
y, x = np.mgrid[:400, :400]
g1 = models.Gaussian2D(15125, 201, 200, sigma_psf, sigma_psf)
g2 = models.Gaussian2D(9485, 199, 201, sigma_psf, sigma_psf)
g3 = models.Gaussian2D(2465, 192, 195, sigma_psf, sigma_psf)
z = g1(x ,y) + g2(x ,y) + g3(x, y)

image = image + z

plt.imshow(image, origin='lower', interpolation='nearest')
plt.title('data')
plt.colorbar()
plt.show()

#=======================================================================================================================

# create and fit models to data

bkg = MMMBackground()(image)
std = MADStdBackgroundRMS()(image)
threshold = bkg + 4*std

daofind = DAOStarFinder(threshold=threshold, fwhm=2)
peak_tbl = daofind(image)
circles = CircularAperture((peak_tbl['xcentroid'], peak_tbl['ycentroid']), r=8)

plt.figure(figsize=(8,6))
plt.imshow(image, origin='lower', interpolation='nearest', cmap='viridis', vmax=np.median(image) + 300)
plt.colorbar()
plt.title('Data')
circles.plot()
plt.show()

PSF_model = models.Gaussian2D + models.Const2D
fitter = fitting.LevMarLSQFitter()
def sigmas(frame, x_mean, y_mean, peak, sky):
    size = 25
    if size < y_mean < image.shape[0]-size and size < x_mean < image.shape[1]-size:
        zoom = frame[y_mean - size:y_mean + size, x_mean - size:x_mean + size]
        fit_init = PSF_model(amplitude_0=peak, x_mean_0=zoom.shape[1]/2, y_mean_0=zoom.shape[0]/2,
                             amplitude_1=sky,
                             fixed={'amplitude_1': True},
                             bounds={'x_mean_0': (zoom.shape[1]/2-5, zoom.shape[1]/2+5),
                                     'y_mean_0': (zoom.shape[1]/2-5, zoom.shape[1]/2+5)}
                             )
        y, x = np.mgrid[:zoom.shape[0], :zoom.shape[1]]
        fit = fitter(fit_init, x, y, zoom)
        sigx = [fit.parameters[3]]
        sigy = [fit.parameters[4]]
        chi = [np.sum(np.square(zoom - fit(x, y)) / fit(x, y))]
        table = Table(data=[sigx, sigy, chi], names=('sigma_x', 'sigma_y', 'Chi^2'))
        return table
    else:
        return None

sig_table = Table(names=('sigma_x', 'sigma_y', 'Chi^2'))
for row in range(len(peak_tbl)):
    row_table = sigmas(image, int(round(peak_tbl[row]['xcentroid'])),
                       int(round(peak_tbl[row]['ycentroid'])),
                       peak_tbl[row]['peak'], bkg)
    if not row_table is None:
        sig_table.add_row(row_table[0])
sig_table.sort('Chi^2')


print sig_table

sigma_x = np.median(np.asarray(sig_table['sigma_x'][:20]))
sigma_y = np.median(np.asarray(sig_table['sigma_y'][:20]))


print [sigma_x, sigma_y]



zoom = 50
target = image[(image.shape[1]-zoom)/2:(image.shape[1]+zoom)/2, (image.shape[0]-zoom)/2:(image.shape[0]+zoom)/2]
ty, tx = np.mgrid[:zoom, :zoom]

target_peaks = find_peaks(target, threshold, box_size=2*zoom)

xpeak = target_peaks['x_peak']
ypeak = target_peaks['y_peak']

model_2 = models.Gaussian2D + models.Gaussian2D + models.Const2D
model_1 = models.Gaussian2D + models.Const2D

bound = 5


fit_init_2 = model_2(amplitude_0=np.max(target), x_mean_0=xpeak, y_mean_0=ypeak, x_stddev_0=sigma_x, y_stddev_0=sigma_y,
                 amplitude_1=np.max(target), x_mean_1=xpeak, y_mean_1=ypeak, x_stddev_1=sigma_x, y_stddev_1=sigma_y,
                 amplitude_2=bkg,
                 fixed={'x_stddev_0': True, 'y_stddev_0': True,
                        'x_stddev_1': True, 'y_stddev_1': True,
                        'amplitude_2': True},
                 bounds={'x_mean_0': (xpeak-bound, xpeak+bound), 'y_mean_0': (ypeak-bound, ypeak+bound),
                         'x_mean_1': (xpeak-bound, xpeak+bound), 'y_mean_1': (ypeak-bound, ypeak+bound)}
                 )

fit_init_1 = model_1(amplitude_0=np.max(target), x_mean_0=xpeak, y_mean_0=ypeak, x_stddev_0=sigma_x, y_stddev_0=sigma_y,
                 amplitude_1=bkg,
                 fixed={'x_stddev_0': True, 'y_stddev_0': True,
                        'amplitude_1': True},
                 bounds={'x_mean_0': (xpeak-bound, xpeak+bound), 'y_mean_0': (ypeak-bound, ypeak+bound)}
                 )


fit1 = fitter(fit_init_1, tx, ty, target)
residual1 = target - fit1(tx, ty)
flux1 = 2 * np.pi * sigma_x * sigma_y * (fit1.parameters[0])

fit2 = fitter(fit_init_2, tx, ty, target)
residual2 = target - fit2(tx, ty)
flux2 = 2 * np.pi * sigma_x * sigma_y * (fit2.parameters[0] + fit2.parameters[6])


print flux1
print flux2


# Plot the data with the best-fit model
plt.figure(figsize=(12, 5))
plt.suptitle('Single 2D Gaussian')
plt.subplot(1, 3, 1)
plt.imshow(target, origin='lower', interpolation='nearest', cmap='viridis')
plt.title("Data")
plt.subplot(1, 3, 2)
plt.imshow(fit1(tx, ty), origin='lower', interpolation='nearest', cmap='viridis')
plt.title("Model")
plt.subplot(1, 3, 3)
plt.imshow(residual1, origin='lower', interpolation='nearest', cmap='viridis')
plt.title("Residuals")

plt.figure(figsize=(12, 5))
plt.suptitle('Double 2D Gaussian')
plt.subplot(1, 3, 1)
plt.imshow(target, origin='lower', interpolation='nearest', cmap='viridis')
plt.title("Data")
plt.subplot(1, 3, 2)
plt.imshow(fit2(tx, ty), origin='lower', interpolation='nearest', cmap='viridis')
plt.title("Model")
plt.subplot(1, 3, 3)
plt.imshow(residual2, origin='lower', interpolation='nearest', cmap='viridis')
plt.title("Residuals")

fig = plt.figure(figsize=(9, 7))
fig.suptitle('Data')
ax = fig.add_subplot(111, projection='3d')
Axes3D.plot_surface(ax, tx, ty, target, cmap='viridis')

fig = plt.figure(figsize=(9, 7))
fig.suptitle('Model')
ax = fig.add_subplot(111, projection='3d')
Axes3D.plot_surface(ax, tx, ty, fit2(tx, ty), cmap='viridis')

fig = plt.figure(figsize=(9, 7))
fig.suptitle('Residuals')
ax = fig.add_subplot(111, projection='3d')
Axes3D.plot_surface(ax, tx, ty, residual2, cmap='viridis')

plt.show()
