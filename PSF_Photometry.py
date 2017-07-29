import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
from photutils import find_peaks, DAOStarFinder, fit_2dgaussian, CircularAperture, IRAFStarFinder
from photutils.background import MADStdBackgroundRMS, MMMBackground
from astropy.table import Table
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D
from astroquery.simbad import Simbad
from astropy.wcs import WCS
import warnings
import os

#=======================================================================================================================

simbad_query = 'SA_104339' #Name on Simbad
directory = 'D:\\UChicago\\Reduced Data\\2015Jun02\\STANDARDS\\104334'  #source directory (line 462 -- save directory)
filter = 'I'

show_figs = 0       # 0 = show no images
save_resids = 0


#=======================================================================================================================

Simbad.add_votable_fields('coo(fk5)','propermotions')
def get_pixels(file_name, simbad_name):

    # Read in relevant information
    hdu = fits.open(file_name)
    wcs = WCS(hdu[0].header)
    targinfo = Simbad.query_object(simbad_name)

    # Determine RA/Dec of target from Simbad query
    targinfo_RA = targinfo['RA_fk5'][0]
    targRA = [float(targinfo_RA[:2]), float(targinfo_RA[3:5]), float(targinfo_RA[6:])]
    RAdeg = targRA[0]*15 + targRA[1]/4 + targRA[2]/240
    dRA = targinfo['PMRA'][0] * (15/3600000.0)      # minor correction for proper motion
    RA = RAdeg + dRA
    targinfo_Dec = targinfo['DEC_fk5'][0]
    targDec = [float(targinfo_Dec[1:3]), float(targinfo_Dec[4:6]), float(targinfo_Dec[7:])]
    Decdeg = (targDec[0]) + targDec[1]/60 + targDec[2]/3600
    if targinfo_Dec[0] == '-':
        Decdeg = np.negative(Decdeg)                # makes negative declinations negative
    dDec = targinfo['PMDEC'][0] * (15/3600000.0)
    Dec = Decdeg + dDec

    # Convert RA/Dec to pixels
    pix = wcs.all_world2pix(RA,Dec,0)
    xpix = int(pix[0])
    ypix = int(pix[1])

    #print xpix, ypix

    return xpix, ypix

#=======================================================================================================================

def create_model(file_name, box_size=100):

    hdu = fits.open(file_name)
    data = hdu[0].data
    error = hdu[2].data

    image = data[500:2500, 500:2500]
    u_image = error[500:2500, 500:2500]

    bkg = MMMBackground()(image)
    std = MADStdBackgroundRMS()(image)
    threshold = bkg + 4 * std

    daofind = DAOStarFinder(threshold=threshold, fwhm=2*np.sqrt(2*np.log(2)))
    iraffind = IRAFStarFinder(threshold=threshold, fwhm=2*np.sqrt(2*np.log(2)), minsep_fwhm=2.5, sharplo=0.5, sharphi=2.0, roundlo=0.0, roundhi=0.2, sky=bkg)
    peak_tbl = iraffind(image)
    circles = CircularAperture((peak_tbl['xcentroid'], peak_tbl['ycentroid']), r=30)

    #=========================================================================
    if show_figs == 1:
        plt.figure(figsize=(8, 6))
        plt.imshow(image, origin='lower', interpolation='nearest', cmap='viridis', vmin=np.median(image), vmax=np.median(image)+500)
        plt.title('Data')
        circles.plot(color='red')
        plt.show()
    #=========================================================================

    PSF_model = models.Gaussian2D + models.Gaussian2D + models.Gaussian2D + models.Const2D
    fitter = fitting.LevMarLSQFitter()

    def sigmas(x_mean, y_mean, peak, image=image, error=u_image, size=box_size):

        if size < y_mean < image.shape[0] - size and size < x_mean < image.shape[1] - size:

            zoom = image[y_mean - size / 2:y_mean + size / 2, x_mean - size / 2:x_mean + size / 2]
            u_zoom = error[y_mean - size / 2:y_mean + size / 2, x_mean - size / 2:x_mean + size / 2]
            sky = MMMBackground()(zoom)
            fit_init = PSF_model(amplitude_0=peak, x_mean_0=zoom.shape[1] / 2, y_mean_0=zoom.shape[0] / 2, theta_0=np.pi / 4,
                                 amplitude_1=0.75 * peak, x_mean_1=zoom.shape[1] / 2, y_mean_1=zoom.shape[0] / 2, theta_1=np.pi / 4,
                                 amplitude_2=0.5 * peak, x_mean_2=zoom.shape[1] / 2, y_mean_2=zoom.shape[0] / 2, theta_2=np.pi / 4,
                                 amplitude_3=sky,
                                 fixed={'amplitude_2': True},
                                 bounds={'x_mean_0': (zoom.shape[1] / 2 - 5, zoom.shape[1] / 2 + 5),
                                         'y_mean_0': (zoom.shape[1] / 2 - 5, zoom.shape[1] / 2 + 5),
                                         'x_stddev_0': (0, 5), 'y_stddev_0': (0, 5),
                                         'theta_0': (0, np.pi),
                                         'x_mean_1': (zoom.shape[1] / 2 - 5, zoom.shape[1] / 2 + 5),
                                         'y_mean_1': (zoom.shape[1] / 2 - 5, zoom.shape[1] / 2 + 5),
                                         'x_stddev_1': (0, 5), 'y_stddev_1': (0, 5),
                                         'theta_1': (0, np.pi),
                                         'x_mean_2': (zoom.shape[1] / 2 - 5, zoom.shape[1] / 2 + 5),
                                         'y_mean_2': (zoom.shape[1] / 2 - 5, zoom.shape[1] / 2 + 5),
                                         'x_stddev_2': (0, 5), 'y_stddev_2': (0, 5),
                                         'theta_2': (0, np.pi)}
                                 )



            def tie_x(model):
                xmean = model.x_mean_0
                return xmean

            def tie_y(model):
                ymean = model.y_mean_0
                return ymean

            def tie_theta(model):
                theta = model.theta_0
                return theta

            #fit_init.x_mean_1.tied = tie_x
            #fit_init.y_mean_1.tied = tie_y
            #fit_init.theta_1.tied = tie_theta

            y, x = np.mgrid[:zoom.shape[0], :zoom.shape[1]]
            fit = fitter(fit_init, x, y, zoom)
            sigx1 = [fit.x_stddev_0.value]
            sigy1 = [fit.y_stddev_0.value]
            sigx2 = [fit.x_stddev_1.value]
            sigy2 = [fit.y_stddev_1.value]
            sigx3 = [fit.x_stddev_2.value]
            sigy3 = [fit.y_stddev_2.value]
            theta1 = [fit.theta_0.value]
            theta2 = [fit.theta_1.value]
            theta3 = [fit.theta_2.value]
            star = zoom[zoom.shape[0] / 2 - 15:zoom.shape[0] / 2 + 15, zoom.shape[1] / 2 - 15:zoom.shape[1] / 2 + 15]
            star_fit = fit(x, y)[zoom.shape[0] / 2 - 15:zoom.shape[0] / 2 + 15, zoom.shape[1] / 2 - 15:zoom.shape[1] / 2 + 15]
            chi = [np.sum(np.square(star - star_fit) / star)]
            table = Table(data=[sigx1, sigy1, sigx2, sigy2, sigx3, sigy3, theta1, theta2, theta3, chi],
                          names=('sigma_x1', 'sigma_y1', 'sigma_x2', 'sigma_y2', 'sigma_x3', 'sigma_y3', 'theta_1', 'theta_2', 'theta_3', 'Chi^2'))
            if show_figs == 2:
                plt.figure()
                plt.imshow(star, origin='lower', interpolation='nearest')
                plt.figure()
                plt.imshow(zoom, origin='lower', interpolation='nearest', vmin=np.median(zoom), vmax=np.median(zoom) + 300)
                plt.figure()
                plt.imshow(zoom - fit(x, y), origin='lower', interpolation='nearest', vmin=np.median(zoom - fit(x,y)), vmax=np.median(zoom - fit(x,y)) + 300)
                table.pprint(max_lines=-1, max_width=-1)
                plt.show()

            return table

        else:
            return None

    sig_table = Table(names=('sigma_x1', 'sigma_y1', 'sigma_x2', 'sigma_y2', 'sigma_x3', 'sigma_y3', 'theta_1', 'theta_2', 'theta_3', 'Chi^2'))
    for row in range(len(peak_tbl)):
        if peak_tbl[row]['peak'] < 100000:
            row_table = sigmas(int(round(peak_tbl[row]['xcentroid'])),
                               int(round(peak_tbl[row]['ycentroid'])),
                               peak_tbl[row]['peak'])
            if not row_table is None:
                sig_table.add_row(row_table[0])
    sig_table.sort('Chi^2')

    sigma_x1 = np.median(np.asarray(sig_table['sigma_x1'][:5]))
    sigma_y1 = np.median(np.asarray(sig_table['sigma_y1'][:5]))
    theta_1 = np.median(np.asarray(sig_table['theta_1'][:5]))
    sigma_x2 = np.median(np.asarray(sig_table['sigma_x2'][:5]))
    sigma_y2 = np.median(np.asarray(sig_table['sigma_y2'][:5]))
    theta_2 = np.median(np.asarray(sig_table['theta_2'][:5]))
    sigma_x3 = np.median(np.asarray(sig_table['sigma_x3'][:5]))
    sigma_y3 = np.median(np.asarray(sig_table['sigma_y3'][:5]))
    theta_3 = np.median(np.asarray(sig_table['theta_3'][:5]))


    #===========================================
    if show_figs == 2:
        sig_table.pprint(max_lines=-1, max_width=-1)
    #===========================================

    return {'sig_x1': sigma_x1, 'sig_y1': sigma_y1, 'theta_1': theta_1,
            'sig_x2': sigma_x2, 'sig_y2': sigma_y2, 'theta_2': theta_2,
            'sig_x3': sigma_x3, 'sig_y3': sigma_y3, 'theta_3': theta_3}

#=======================================================================================================================

def psf_fit(file_name, x_coord, y_coord, sigma_x1, sigma_y1, theta_1,
            sigma_x2=1, sigma_y2=1, theta_2=1, sigma_x3=1, sigma_y3=1, theta_3=1, box_size=100):

    hdu = fits.open(file_name)
    data = hdu[0].data
    error = hdu[2].data

    epoch = hdu[0].header['OBSERNO']

    centzoom = data[y_coord - box_size / 2:y_coord + box_size / 2,
               x_coord - box_size / 2:x_coord + box_size / 2]
    centroid = fit_2dgaussian(centzoom)
    xcent = x_coord - box_size / 2 + int(centroid.x_mean.value)
    ycent = y_coord - box_size / 2 + int(centroid.y_mean.value)

    image = data[ycent - 100:ycent + 100, xcent - 100:xcent + 100]
    u_image = error[ycent - 100:ycent + 100, xcent - 100:xcent + 100]

    target = image[(image.shape[1] - box_size) / 2:(image.shape[1] + box_size) / 2,
             (image.shape[0] - box_size) / 2:(image.shape[0] + box_size) / 2]
    u_target = u_image[(image.shape[1] - box_size) / 2:(image.shape[1] + box_size) / 2,
             (image.shape[0] - box_size) / 2:(image.shape[0] + box_size) / 2]
    ty, tx = np.mgrid[:box_size, :box_size]

    bkg = MMMBackground()(target)
    std = MADStdBackgroundRMS()(target)
    threshold = bkg + 4 * std

    target_peaks = find_peaks(target, threshold, box_size=box_size)

    xpeak = target_peaks['x_peak']
    ypeak = target_peaks['y_peak']


    fitter = fitting.LevMarLSQFitter()
    bound = 5
    model_3 = models.Gaussian2D(amplitude=np.max(target), x_mean=xpeak, y_mean=ypeak,
                                x_stddev=sigma_x1, y_stddev=sigma_y1, theta=theta_1,
                                fixed={'x_stddev': True, 'y_stddev': True, 'theta': True},
                                bounds={'x_mean': (xpeak - bound, xpeak + bound),
                                        'y_mean': (ypeak - bound, ypeak + bound)})\
              + models.Gaussian2D(amplitude=0.75*np.max(target), x_mean=xpeak, y_mean=ypeak,
                                  x_stddev=sigma_x2, y_stddev=sigma_y2, theta=theta_2,
                                fixed={'x_stddev': True, 'y_stddev': True, 'theta': True},
                                bounds={'x_mean': (xpeak - bound, xpeak + bound),
                                        'y_mean': (ypeak - bound, ypeak + bound)})\
              + models.Gaussian2D(amplitude=0.5*np.max(target), x_mean=xpeak, y_mean=ypeak,
                                  x_stddev=sigma_x3, y_stddev=sigma_y3, theta=theta_3,
                                fixed={'x_stddev': True, 'y_stddev': True, 'theta': True},
                                bounds={'x_mean': (xpeak - bound, xpeak + bound),
                                        'y_mean': (ypeak - bound, ypeak + bound)})\
              + models.Const2D(amplitude=bkg, fixed={'amplitude': True})

    model_2 = models.Gaussian2D(amplitude=np.max(target), x_mean=xpeak, y_mean=ypeak,
                                x_stddev=sigma_x1, y_stddev=sigma_y1, theta=theta_1,
                                fixed={'x_stddev': True, 'y_stddev': True, 'theta': True},
                                bounds={'x_mean': (xpeak - bound, xpeak + bound),
                                        'y_mean': (ypeak - bound, ypeak + bound)})\
              + models.Gaussian2D(amplitude=0.5*np.max(target), x_mean=xpeak, y_mean=ypeak,
                                  x_stddev=sigma_x2, y_stddev=sigma_y2, theta=theta_2,
                                fixed={'x_stddev': True, 'y_stddev': True, 'theta': True},
                                bounds={'x_mean': (xpeak - bound, xpeak + bound),
                                        'y_mean': (ypeak - bound, ypeak + bound)})\
              + models.Const2D(amplitude=bkg, fixed={'amplitude': True})

    model_1 = models.Gaussian2D(amplitude=np.max(target), x_mean=xpeak, y_mean=ypeak,
                                x_stddev=sigma_x1, y_stddev=sigma_y1, theta=theta_1,
                                fixed={'x_stddev': True, 'y_stddev': True, 'theta': True},
                                bounds={'x_mean': (xpeak - bound, xpeak + bound),
                                        'y_mean': (ypeak - bound, ypeak + bound)})\
              + models.Const2D(amplitude=bkg, fixed={'amplitude': True})

    def tie_x(model):
        xmean = model.x_mean_0
        return xmean

    def tie_y(model):
        ymean = model.y_mean_0
        return ymean

    def tie_theta(model):
        theta = model.theta_0
        return theta

    #model_2.x_mean_1.tied = tie_x
    #model_2.y_mean_1.tied = tie_y
    #model_2.theta_1.tied = tie_theta


    fit1 = fitter(model_1, tx, ty, target)
    residual1 = target - fit1(tx, ty)
    flux1 = round(2 * np.pi * sigma_x1 * sigma_y1 * fit1.amplitude_0.value, 4)

    fit2 = fitter(model_2, tx, ty, target)
    residual2 = target - fit2(tx, ty)
    flux2 = round(2 * np.pi * (sigma_x1 * sigma_y1 * fit2.amplitude_0.value +
                               sigma_x2 * sigma_y2 * fit2.amplitude_1.value), 4)

    fit3 = fitter(model_3, tx, ty, target)
    residual3 = target - fit3(tx, ty)
    flux3 = round(2 * np.pi * (sigma_x1 * sigma_y1 * fit3.amplitude_0.value +
                               sigma_x2 * sigma_y2 * fit3.amplitude_1.value +
                               sigma_x3 * sigma_y3 * fit3.amplitude_2.value), 4)

    #====================================================================================

    if show_figs == 4:
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

        plt.figure(figsize=(12, 5))
        plt.suptitle('Triple 2D Gaussian')
        plt.subplot(1, 3, 1)
        plt.imshow(target, origin='lower', interpolation='nearest', cmap='viridis')
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(fit3(tx, ty), origin='lower', interpolation='nearest', cmap='viridis')
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(residual3, origin='lower', interpolation='nearest', cmap='viridis')
        plt.title("Residuals")

        #==================================================================================

        fig1 = plt.figure(figsize=(9, 7))
        fig1.suptitle('Single 2D Gaussian Residuals')
        ax = fig1.add_subplot(111, projection='3d')
        ax.set_zlim(-4000, 4000)
        Axes3D.plot_surface(ax, tx, ty, residual1, cmap='viridis')

        fig2 = plt.figure(figsize=(9, 7))
        fig2.suptitle('Double 2D Gaussian Residuals')
        ax = fig2.add_subplot(111, projection='3d')
        ax.set_zlim(-4000, 4000)
        Axes3D.plot_surface(ax, tx, ty, residual2, cmap='viridis')

        fig3 = plt.figure(figsize=(9, 7))
        fig3.suptitle('Triple 2D Gaussian Residuals')
        ax = fig3.add_subplot(111, projection='3d')
        ax.set_zlim(-4000, 4000)
        Axes3D.plot_surface(ax, tx, ty, residual3, cmap='viridis')

        plt.show()

    if save_resids == 1:
        fig3 = plt.figure(figsize=(9, 7))
        fig3.suptitle(simbad_query + ' Triple Gaussian Model Residuals')
        ax = fig3.add_subplot(111, projection='3d')
        ax.set_zlim(-4000, 4000)
        Axes3D.plot_surface(ax, tx, ty, residual3, cmap='viridis')
        plt.savefig('Data\\LMI_2015Jun02\\Residuals\\' + simbad_query + '_obsno' + str(epoch) + '.png')

    #====================================================================================

    chi_fit1 = round(np.sum(np.square(residual1 / u_target)) / (box_size**2), 4)
    chi_fit2 = round(np.sum(np.square(residual2 / u_target)) / (box_size**2), 4)
    chi_fit3 = round(np.sum(np.square(residual3 / u_target)) / (box_size**2), 4)

    return {'flux1': flux1, 'flux2': flux2, 'flux3': flux3,
            'chi1': chi_fit1, 'chi2': chi_fit2, 'chi3': chi_fit3}

#=======================================================================================================================

name_list = []
number_list = []
time_list = []
HA_list = []
ZA_list = []
air_list = []
filter_list = []
exp_list = []

f1_list = []
f2_list = []
f3_list = []

i1_list = []
i2_list = []
i3_list = []

c1_list = []
c2_list = []
c3_list = []


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for name in os.listdir(directory):

        ext = os.path.splitext(name)[1]
        if ext == '.fits':

            file = directory + '\\' + name
            hdu = fits.open(file)[0]

            if hdu.header['FILTERS'] == filter:

                pixels = get_pixels(file_name=file, simbad_name=simbad_query)
                model = create_model(file_name=file, box_size=200)
                photometry = psf_fit(file_name=file, x_coord=pixels[0], y_coord=pixels[1],
                                     sigma_x1=model['sig_x1'], sigma_y1=model['sig_y1'], theta_1=model['theta_1'],
                                     sigma_x2=model['sig_x2'], sigma_y2=model['sig_y2'], theta_2=model['theta_2'],
                                     sigma_x3=model['sig_x3'], sigma_y3=model['sig_y3'], theta_3=model['theta_3'],
                                     box_size=80)

                targname = simbad_query
                number = hdu.header['OBSERNO']
                time = hdu.header['UT']
                hour_angle = hdu.header['HA']
                zenith_angle = hdu.header['ZA']
                airmass = hdu.header["AIRMASS"]
                exposure = hdu.header['EXPTIME']

                i_mag1 = round(-2.5 * np.log10(photometry['flux1'] / exposure) + 25, 4)
                i_mag2 = round(-2.5 * np.log10(photometry['flux2'] / exposure) + 25, 4)
                i_mag3 = round(-2.5 * np.log10(photometry['flux3'] / exposure) + 25, 4)

                name_list.append(str(targname))
                number_list.append(number)
                time_list.append(time)
                HA_list.append(hour_angle)
                ZA_list.append(zenith_angle)
                air_list.append(airmass)
                filter_list.append(filter)
                exp_list.append(exposure)

                f1_list.append(photometry['flux1'])
                f2_list.append(photometry['flux2'])
                f3_list.append(photometry['flux3'])

                i1_list.append(i_mag1)
                i2_list.append(i_mag2)
                i3_list.append(i_mag3)

                c1_list.append(photometry['chi1'])
                c2_list.append(photometry['chi2'])
                c3_list.append(photometry['chi3'])

                print i_mag3


columns = 'Target', 'ObsNo', 'UT', 'HA', 'ZA', 'AirMass', 'Filter', 'IntTime', 'Flux', 'Imag'
data = [name_list, number_list, time_list, HA_list, ZA_list, air_list, filter_list, exp_list, f3_list, i3_list]
data_table = Table(data=data, names=columns, meta={'name': simbad_query})

table_name = 'Data\\LMI_2015Jun02\\count_data\\standards\\%s-band\\%s_%s_data.txt' % (filter, simbad_query, filter)
if os.path.isfile(table_name) is True:
    print 'Data table already exists for the target \'%s\'' % simbad_query
else:
    data_table.write(table_name, format='ascii')

data_table.show_in_browser(jsviewer=True)
