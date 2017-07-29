from astroquery.simbad import Simbad
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
from astropy.table import Table
from photutils import fit_2dgaussian, CircularAperture, CircularAnnulus, aperture_photometry, find_peaks
import numpy as np
import os
import matplotlib.pyplot as plt



target_simbad = 'SA_104339'     # name in Simbad
directory = 'D:\\UChicago\\Reduced Data\\2015Jun02\\STANDARDS\\104334'
filter = 'R'        # read from FITS header


########################################################################################################################

Simbad.add_votable_fields('coo(fk5)','propermotions')
def calc_electrons(file, simbad):

    # Read in relevant information
    hdu = fits.open(file)
    data = hdu[0].data
    error = hdu[2].data
    wcs = WCS(hdu[0].header)
    targinfo = Simbad.query_object(simbad)

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
    xpix = int(pix[0])      # adding constants because reference pixel appears to be off (regardless of target)
    ypix = int(pix[1])

    # Trim data to 100x100 pixels near target; fit 2D Gaussian to find center pixel of target
    centzoom = data[ypix-45:ypix+45, xpix-45:xpix+45]
    centroid = fit_2dgaussian(centzoom)
    xcent = xpix - 45 + int(centroid.x_mean.value)
    ycent = ypix - 45 + int(centroid.y_mean.value)

    #plt.figure()
    #plt.imshow(centzoom, origin='lower', cmap='gray', vmin=np.median(data), vmax=np.median(data) + 200)
    #plt.colorbar()
    #plt.show()

    # Find max pixel value in zoomed area, median of data, and median absolute deviation of data
    peak = np.max(centzoom)
    median = np.median(data)            # global background estimate
    sigma = mad_std(data)               # like std of background, but more resilient to outliers

    # Find an appropriate aperture radius
    radius = 1
    an_mean = peak          # peak is just a starting value that will always be greater than the median
    while an_mean > median + sigma:
        annulus = CircularAnnulus((xcent, ycent), r_in=radius, r_out=radius + 1)
        an_sum = aperture_photometry(data,annulus)
        an_mean = an_sum['aperture_sum'][0] / annulus.area()
        radius += 1         # radius is selected once mean pix value inside annulus is within 2 sigma of median

    radius = 35

    # Draw aperture around target, sum pixel values, and calculate error
    aperture = CircularAperture((xcent, ycent), r=radius)
    ap_table = aperture_photometry(data, aperture, error=error)
    ap_sum = ap_table['aperture_sum'][0]
    ap_error = ap_table['aperture_sum_err'][0]

    #print ap_table
    #plt.imshow(data, origin='lower', interpolation='nearest', vmin=np.median(data), vmax=np.median(data) + 200)
    #aperture.plot()
    #plt.show()

    # Find appropriate sky aperture, sum pixel values, calculate error
    def find_sky():

        apzoom = data[ycent-250:ycent+250, xcent-250:xcent+250]     # trim data to 400x400 region centered on target
        errorzoom = error[ycent-250:ycent+250, xcent-250:xcent+250]
        rand_x = np.random.randint(0,500)                           # randomly select pixel in region
        rand_y = np.random.randint(0,500)

        if rand_x in range(250-3*radius,250+3*radius)\
            or rand_y in range(250-3*radius,250+3*radius):
            return find_sky()                                       # reselect pixels if aperture overlaps target

        elif rand_x not in range(2*radius, 500-2*radius)\
            or rand_y not in range(2*radius, 500-2*radius):
            return find_sky()

        else:
            sky = CircularAperture((rand_x,rand_y), r=radius)
            sky_table = aperture_photometry(apzoom, sky, error=errorzoom)
            sky_sum = sky_table['aperture_sum'][0]
            sky_error = sky_table['aperture_sum_err'][0]

            sky_x = int(sky_table['xcenter'][0].value)
            sky_y = int(sky_table['ycenter'][0].value)
            sky_zoom = apzoom[sky_y-radius:sky_y+radius, sky_x - radius:sky_x + radius]

            sky_avg = sky_sum/sky.area()                            # reselect pixels if bright source is in aperture
            if np.max(sky_zoom) < median + 5*sigma and sky_avg > 0:

                #plt.imshow(apzoom, origin='lower', interpolation='nearest', vmin=np.median(data), vmax=np.median(data) + 200)
                #sky.plot()
                #plt.show()

                return [sky_sum, sky_error]
            else:
                return find_sky()

    # Calculate final electron count with uncertainty
    sample_size = 100
    list = np.arange(0,sample_size)
    sums = []                           # list source-sky value of each iteration
    errors = []
    for i in list:
        final_sum = ap_sum - find_sky()[0]
        sums.append(final_sum)
        final_error = ap_error + find_sky()[1]
        errors.append(final_error)

    electron_counts = np.mean(sums)
    uncert = np.std(sums)

    return [electron_counts, uncert, sums, errors]    # return mean value of source-sky and propagated error


name_list = []
number_list = []
time_list = []
HA_list = []
ZA_list = []
air_list = []
filter_list = []
exp_list = []
electon_list = []
elec_error_list = []
mag_list = []
mag_error_list = []

all_counts_I = []
all_errors_I = []
Iexp = []
all_counts_R = []
all_errors_R = []
Rexp = []

# Return counts and propagated error from each frame in target directory
for name in os.listdir(directory):
    ext = os.path.splitext(name)[1]
    if ext == '.fits':
        file = directory + '\\' + name
        hdu = fits.open(file)[0]

        if hdu.header['FILTERS'] == filter:

            targname = target_simbad
            number = hdu.header['OBSERNO']
            time = hdu.header['UT']
            hour_angle = hdu.header['HA']
            zenith_angle = hdu.header['ZA']
            airmass = hdu.header["AIRMASS"]
            filter = hdu.header['FILTERS']
            exposure = hdu.header['EXPTIME']
            result = calc_electrons(file=file, simbad=target_simbad)
            electrons = int(round(result[0]))
            elec_error = int(round(result[1]))
            ins_magnitude = round(-2.5*np.log10(electrons/exposure) + 25, 5)
            ins_magnitude_error = round((2.5*elec_error)/(electrons*np.log(10)), 5)

            name_list.append(str(targname))
            number_list.append(number)
            time_list.append(time)
            HA_list.append(hour_angle)
            ZA_list.append(zenith_angle)
            air_list.append(airmass)
            filter_list.append(filter)
            exp_list.append(exposure)
            electon_list.append(electrons)
            elec_error_list.append(elec_error)
            mag_list.append(ins_magnitude)
            mag_error_list.append(ins_magnitude_error)

            print number, filter, electrons, elec_error, ins_magnitude, ins_magnitude_error

# Put data in table and save in target directory
columns = 'Target', 'ObsNo', 'UT', 'HA', 'ZA', 'AirMass', 'Filter', 'IntTime', 'IntCounts', 'IC_error', 'Imag', 'IM_error'
data = [name_list, number_list, time_list, HA_list, ZA_list, air_list, filter_list, exp_list, electon_list,
        elec_error_list, mag_list, mag_error_list]
data_table = Table(data=data, names=columns, meta={'name': target_simbad})
data_table.show_in_browser(jsviewer=True)
table_name = 'C:\\Users\\caleb\\documents\\PYTHON\\U_Chicago_Research\\Data\\LMI_2015Jun02\\count_data\\standards\\aperture\\%s-band\\%s_%s_data.txt' % (filter, target_simbad, filter)
if os.path.isfile(table_name) is True:
    print 'Data table already exists for the target \'%s\'' % target_simbad
else:
    data_table.write(table_name, format='ascii')









'''
import matplotlib.pyplot as plt
RA = 139.4375
Dec = 46.2069
file = 'D:\\UChicago\\Reduced Data\\2015Jun02\\Part 1\\1RXS_J091744.5+461229\\lmi.1RXS_J091744.5+461229_91_I.fits'
hdu = fits.open(file)
data = hdu[0].data
wcs = WCS(hdu[0].header)
pix = wcs.all_world2pix(RA,Dec,0)
xpix = int(pix[0]) + 5      # adding constants because reference pixel appears to be off (regardless of target)
ypix = int(pix[1]) - 30

# Trim data to 100x100 pixels near target; fit 2D Gaussian to find center pixel of target
centzoom = data[ypix-50:ypix+50, xpix-50:xpix+50]

plt.figure()
plt.imshow(centzoom, origin='lower', cmap='gray', vmin=np.median(data), vmax=np.median(data)+200)
plt.colorbar()
plt.show()
'''

