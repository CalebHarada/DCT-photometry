
"""

LMI_Photometry
Caleb K. Harada, 2017
(The University of Maryland; The University of Chicago)

This module contains three functions that produce useful photometry from raw .FITS images from LMI at DCT:
-   'Data_Reduction' creates and applies a master bias, flat, and dark (optional) frame to science images, and updates
    the .FITS header to make targets Simbad-compatible.
-   'Aperture_Photometry' measures raw electron counts for a target star and utilizes the .FITS header to calculate and
    save fluxes and instrumental magnitudes. It flags specified stars as standards to be used for standard magnitude
    transformations.
-   'Convert_Magnitudes' reads the magnitudes and airmass values saved in the .FITS headers of standard stars,
    calculates a magnitude transformation for each filter used, then applies the transformation to science images to
    convert their instrumental magnitudes to standard magnitudes. It saves measurements and uncertainties in a .txt
    table in 'ascii' format.

contact: charada@umd.edu
"""


import ccdproc as cp
from ccdproc import ImageFileCollection, Combiner
from astropy import units as u
import os
from astroquery.simbad import Simbad
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from photutils import fit_2dgaussian, CircularAperture, aperture_photometry, MADStdBackgroundRMS, MMMBackground
import numpy as np
import matplotlib.pyplot as plt
import ast
from scipy.optimize import curve_fit

#######################################################################################################################

def Data_Reduction(directory, filters, targets, save_to=None, dark_exp=1.0, subtract_dark=False):

    '''

    :param directory: str
            A directory containing raw .FITS images and calibration frames
    :param filters: dict
            Filters used and corresponding flat exposures
            {'filter' : flat exposure}
    :param targets: dict
            "SCITARG" name in .FITS header and corresponding name in Simbad
            {'FITS target name' : 'Simbad target name'}
    :param save_to: str, opt
            Optional second directory to save calibrated frames to
    :param dark_exp: float
            Exposure time for dark frames, default is 1.0 sec
    :param subtract_dark: bool
            Set to True to subtract dark frame, default is False
            Note: LMI has negligible dark current
    :return: None
    '''

    pipeline_out = '%s\\pipeline_out' % directory
    if not os.path.exists(pipeline_out) and save_to == None:
        os.makedirs(pipeline_out)
        print '\'pipeline_out\' directory created'

    ifc = ImageFileCollection(location=directory, keywords='*')

    def make_bias():
        bias_frames = ifc.files_filtered(IMAGETYP='bias')
        print 'Calculating master bias frame...'
        bias_data = []
        for frame in bias_frames:
            bias_read = cp.CCDData.read('%s\\%s' % (directory, frame))
            bias_dev = cp.create_deviation(bias_read,
                                           gain=bias_read.header['GAIN'] * u.electron / u.adu,
                                           readnoise=bias_read.header['RDNOISE'] * u.electron)
            bias_gained = cp.gain_correct(bias_dev, bias_read.header['GAIN'] * u.electron / u.adu)
            bias_data.append(bias_gained)
        bias_comb = Combiner(bias_data)
        master_bias = bias_comb.median_combine()
        master_bias.header['IMAGETYP'] = 'MASTER BIAS'
        return master_bias

    if save_to == None:
        bias_name = '%s\\Master_Bias.fits' % pipeline_out
        if os.path.isfile(bias_name) is False:
            cp.fits_ccddata_writer(ccd_data=make_bias(), filename=bias_name)
            print 'Master bias created'
        else:
            print 'Master bias already exists'
    else:
        bias_name = '%s\\Master_Bias.fits' % save_to
        if os.path.isfile(bias_name) is False:
            cp.fits_ccddata_writer(ccd_data=make_bias(), filename=bias_name)
            print 'Master bias created'
        else:
            print 'Master bias already exists'

    def make_dark():
        dark_frames = ifc.files_filtered(IMAGETYP='dark', EXPTIME=dark_exp)
        print 'Calculating master dark frame...'
        dark_data = []
        for frame in dark_frames:
            dark_read = cp.CCDData.read('%s\\%s' % (directory, frame))
            dark_dev = cp.create_deviation(dark_read,
                                           gain=dark_read.header['GAIN'] * u.electron / u.adu,
                                           readnoise=dark_read.header['RDNOISE'] * u.electron)
            dark_gained = cp.gain_correct(dark_dev, dark_read.header['GAIN'] * u.electron / u.adu)
            master_bias = cp.CCDData.read(bias_name)
            bias_sub = cp.subtract_bias(dark_gained, master_bias)
            dark_data.append(bias_sub)
        dark_comb = Combiner(dark_data)
        master_dark = dark_comb.median_combine()
        master_dark.header['EXPTIME'] = (dark_exp, 'integration time, seconds')
        master_dark.header['IMAGETYP'] = 'MASTER DARK'
        return master_dark

    if subtract_dark == True:
        if save_to == None:
            dark_name = '%s\\Master_Dark.fits' % pipeline_out
            if os.path.isfile(dark_name) is False:
                cp.fits_ccddata_writer(ccd_data=make_dark(), filename=dark_name)
                print 'Master dark created'
            else:
                print 'Master dark already exists'
        else:
            dark_name = '%s\\Master_Dark.fits' % save_to
            if os.path.isfile(dark_name) is False:
                cp.fits_ccddata_writer(ccd_data=make_dark(), filename=dark_name)
                print 'Master dark created'
            else:
                print 'Master dark already exists'

    def make_flat(filter, flat_exp):
        flat_frames = ifc.files_filtered(IMAGETYP='sky flat', EXPTIME=flat_exp, FILTERS=filter)
        print 'Calculating %s flat field...' % filter
        flat_data = []
        for frame in flat_frames:
            flat_read = cp.CCDData.read('%s\\%s' % (directory, frame))
            flat_dev = cp.create_deviation(flat_read,
                                           gain=flat_read.header['GAIN'] * u.electron / u.adu,
                                           readnoise=flat_read.header['RDNOISE'] * u.electron)
            flat_gained = cp.gain_correct(flat_dev, flat_read.header['GAIN'] * u.electron / u.adu)
            master_bias = cp.CCDData.read(bias_name)
            bias_sub = cp.subtract_bias(flat_gained, master_bias)
            flat_data.append(bias_sub)
        flat_comb = Combiner(flat_data)
        flat_medcomb = flat_comb.median_combine()
        if subtract_dark == True:
            master_dark = cp.CCDData.read(dark_name)
            master_flat = cp.subtract_dark(flat_medcomb, master_dark,
                                           dark_exposure=master_dark.header['EXPTIME'] * u.second,
                                           data_exposure=flat_exp * u.second,
                                           scale=False)
        else:
            master_flat = flat_medcomb
        master_flat.header['EXPTIME'] = (flat_exp, 'integration time, seconds')
        master_flat.header['IMAGETYP'] = 'MASTER FLAT'
        master_flat.header['FILTERS'] = (filter, 'Composite Filter Name')
        return master_flat

    for filter, flat_exp in filters.items():
        if save_to == None:
            flat_name = '%s\\Master_Flat_%s.fits' % (pipeline_out, filter)
            if os.path.isfile(flat_name) is False:
                cp.fits_ccddata_writer(ccd_data=make_flat(filter, flat_exp), filename=flat_name)
                print 'Master %s flat created' % filter
            else:
                print 'Master %s flat already exists' % filter
        else:
            flat_name = '%s\\Master_Flat_%s.fits' % (save_to, filter)
            if os.path.isfile(flat_name) is False:
                cp.fits_ccddata_writer(ccd_data=make_flat(filter, flat_exp), filename=flat_name)
                print 'Master %s flat created' % filter
            else:
                print 'Master %s flat already exists' % filter


    def reduce(frame, target, simbadref, filter, epoch):
        print 'Processing %s (%s)...' % (simbadref, target)
        target_read = cp.CCDData.read('%s\\%s' % (directory, frame))
        target_dev = cp.create_deviation(target_read,
                                    gain=target_read.header['GAIN'] * u.electron / u.adu,
                                    readnoise=target_read.header['RDNOISE'] * u.electron)
        target_gained = cp.gain_correct(target_dev, target_read.header['GAIN'] * u.electron/u.adu)
        master_bias = cp.CCDData.read(bias_name)
        bias_sub = cp.subtract_bias(target_gained, master_bias)
        if subtract_dark == True:
            master_dark = cp.CCDData.read(dark_name)
            dark_sub = cp.subtract_dark(bias_sub, master_dark,
                                        dark_exposure=master_dark.header['EXPTIME'] * u.second,
                                        data_exposure=target_gained.header['EXPTIME'] * u.second,
                                        scale=False)
        else:
            dark_sub = bias_sub
        if save_to == None:
            flat_name = '%s\\Master_Flat_%s.fits' % (pipeline_out, filter)
            master_flat = cp.CCDData.read(flat_name)
        else:
            flat_name = '%s\\Master_Flat_%s.fits' % (save_to, filter)
            master_flat = cp.CCDData.read(flat_name)
        target_reduced = cp.flat_correct(dark_sub, master_flat,
                                          min_value=0.1)
        target_reduced.header['QUERNAM'] = (simbadref, 'Use to query Simbad')
        target_reduced.header['REDUCED'] = ('True', 'Indicates this CCD frame has been reduced')
        target_reduced.header['EPOCH'] = (epoch, 'Arbitrary epoch')
        target_reduced.header.rename_keyword('RADECSYS', 'RADESYSa')
        return target_reduced

    for filter, flat_exp in filters.items():
        for target, simbadref in targets.items():
            target_frames = ifc.files_filtered(SCITARG=target, FILTERS=filter)
            epoch = 0
            for frame in target_frames:
                if save_to == None:
                    target_name = '%s\\%s_%s_%s.fits' % (pipeline_out, simbadref, filter, epoch)
                    if os.path.isfile(target_name) is False:
                        cp.fits_ccddata_writer(ccd_data=reduce(frame, target, simbadref, filter, epoch), filename=target_name)
                        print '%s_%s_%s done' % (simbadref, filter, epoch)
                    else:
                        print '%s_%s_%s has already been reduced' % (simbadref, filter, epoch)
                else:
                    target_name = '%s\\%s_%s_%s.fits' % (save_to, simbadref, filter, epoch)
                    if os.path.isfile(target_name) is False:
                        cp.fits_ccddata_writer(ccd_data=reduce(frame, target, simbadref, filter, epoch), filename=target_name)
                        print '%s_%s_%s done' % (simbadref, filter, epoch)
                    else:
                        print '%s_%s_%s has already been reduced' % (simbadref, filter, epoch)
                epoch += 1

    print 'Data reduction done.'

#######################################################################################################################

def Aperture_Photometry(directory, ap_radius, standards, show_figures=False):

    '''

    :param directory: str
            A directory containing reduced .FITS images
    :param ap_radius: int
            Radius of aperture used for photometry
    :param standards: dict
            Simbad-compatible name with list of standard star names in the field
            {'Query Name' : ['Standard Query Names']}
    :param show_figures: bool
            Display optional figures that are relevant, default is False
    :return: None
    '''

    Simbad.add_votable_fields('coo(fk5)', 'propermotions')

    def count_electrons(file, query_name):

        hdu = fits.open(file)
        data = hdu[0].data
        wcs = WCS(hdu[0].header)
        targinfo = Simbad.query_object(query_name)

        print 'Finding target...'

        # Determine RA/Dec of target from Simbad query
        targinfo_RA = targinfo['RA_fk5'][0]
        targRA = [float(targinfo_RA[:2]), float(targinfo_RA[3:5]), float(targinfo_RA[6:])]
        RA = targRA[0] * 15 + targRA[1] / 4 + targRA[2] / 240
        dRA = targinfo['PMRA'][0] * (15 / 3600000.0)
        if dRA > 0:
            RA = RA + dRA
        targinfo_Dec = targinfo['DEC_fk5'][0]
        targDec = [float(targinfo_Dec[1:3]), float(targinfo_Dec[4:6]), float(targinfo_Dec[7:])]
        Dec = (targDec[0]) + targDec[1] / 60 + targDec[2] / 3600
        if targinfo_Dec[0] == '-':
            Dec = np.negative(Dec)
        dDec = targinfo['PMDEC'][0] * (15 / 3600000.0)
        if dDec > 0:
            Dec = Dec + dDec

        # Convert RA/Dec to pixels
        pix = wcs.all_world2pix(RA, Dec, 0)
        xpix = int(pix[0])
        ypix = int(pix[1])

        # Trim data to 90x90 pixels near target; fit 2D Gaussian to find center pixel of target
        centzoom = data[ypix - 45:ypix + 45, xpix - 45:xpix + 45]
        centroid = fit_2dgaussian(centzoom)
        xcent = xpix - 45 + int(centroid.x_mean.value)
        ycent = ypix - 45 + int(centroid.y_mean.value)

        if show_figures == True:
            plt.figure()
            plt.imshow(centzoom, origin='lower', vmin=np.median(data) - 100, vmax=np.median(data) + 400)
            plt.colorbar()
            plt.show()

        print 'Calculating flux...'

        # Draw aperture around target, sum pixel values, and calculate error
        aperture = CircularAperture((xcent, ycent), r=ap_radius)
        ap_table = aperture_photometry(data, aperture)
        ap_sum = ap_table['aperture_sum'][0]

        # print ap_table
        if show_figures == True:
            plt.figure()
            plt.imshow(data, origin='lower', interpolation='nearest', vmin=np.median(data) - 100, vmax=np.median(data) + 400)
            aperture.plot(color='red')
            plt.show()

        apzoom = data[ycent - 250:ycent + 250,
                 xcent - 250:xcent + 250]  # trim data to 400x400 region centered on target

        bkg = MMMBackground()(apzoom)
        std = MADStdBackgroundRMS()(apzoom)
        threshold = bkg + 5 * std

        # Find appropriate sky aperture, sum pixel values, calculate error
        def find_sky():

            rand_x = np.random.randint(0, 500)  # randomly select pixel in region
            rand_y = np.random.randint(0, 500)

            if rand_x in range(250 - 3 * ap_radius, 250 + 3 * ap_radius) \
                    or rand_y in range(250 - 3 * ap_radius, 250 + 3 * ap_radius):
                return find_sky()  # reselect pixels if aperture overlaps target

            elif rand_x not in range(2 * ap_radius, 500 - 2 * ap_radius) \
                    or rand_y not in range(2 * ap_radius, 500 - 2 * ap_radius):
                return find_sky()

            else:
                sky = CircularAperture((rand_x, rand_y), r=ap_radius)
                sky_table = aperture_photometry(apzoom, sky)
                sky_sum = sky_table['aperture_sum'][0]

                sky_x = int(sky_table['xcenter'][0].value)
                sky_y = int(sky_table['ycenter'][0].value)
                sky_zoom = apzoom[sky_y - ap_radius : sky_y + ap_radius, sky_x - ap_radius : sky_x + ap_radius]

                sky_avg = sky_sum / sky.area()  # reselect pixels if bright source is in aperture
                if np.max(sky_zoom) < threshold and sky_avg > 0:

                    if show_figures == True:
                        plt.figure()
                        plt.imshow(apzoom, origin='lower', interpolation='nearest', vmin=np.median(data) - 100,
                                   vmax=np.median(data) + 400)
                        sky.plot()
                        plt.show()

                    return sky_sum
                else:
                    return find_sky()

        print 'Estimating sky background...'

        # Calculate final electron count
        sample_size = 100
        list = np.arange(0, sample_size)
        sums = []
        for i in list:
            sky_sum = find_sky()
            final_sum = ap_sum - sky_sum
            sums.append(final_sum)

        electron_counts = np.mean(sums)
        elec_uncert = np.std(sums)

        hdu.close()
        return electron_counts, elec_uncert

    for file_name in os.listdir(directory):
        ext = os.path.splitext(file_name)[1]
        if ext == '.fits':
            file = '%s\\%s' % (directory, file_name)
        hdu = fits.open(file, mode='update')
        type = hdu[0].header['IMAGETYP']
        if type == 'OBJECT':
            exposure = hdu[0].header['EXPTIME']
            query_name = hdu[0].header['QUERNAM']
            if query_name in standards.keys():
                standard_list = []
                stdmag_list = []
                uncert_list = []
                hdu[0].header['STANDARD'] = ('YES', 'Indicates standard star')
                for standard in standards[query_name]:
                    standard_list.append(standard)
                    electron_counts, elec_uncert = count_electrons(file, standard)

                    print 'Calculating instrumental magnitude...'

                    ins_magnitude = -2.5 * np.log10(electron_counts / exposure) + 25
                    ins_magnitude_error = (2.5 * elec_uncert) / (electron_counts * np.log(10))
                    stdmag_list.append(round(ins_magnitude, 5))
                    uncert_list.append(round(ins_magnitude_error, 5))

                    print standard + ' done.'

                print 'Updating FITS header...'

                hdu[0].header['STDSTARS'] = (str(standard_list),
                                            'List of standard stars')
                hdu[0].header['OBJMAGS'] = (str(stdmag_list),
                                            'Instrumental mags')
                hdu[0].header['UNCERTS'] = (str(uncert_list),
                                            'Estimated uncertainties')

                hdu.close()

            else:
                hdu[0].header['STANDARD'] = ('NO', 'Indicates standard star')
                electron_counts, elec_uncert = count_electrons(file, query_name)

                print 'Calculating instrumental magnitude...'

                ins_magnitude = -2.5*np.log10(electron_counts/exposure) + 25
                ins_magnitude_error = (2.5*elec_uncert)/(electron_counts*np.log(10))

                print 'Updating FITS header...'

                hdu[0].header['OBJMAG'] = (round(ins_magnitude, 5), 'Target object\'s instrumental mag')
                hdu[0].header['u_OBJMAG'] = (round(ins_magnitude_error, 5), 'Estimated instrumental mag uncertainty')

                hdu.close()

                print query_name + ' done.'

    print 'Aperture photometry done.'

#######################################################################################################################

def Convert_Magnitudes(directory, filters, bin_size=10, show_figures=False):

    '''

    :param directory: str
            A directory containing reduced .FITS images instrumental magnitudes appended to the .FITS header
    :param filters: list
            Filters used
            ['filter 1', 'filter 2']
    :param bin_size: int
            Number of epochs target is observed, default is 10
    :param show_figures: bool
            Display optional figures that are relevant, default is False
    :return: None
    '''

    ifc = ImageFileCollection(location=directory, keywords='*')

    def read_science(file):
        hdu = fits.open(file)[0]
        ZA = hdu.header['ZA']

        target = hdu.header['QUERNAM']
        print 'Reading science target info...'
        epoch = hdu.header['EPOCH']
        filt = hdu.header['FILTERS']
        obsno = hdu.header['OBSERNO']
        airmass = 1/np.cos(ZA*(np.pi/180))
        mag_ins = hdu.header['OBJMAG']
        return target, epoch, filt, obsno, airmass, mag_ins

    table_out = Table(names=('Target', 'Filter', 'Standard_Mag', 'Error'),
                      dtype=('S24', 'S1', 'f8', 'f8'))
    for filter in filters:

        std_targs = ifc.files_filtered(IMAGETYP='OBJECT', FILTERS=filter, STANDARD='YES')
        Simbad.add_votable_fields('fluxdata(%s)' % filter)

        X = []
        uX = []
        Y = []
        uY = []

        fig = plt.figure(figsize=(10, 8))
        fig.suptitle('%s-Band Magnitude Transformation' % filter)
        ax = fig.add_subplot(111)
        ax.set_xlabel('Airmass')
        ax.set_ylabel('Instrumental - Standard (mag)')
        ax.grid(color='white', linestyle='-', linewidth=1, alpha=1)
        ax.patch.set_facecolor('#2B3856')
        ax.patch.set_alpha(0.1)
        for spine in ax.spines:
            ax.spines[spine].set_color('white')
        for axis in ('x', 'y'):
            ax.tick_params(axis=axis, color='white')

        standard_table = Table(names=('Target', 'Airmass', 'Airmass_uncert', 'Mag_ins', 'Mi_uncert', 'Mag_real', 'Mr_uncert'),
                               dtype=('S24', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

        std_table = Table(names=('target', 'obsno', 'airmass', 'mag_ins'),
                          dtype=('S24', 'f8', 'f8', 'f8'))
        for targ in std_targs:
            print 'Reading standard magnitudes...'
            file = '%s\\%s' % (directory, targ)
            hdu = fits.open(file)[0]
            ZA = hdu.header['ZA']
            obsno = hdu.header['OBSERNO']
            airmass = 1 / np.cos(ZA * (np.pi / 180))
            standards = ast.literal_eval(hdu.header['STDSTARS'])
            stdmags = ast.literal_eval(hdu.header['OBJMAGS'])
            i = 0
            while i < len(standards):
                std_table.add_row([standards[i], obsno, airmass, stdmags[i]])
                i += 1

        std_table_grouped = std_table.group_by('target')
        count = 1
        for group in std_table_grouped.groups:
            group_length = len(group)
            num_obs = group_length / bin_size
            while num_obs > 0:
                bin_group = group[0: bin_size]
                mean_airmass = np.mean(bin_group['airmass'])
                unc_airmass = np.std(bin_group['airmass'])
                mean_mag_ins = np.mean(bin_group['mag_ins'])
                unc_mag_ins = np.std(bin_group['mag_ins'])
                target = bin_group['target'][0]
                sim_table = Simbad.query_object(target)
                mag_real = sim_table['FLUX_%s' % filter].quantity.value[0]
                unc_mag_real = sim_table['FLUX_ERROR_%s' % filter].quantity.value[0]
                standard_table.add_row([target, mean_airmass, unc_airmass, mean_mag_ins, unc_mag_ins, mag_real, unc_mag_real])
                group.remove_rows(slice(0, bin_size))
                if count == 1:
                    ax.errorbar(x=mean_airmass, y=mean_mag_ins - mag_real, yerr=np.sqrt(unc_mag_ins**2 + unc_mag_real**2),
                                xerr=unc_airmass, fmt='o', color='#4863A0', ms=4, elinewidth=1, capsize=2, label='Data')
                else:
                    ax.errorbar(x=mean_airmass, y=mean_mag_ins - mag_real, yerr=np.sqrt(unc_mag_ins**2 + unc_mag_real**2),
                                xerr=unc_airmass, fmt='o', color='#4863A0', ms=4, elinewidth=1, capsize=2)
                count += 1
                X.append(mean_airmass)
                uX.append(unc_airmass)
                Y.append(mean_mag_ins - mag_real)
                uY.append(np.sqrt(unc_mag_ins ** 2 + unc_mag_real ** 2))
                num_obs -= 1

        X = np.asarray(X, dtype='float64')
        Y = np.asarray(Y, dtype='float64')
        uY = np.asarray(uY, dtype='float64')

        print 'Calculating %s-band standard magnitude transformation...' % filter
        def func(x, A, B):
            return A + B * x

        popt, pcov = curve_fit(f=func, xdata=X, ydata=Y, sigma=uY)

        xline = np.linspace(-1, 10, 1000)
        ax.plot(xline, func(xline, *popt), 'k--', label='Fit')

        ps = np.random.multivariate_normal(popt, pcov, 10000)
        ysample = np.asarray([func(xline, *pi) for pi in ps])
        lower1 = np.percentile(ysample, 15.87, axis=0)
        upper1 = np.percentile(ysample, 84.13, axis=0)
        lower2 = np.percentile(ysample, 2.28, axis=0)
        upper2 = np.percentile(ysample, 97.72, axis=0)

        ax.fill_between(xline, lower1, upper1, facecolor='#4863A0', alpha=0.15, zorder=2)
        ax.fill_between(xline, lower2, upper2, facecolor='#4863A0', alpha=0.15, zorder=3)

        plt.legend(loc=2)
        plt.xlim(1, 1.5)
        plt.ylim(np.min(Y) - 0.1, np.max(Y) + 0.1)

        if show_figures == True:
            print 'STANDARD STARS (%s):' % filter
            standard_table.pprint(max_lines=-1, max_width=-1)


        mag_ins_table = Table(names=('Target', 'Airmass', 'Airmass_std', 'Mag_ins', 'Mag_ins_std'),
                          dtype=('S24', 'f8', 'f8', 'f8', 'f8'))
        science_targs = ifc.files_filtered(IMAGETYP='OBJECT', FILTERS=filter, STANDARD='NO')
        science_table = Table(names=('Target', 'Epoch', 'Filter', 'OBSERNO', 'Airmass', 'Mag_ins'),
                      dtype=('S24', 'f8', 'S1', 'f8', 'f8', 'f8'))
        for targ in science_targs:
            file = '%s\\%s' % (directory, targ)
            target, epoch, filt, obsno, airmass, mag_ins = read_science(file)
            science_table.add_row([target, epoch, filt, obsno, airmass, mag_ins])
        table_by_target = science_table.group_by('Target')
        for group in table_by_target.groups:
            group_length = len(group)
            num_obs = group_length / bin_size
            while num_obs > 0:
                bin_group = group[0 : bin_size]
                mean_airmass = np.mean(bin_group['Airmass'])
                std_airmass = np.std(bin_group['Airmass'])
                mean_mag_ins = np.mean(bin_group['Mag_ins'])
                std_mag_ins = np.std(bin_group['Mag_ins'])
                target = bin_group['Target'][0]
                mag_ins_table.add_row([target, mean_airmass, std_airmass, mean_mag_ins, std_mag_ins])
                group.remove_rows(slice(0, bin_size))
                num_obs -= 1

        for row in mag_ins_table:
            target = row['Target']
            print 'Applying standard magnitude transformation to %s...' % target
            mi = row['Mag_ins']
            am = row['Airmass']
            mag_real = mi - func(am, *popt)
            samplemag = np.asarray([func(am, *pi) for pi in ps])
            uncert = func(am, *popt) - np.percentile(samplemag, 15.87, axis=0)
            mi_error = row['Mag_ins_std']
            mr_error = np.sqrt(mi_error**2 + uncert**2)
            table_out.add_row([target, filter, mag_real, mr_error])

        if show_figures == True:
            plt.show()

    print 'Saving measurements table...'
    name = '%s\\Measurements_Table.txt' % (directory)
    if os.path.isfile(name) is True:
        print 'Measurements table already exists.'
    else:
        table_out.write(name, format='ascii')
        print 'Measurements table saved.'

    if show_figures == True:
        print 'CONVERTED SCIENCE MAGS: '
        print table_out



