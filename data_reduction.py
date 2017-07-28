import ccdproc as cp
from ccdproc import ImageFileCollection
from astropy import units as u
import os


source_directory = '...'
save_directory = '...'
bias = cp.CCDData.read(...)
dark = cp.CCDData.read(...)
Iflat = cp.CCDData.read(...)
Rflat = cp.CCDData.read(...)
# must have master flats, darks, and bias already in directory

################################################################################

ifc = ImageFileCollection(location=source_directory,keywords=['SCITARG'])
    # Useful keywords: 'OBSERNO','IMAGETYP','SCITARG','FILTERS','EXPTIME'
targ_names = ifc.values(keyword='SCITARG')
targets = []
for name in targ_names:
    if not name in targets:
        targets.append(name)
for target in targets:
    if not target == 'unspecified'\
            and not target == 'sky flats'\
            and not target == 'moon':
        science_fits = ifc.files_filtered(SCITARG=target)
        for fits in science_fits:
            data = cp.CCDData.read('%s\\%s' % (source_directory, fits))
            data = cp.create_deviation(data,
                                        gain=data.header['GAIN'] * u.electron / u.adu,
                                        readnoise=data.header['RDNOISE'] * u.electron)
            data_gained = cp.gain_correct(data, data.header['GAIN'] * u.electron/u.adu)
            bias_sub = cp.subtract_bias(data_gained, bias)
            dark_sub = cp.subtract_dark(bias_sub, dark,
                                        dark_exposure=dark.header['EXPTIME'] * u.second,
                                        data_exposure=data.header['EXPTIME'] * u.second,
                                        scale=False)

            if not os.path.isdir('%s\\%s' % (save_directory, target)):
                os.makedirs('%s\\%s' % (save_directory, target))

            if data.header['FILTERS'] == 'I':
                title = '%s\\%s\\lmi.%s_%s_I.fits' % (save_directory, target, target, data.header['OBSERNO'])
                if os.path.isfile(title) is False:
                    science_reduced = cp.flat_correct(dark_sub, Iflat,
                                                      min_value=0.1)
                    science_reduced.header['REDUCED'] = ('True', 'Indicates this CCD frame has been reduced')
                    science_reduced.header.rename_keyword('RADECSYS', 'RADESYSa')
                    cp.fits_ccddata_writer(science_reduced, title)
                else:
                    pass
            elif data.header['FILTERS'] == 'R':
                title = '%s\\%s\\lmi.%s_%s_R.fits' % (save_directory, target, target, data.header['OBSERNO'])
                if os.path.isfile(title) is False:
                    science_reduced = cp.flat_correct(dark_sub, Rflat,
                                                      min_value=0.1)
                    science_reduced.header['REDUCED'] = ('True', 'Indicates this CCD frame has been reduced')
                    science_reduced.header.rename_keyword('RADECSYS', 'RADESYSa')
                    cp.fits_ccddata_writer(science_reduced, title)
                else:
                    pass
            else:
                pass
