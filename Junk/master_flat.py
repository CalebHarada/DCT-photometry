import matplotlib.pyplot as plt
import ccdproc as cp
from ccdproc import ImageFileCollection, Combiner
from astropy import units as u

directory = '...'
save_directory = '...'
filter = '...'
# set desired filter
master_bias = cp.CCDData.read('...')
master_dark = cp.CCDData.read('...')
# must already have these files created

#########################################################################

exp_time = 5                # Should not change from 1

ic = ImageFileCollection(location=directory,keywords='*')
    # Useful keywords: 'OBSERNO','IMAGETYP','SCITARG','FILTERS','EXPTIME'
flat_frames = ic.files_filtered(IMAGETYP='sky flat',EXPTIME='*',FILTERS=filter)
flat_data = []
for frame in flat_frames:
    data = cp.CCDData.read(directory + '\\' + frame)
    data_dev = cp.create_deviation(data,
                                   gain=data.header['GAIN'] * u.electron / u.adu,
                                   readnoise=data.header['RDNOISE'] * u.electron)
    data_gained = cp.gain_correct(data_dev, data.header['GAIN'] * u.electron / u.adu)
    bias_sub = cp.subtract_bias(data_gained, master_bias)
    flat_data.append(bias_sub)
flat_comb = Combiner(flat_data)
flat_comb = flat_comb.median_combine()
master_flat = cp.subtract_dark(flat_comb, master_dark,
                          dark_exposure=master_dark.header['EXPTIME'] * u.second,
                          data_exposure=exp_time * u.second,
                          scale=False)
master_flat.header['EXPTIME'] = (exp_time, 'integration time, seconds')
master_flat.header['IMAGETYP'] = 'MASTER FLAT'
master_flat.header['FILTERS'] = (filter, 'Composite Filter Name')

cp.fits_ccddata_writer(master_flat, save_directory)

plt.figure()
plt.imshow(master_flat, cmap='gray')
plt.colorbar()
plt.show()
