import matplotlib.pyplot as plt
import ccdproc as cp
from ccdproc import ImageFileCollection, Combiner
from astropy import units as u


### Set desired exposure time and the source directory to create master dark
    # NOTE: Since LMI dark current is negligible, only need to use one master dark (regardless of integration time)
exp_time = 10

directory = 'Data\\LMI_2015Jun02\\LMI.20150602'

master_bias = cp.CCDData.read('Data\\LMI_2015Jun02\\Master_Bias.fits')               # There is only one master bias frame
ic = ImageFileCollection(location=directory,keywords='*')
    # Useful keywords: 'OBSERNO','IMAGETYP','SCITARG','FILTERS','EXPTIME'
dark_frames = ic.files_filtered(IMAGETYP='dark', EXPTIME=exp_time)
dark_data = []
for frame in dark_frames:
    data = cp.CCDData.read(directory + '\\' + frame)
    data_dev = cp.create_deviation(data,
                                   gain=data.header['GAIN'] * u.electron / u.adu,
                                   readnoise=data.header['RDNOISE'] * u.electron)
    data_gained = cp.gain_correct(data_dev, data.header['GAIN'] * u.electron / u.adu)
    bias_sub = cp.subtract_bias(data_gained, master_bias)
    dark_data.append(bias_sub)
dark_comb = Combiner(dark_data)
master_dark = dark_comb.median_combine()
master_dark.header['EXPTIME'] = (exp_time, 'integration time, seconds')
master_dark.header['IMAGETYP'] = 'MASTER DARK'

name = 'Data\\LMI_2015Jun02\\Master_Dark_e%s.fits' % str(exp_time)
cp.fits_ccddata_writer(master_dark, name)

plt.figure()
plt.imshow(master_dark, cmap='gray')
plt.colorbar()
plt.show()
