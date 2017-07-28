import matplotlib.pyplot as plt
import ccdproc as cp
from ccdproc import ImageFileCollection, Combiner
from astropy import units as u

directory = 'Data\\LMI_2015Jun02\\LMI.20150602'


ic = ImageFileCollection(location=directory, keywords='*')
    # Useful keywords: 'OBSERNO','IMAGETYP','SCITARG','FILTERS','EXPTIME'
bias_frames = ic.files_filtered(IMAGETYP='bias')
bias_data = []
for frame in bias_frames:
    data = cp.CCDData.read(directory + "\\" + frame)
    data_dev = cp.create_deviation(data,
                                   gain=data.header['GAIN'] * u.electron/u.adu,
                                   readnoise=data.header['RDNOISE'] * u.electron)
    data_gained = cp.gain_correct(data_dev, data.header['GAIN'] * u.electron/u.adu)
    bias_data.append(data_gained)
bias_comb = Combiner(bias_data)
master_bias = bias_comb.median_combine()

cp.fits_ccddata_writer(master_bias,'Data\\LMI_2015Jun02\\Master_Bias.fits')

plt.figure()
plt.imshow(master_bias, cmap='gray')
plt.colorbar()
plt.show()