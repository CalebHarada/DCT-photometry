import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import os


filter = 'I'
directory = 'Data\\LMI_2015Jun02\\count_data\\standards\\aperture'

standards = Table.read('%s\\standards.csv' % directory)     # read standard magnitudes
#standards.pprint()


########################################################################################################################

def get_target_info(target_file):
    mags = []
    AMs = []
    table = Table.read('%s\\%s-band\\%s' % (directory, filter, target_file), format='ascii')
    target = table['Target'][0]
    for row in range(len(table)):
        mag = table[row]['Imag']
        ZA = table[row]['ZA']
        AM = 1/np.cos(ZA*(np.pi/180))
        mags.append(mag)
        AMs.append(AM)

    return [str(target), mags, AMs]

target_info = []
for file in os.listdir('%s\\%s-band' % (directory, filter)):
    ext = os.path.splitext(file)[1]
    if ext == '.txt':
        target_info.append(get_target_info(file))


x_values = []
y_values = []

plt.figure(figsize=(10,8))
fig = plt.axes()
plt.title(filter + '-band')
plt.xlabel('airmass')
plt.ylabel('m (standard) - m (instrumental)')

for info in target_info:
    standard_name = info[0]
    for row in range(len(standards)):
        if standards[row]['Simbad Name'] == standard_name:
            std_mag = standards[row][filter]
            V = standards[row]['V']

    ins_mag = info[1]
    airmass = info[2]

    x = airmass
    y = std_mag - ins_mag

    fig.plot(x, y, '.', label=standard_name)
    plt.legend()

    for i in x:
        x_values.append(i)
    for i in y:
        y_values.append(i)



x_array = np.asarray(x_values, dtype='float64')
y_array = np.asarray(y_values, dtype='float64')


delta = len(x_array)*np.sum(np.square(x_array)) - np.square(np.sum(x_array))
A = (np.sum(np.square(x_array))*np.sum(y_array) - np.sum(x_array)*np.sum(x_array*y_array)) / delta
B = (len(x_array)*np.sum(x_array*y_array) - np.sum(x_array)*np.sum(y_array)) / delta


yfit = A + B*x_array
plt.plot(x_array, yfit, 'm--')

print A
print B


plt.show()



