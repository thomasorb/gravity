# Grossir, Chazot 2019
import scipy
import numpy as np

altitude = [
    1e3,
    20e3,
    40e3,
    60e3,
    80e3,
    100e3,
    120e3
    ]

density = [
    1.2e-2,
    3e-3,
    3e-4,
    3e-5,
    2.5e-6,
    2e-7,
    1.5e-8
]

temperature = [
    240,
    190,
    160,
    165,
    155,
    140,
    120
    ]

pressure = [
    5e2,
    1e2,
    1e1,
    8e-1,
    7e-2,
    4e-3,
    3e-4
]

def log_interp1d(xx, yy, kind='linear'):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind, bounds_error=False, fill_value='extrapolate')
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

fdens = log_interp1d(altitude, density, kind='linear')
ftemp = scipy.interpolate.interp1d(altitude, np.array(temperature) - 273.15, kind='linear', bounds_error=False, fill_value='extrapolate')
fpres = log_interp1d(altitude, pressure, kind='linear')