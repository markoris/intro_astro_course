import numpy as np
import blackbody as bb

wavelength = np.linspace(0.1, 2, 500) * 10**-6
temperature = [6000, 7000, 8000]
bb.bb_plot(wavelength, temperature)

