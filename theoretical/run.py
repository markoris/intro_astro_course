import numpy as np
import star

wavelength = np.linspace(100, 2000, 500) * 10**-9 # nanometers

star1 = star.Star('K')
star2 = star.Star('A')
star.bb_plot([star1,star2], wavelength)
