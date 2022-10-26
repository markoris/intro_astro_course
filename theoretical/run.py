import numpy as np
import star

wavelength = np.linspace(1, 50, 500) * 10**-7

star1 = star.Star('K')
star2 = star.Star('A')
star.bb_plot([star1,star2], wavelength)
