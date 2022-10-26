import numpy as np
import matplotlib.pyplot as plt

## constants                                                                                                                                                    
h = 6.626*10**-34 # [ J s]
c = 3.0*10**8 # [m s**-2]
k = 1.38*10**-23 # [J K]
sig = 5.67*10**-8 # [W m**-2 K**-4]
b = 2.897*10**-3 # [m K]

class Star:

    def __init__(self, spec_type):
        assert spec_type in 'OBAFGKM'
        self.spec_type = spec_type

    @property
    def temp(self):
        spec_class = {'O':41000, 'B':31000, 'A':9500, 'F':7240, 'G':5290, 'K':5300, 'M':3850}
        return spec_class[self.spec_type]
        
    ## Wien's law - turns temp to peak wavelength
    def wien(self):
        return b/self.temp
        
    ## Stefan-Boltazmann - turns temp to bb radiance
    def energy(self):
        return sig * self.temp**4

    ## Blackbody/Planck function
    def planck(self, wav):
        a = (2.0 * h * c**2) / wav**5
        b = h * c / (wav * k * float(self.temp))
        rad = a / (np.exp(b) - 1.0)
        return rad

## plotting function
def bb_plot(star, wav):
    if len(star) == 1:
        plt.plot(wav, star[0].planck(wav))
    else:
        for i in range(len(star)):
            plt.plot(wav, star[i].planck(wav))
    plt.xlabel('Wavelength')
    plt.ylabel('Intensity')
    plt.show()
