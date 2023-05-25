import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import block_reduce

class Image():

    def __init__(self, n_pixels):

        self.n_pixels = n_pixels

        # load the original image
        image_orig = plt.imread('jwst_deep_field.png')
        # get original image size
        shape_orig = image_orig.shape[0]
        # calculate downsampling window size given user-input image size
        bs = (int(shape_orig/self.n_pixels), int(shape_orig/self.n_pixels), 1)
        # generate downsampled image
        self.image = block_reduce(image_orig, block_size=bs, func=np.mean)

        self.eff_dict = {#"dark_current": self.dark_current,
                 "satellite_transits": self.satellite_transits,
                 "dead_pixels": self.dead_pixels,
                 "dead_arrays": self.dead_arrays,
                 "telescope_cover": self.telescope_cover}

#    def __init__(self, n_pixels):
#
#        self.n_pixels = n_pixels
#
#        x = np.linspace(0, 1, n_pixels, endpoint=True)
#        y = np.linspace(0, 1, n_pixels, endpoint=True)
#        X, Y = np.meshgrid(x, y)
#        self.pos = np.empty(X.shape + (2,))
#        self.pos[:, :, 0] = X
#        self.pos[:, :, 1] = Y
#
#        self.eff_dict = {"dark_current": self.dark_current,
#                 "dead_pixels": self.dead_pixels,
#                 "dead_arrays": self.dead_arrays,
#                 "telescope_cover": self.telescope_cover}

#    def add_objects(self, n_objects=20, prob=0.3):
#        '''
#        write me
#        '''
#
#        probabilities = np.random.uniform(size=n_objects)
#        mus = np.random.uniform(0.05, 0.95, size=(n_objects,2))
#        r = np.array([np.sqrt(np.sum((mu1-mu2)**2)) for mu1 in mus for mu2 in mus])
#        r = np.unique(r[np.where(r > 0)[0]])
#        criterion = 0.1-((n_objects-10)*0.0025)
#        while np.any(r < criterion):
#            mus = np.random.uniform(0.05, 0.95, size=(n_objects,2))
#            r = np.array([np.sqrt((mu1[0]-mu2[0])**2 + (mu1[1]-mu2[1])**2) for mu1 in mus for mu2 in mus])
#            r = np.unique(r[np.where(r > 0)[0]])
#
#        for i in range(n_objects):
#
#            if probabilities[i] < prob: obj = 'galaxy'
#            else: obj = 'star'
#
#            scale_factor = 0.5
#            if i > 1: scale_factor = 3
#
#            mu = mus[i]
#            if obj == 'galaxy':
#                sigma = np.array([[np.random.uniform(0.2, 0.3), np.random.uniform(-0.2, 0.2)],
#                     [np.random.uniform(-0.2, 0.2), np.random.uniform(0.2, 0.3)]])
#            else:
#                sigma = np.random.uniform(0.4, 0.8)
#                sigma = np.array([[sigma, 0],
#                     [0, sigma]])
#                
#            sigma /= (2*self.n_pixels)*scale_factor
#            
#            try:
#                self.image += self.multivariate_gaussian(self.pos, mu, sigma)
#            except:
#                self.image = self.multivariate_gaussian(self.pos, mu, sigma)

    def add_effects(self, effects):
    
        #if "telescope_cover" in effects and "dark_current" not in effects:
        if "telescope_cover" in effects and "satellite_transits" not in effects:
            effects.remove("telescope_cover")
            effects.insert(0, "telescope_cover")

        #if "dark_current" in effects and "telescope_cover" not in effects:
        if "satellite_transits" in effects and "telescope_cover" not in effects:
            #effects.remove("dark_current")
            #effects.insert(0, "dark_current")
            effects.remove("satellite_transits")
            effects.insert(0, "satellite_transits")

        #if "telescope_cover" in effects and "dark_current" in effects:
        if "telescope_cover" in effects and "satellite_transits" in effects:
            effects.remove("telescope_cover")
            effects.remove("satellite_transits")
            effects.insert(0, "telescope_cover")
            # do not reinsert satellite transits, as they should not appear!
            #effects.remove("dark_current")
            #effects.insert(0, "dark_current")
    
        for effect in effects:
            self.eff_dict[effect]()

        if np.max(self.image) > 0.5: # if telescope cover is on, we want extremely low photon counts
            self.image /= np.max(self.image)

#    def dark_current(self):
#
#        self.image += np.random.uniform(0.1, 0.2, size=(self.n_pixels, self.n_pixels))

    def satellite_transits(self, n_transits=5):
        from skimage.draw import line
        transits = np.zeros((self.n_pixels, self.n_pixels, 3))
        transit_paths = np.random.uniform(0, self.n_pixels-1, size=(n_transits, 2)).astype('int')
        print(transit_paths)
        for i in range(n_transits):
            rr, cc = line(transit_paths[i, 0], 0, transit_paths[i, 1], self.n_pixels-1)
            self.image[rr, cc, :] = 1

    def dead_pixels(self):

        dead_pixels = np.random.randint(0, self.n_pixels-1, size=(10000, 2)) # ~20% of total pixels?
        for i in range(dead_pixels.shape[0]):
            self.image[dead_pixels[i, 0], dead_pixels[i, 1]] = 0

    def dead_arrays(self):

        dead_arrays = np.random.randint(0, self.n_pixels-1, size=(30))
        self.image[dead_arrays, :] = 0
    
    def telescope_cover(self):
        self.image = np.zeros((self.n_pixels, self.n_pixels), dtype=np.int8)

    def plot(self):

        import matplotlib.pyplot as plt
        
        # Gaussians representing stars/galaxies will have floating point noise far away from objects
        # Anything under a value of 0.01 should just be floored to 0 for highlighting effects of dark current
        
        #self.image = np.where(((self.image < 0.01) & (self.image != 0)), 0, self.image)
        plt.figure(figsize=(7, 7))
        plt.imshow(self.image)
        plt.tight_layout()

    def write_fits(self, filename):

         # placeholder if using FITS files in future
         # currently returns nothing
         return
        
#        from astropy.io import fits
#        hdu = fits.PrimaryHDU(self.image)
#        hdul = fits.HDUList([hdu])
#        try:
#            hdul.writeto(filename)
#        except IOError:
#            import os
#            os.remove(filename)
#            hdul.writeto(filename)
        
#    def multivariate_gaussian(self, pos, mu, Sigma):
#        """
#        Return the multivariate Gaussian distribution on array pos.
#
#        pos is an array constructed by packing the meshed arrays of variables
#        x_1, x_2, x_3, ..., x_k into its _last_ dimension.
#
#        """
#        n = mu.shape[0]
#        Sigma_det = np.linalg.det(Sigma)
#        Sigma_inv = np.linalg.inv(Sigma)
#        N = np.sqrt((2*np.pi)**n * Sigma_det)
#        # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
#        # way across all the input variables.
#        fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)
#        fac = np.exp(-fac / 2) / N
#        fac /= np.max(fac)
#        return fac

