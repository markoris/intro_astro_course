import numpy as np

def image_grid(n_pixels):

	x = np.linspace(0, 1, n_pixels, endpoint=True)
	y = np.linspace(0, 1, n_pixels, endpoint=True)
	X, Y = np.meshgrid(x, y)
	pos = np.empty(X.shape + (2,))
	pos[:, :, 0] = X
	pos[:, :, 1] = Y
	return pos

def multivariate_gaussian(pos, mu, Sigma):
    """Return the multivariate Gaussian distribution on array pos.

    pos is an array constructed by packing the meshed arrays of variables
    x_1, x_2, x_3, ..., x_k into its _last_ dimension.

    """
    n = mu.shape[0]
    Sigma_det = np.linalg.det(Sigma)
    Sigma_inv = np.linalg.inv(Sigma)
    N = np.sqrt((2*np.pi)**n * Sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)

    return np.exp(-fac / 2) / N

def make_star_image(n_pixels, multiple=False):
	
	'''
	We want 0 covariance, same-sigma matrix for stars since they are, indeed, spherical
	'''

	mu = np.random.uniform(0.2, 0.8, size=2)
	sigma = np.random.uniform(1, 5)/(2*n_pixels)
	sigma = np.array([[sigma, 0],
			 [0, sigma]])
	pos = image_grid(n_pixels)
	star_image = multivariate_gaussian(pos, mu, sigma)
	return star_image
		

def make_galaxy_image(n_pixels, multiple=False):
	'''
	write me
	'''
	galaxy_image = 0
	#sigma = np.array([[np.random.uniform(0.1, 0.4), np.random.uniform(0, 0.2)],
	#		 [np.random.uniform(0, 0.2,), np.random.uniform(0.1, 0.4)]])
	#sigma = np.array([[np.random.uniform(0.01, 0.1), 0],
	#		 [0, np.random.uniform(0.01, 0.1)]])
	return galaxy_image

def plot_image(image):
	import matplotlib.pyplot as plt
	plt.imshow(image)
	plt.show()

def write_fits(image, filename):
	from astropy.io import fits
	hdu = fits.PrimaryHDU(image)
	hdul = fits.HDUList([hdu])
	hdul.writeto(filename)
