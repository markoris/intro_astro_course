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
    fac = np.exp(-fac / 2) / N
    fac /= np.max(fac)
    return fac

def make_image(n_pixels, n_objects=1):
	'''
	write me
	'''

	probabilities = np.random.uniform(size=n_objects)

	for i in range(n_objects):

		if probabilities[i] > 0.5: obj = 'galaxy'
		else: obj = 'star'

		scale_factor = 0.2
		if i > 1: scale_factor = 3

		mu = np.random.uniform(0, 1, size=2)
		if obj == 'galaxy':
			sigma = np.array([[np.random.uniform(0.2, 0.3), np.random.uniform(-0.2, 0.2)],
				 [np.random.uniform(-0.2, 0.2), np.random.uniform(0.2, 0.3)]])
		else:
			sigma = np.random.uniform(0.4, 0.8)
			sigma = np.array([[sigma, 0],
				 [0, sigma]])
			
		sigma /= (2*n_pixels)*scale_factor
		pos = image_grid(n_pixels)
		try:
			image += multivariate_gaussian(pos, mu, sigma)
		except NameError:
			image = multivariate_gaussian(pos, mu, sigma)

	return image

def plot_image(image):

	import matplotlib.pyplot as plt
	plt.imshow(image)
	plt.show()

def write_fits(image, filename):

	from astropy.io import fits
	hdu = fits.PrimaryHDU(image)
	hdul = fits.HDUList([hdu])
	try:
		hdul.writeto(filename)
	except IOError:
		import os
		os.remove(filename)
		hdul.writeto(filename)
		
