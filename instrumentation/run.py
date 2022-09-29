import make_image as mi

image = mi.make_image(1024, n_objects=25)
mi.plot_image(image)
mi.write_fits(image, 'image.fits')
