import make_image as mi

star_image = mi.make_star_image(1024)
mi.plot_image(star_image)
mi.write_fits(star_image, "star.fits")
