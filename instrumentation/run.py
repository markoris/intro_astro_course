import image as img

image = img.Image(1024)
image.add_objects(25)
image.add_effects(["telescope_cover", "dark_current"])
image.plot()
image.write_fits("my_image.fits")
