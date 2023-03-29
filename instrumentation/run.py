import image as img

image = img.Image(1024)
#image = img.Image(256)
image.add_objects(40, prob=0.3)
#image.add_effects([])
#image.add_effects(["dark_current"])
#image.add_effects(["telescope_cover"])
#image.add_effects(["dead_pixels"])
#image.add_effects(["dead_arrays"])
#image.add_effects(["dark_current", "telescope_cover"])
image.plot()
image.write_fits("my_image.fits")
