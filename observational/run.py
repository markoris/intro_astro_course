import lightcurve as lc

known_period_108 = 2.5255917769                                                                                              
known_period_105 = 9.27358173
known_period_565 = 3.040330459

data = lc.load_data('5651104/plot.tbl')
time = data[0]
flux = data[1]
lc.plot_lc(time[0:1500],flux[0:1500])

t_fold = lc.fold(time[0:5000], known_period_565)
lc.plot_phase(t_fold, flux[0:5000])
