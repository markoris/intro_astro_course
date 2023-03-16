import lightcurve as lc

# Kepler-695b: 5651104
known_period_565 = 3.040330459

data = lc.load_data('5651104/plot.tbl')
time = data[0]
flux = data[1]
lc.plot_lc(time[0:500],flux[0:500])

t_fold = lc.fold(time[0:500], known_period_565)
lc.plot_phase(t_fold, flux[0:500])
exit(0)
# Kepler-18c: 8644288
known_period_864 = 7.641567589 #14.8589085

data = lc.load_data('8644288/plot.tbl')
time = data[0]
flux = data[1]
lc.plot_lc(time[500:1000],flux[500:1000])

t_fold = lc.fold(time[0:5000], known_period_565)
lc.plot_phase(t_fold, flux[0:5000])
