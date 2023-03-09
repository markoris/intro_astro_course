import numpy as np
import matplotlib.pyplot as plt
# -- hand-wavy approximation    

def A(t, t_merg):
    return np.where(t <= t_merg, t**1.4, np.exp(-10*(t-t_merg)))

dt = 0.01
t_merg = 1
f = 10
t = np.linspace(dt, 2, int(1e4))

h_t = A(t, t_merg)*np.cos(2*np.pi*f*t)

plt.figure(figsize=(19.2, 10.8))
plt.rc('font', size=30)
plt.rc('lines', lw=3)
plt.plot(t, h_t, c='k')
plt.savefig('wf.png')
