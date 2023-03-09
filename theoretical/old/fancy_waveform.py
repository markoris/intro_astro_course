# -- realistic amplitude

G = 4.3e-3 # pc (km/s)^2 / M_sun
c = 3e5 # km/s

def two_body_eom(pos_vel, time, m1):
    r_mag = np.linalg.norm(pos_vel[0:2])
    vel = pos_vel[2:]
    a = -m1*G*pos_vel[0:2] / (r_mag**3)
    return np.concatenate((vel, a))

def m_chirp(m1, m2):
    return (m1*m2)**(3/5)/(m1 + m2)**(1/5)

def h_0(m_c, D, f):
    return 4*G/c**2*m_c/D*(G/c**3*np.pi*f*m_c)**(2/3)

def chirp(f, m_c):
    return 96/5*c**3/G*f/m_c*(G/c**3*np.pi*f*m_c)**(8/3)

dt = 0.001
t = np.linspace(dt, 0.005, int(1e4))

m1 = 10
m2 = 5
D = 1e6 # 1 Mpc
f = np.zeros_like(t)
pos_vel = np.zeros(len(t), 4)
pos_vel[0, 0] = 1e5 # x-position (km)
pos_vel[0, 1] = 0 # y-position (km)
pos_vel[0, 2] = 0 # x-velocity (km/s)
pos_vel[0, 3] = 1e3 # y-velocity (km/s)

for idx in range(len(t)):
# need to add some kind of energy-loss term to the two_body_EOM
# this will cause a decaying orbit for which we calculate the radius at each time step
# from the radius/velocity we get a new orbital period at each time
# the orbital period implies an orbital frequency
# the orbital frequency relates to the GW frequency as f = 2*f_orb
# we then have the frequency vector evolve as a function of time, should create increasingly closer-spaced waveform        


f = 1/t # Hz
print(f)
m_c = m_chirp(m1, m2)
h0 = h_0(m_c, D, f)
f_dot = chirp(f, m_c)

h_t = h0*np.cos(2*np.pi*f*t + np.pi*f_dot*t**2)

plt.figure(figsize=(19.2, 10.8))
plt.rc('font', size=30)
plt.rc('lines', lw=3)
plt.plot(t, h_t, c='k')
