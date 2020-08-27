import numpy as np

# Define units, with GeV as base unit
GeV = 10**6
eV = 10**-9*GeV
KeV = 10**-6*GeV
MeV = 10**-3*GeV
TeV = 10**3*GeV

Sec = (1/(6.582119*10**-16))/eV
Kmps = 3.3356*10**-6
Centimeter = 5.067730716156396*10**13/GeV
Meter = 100*Centimeter
Km = 10**5*Centimeter
Kilogram = 5.6095883571872*10**35*eV
Gram = 1e-3 * Kilogram
Day = 86400*Sec
Year = 365*Day
KgDay = Kilogram*Day
amu = 1.66053892*10**-27*Kilogram
Mpc = 3.086*10**24*Centimeter
joule = Kilogram*Meter**2/Sec**2
erg = 1e-7*joule
Angstrom = 1e-10*Meter
Kelv = 8.62e-14 * GeV
k_B = 1.3806488e-16 * erg / Kelv
barn = 1e-24 * Centimeter ** 2
Jy = 1e-23 * erg / Sec / Centimeter ** 2 / Sec ** -1
Hz = 1 / Sec

# Electromagnetic units
# Mx = Centimeter ** 1.5 * Gram ** 0.5 * Sec ** -1
# Gauss = Mx / Centimeter ** 2
# nGauss = 1e-9 * Gauss
alpha = 1 / 137
q_e = np.sqrt(4 * np.pi * alpha)
Coulomb = 6.2415090741e18 * q_e
Amp = Coulomb / Sec
Wb = Kilogram * Meter ** 2 / Sec ** 2 / Amp
Mx = 1e-8 * Wb
Gauss = Mx / Centimeter ** 2
nGauss = 1e-9 * Gauss

# Particle and astrophysics parameters
M_s = 1.99*10**30*(Kilogram)

# Some conversions
kpc = 1.e-3*Mpc
pc = 1e-3*kpc
asctorad = np.pi/648000.
radtoasc = 648000./np.pi

# Constants
GN = 6.67e-11*Meter**3/Kilogram/Sec**2
h = 0.7
H_0 = 100*h*(Kmps/Mpc)
rho_c = 3*H_0**2/(8*np.pi*GN)
Rsun = 8.*kpc
R200_MW = 200.*kpc
m_p = 1.67261777e-27 * Kilogram  # Proton mass
m_e = 0.5109989461 * MeV  # Electron mass

# To get velocity of Earth
vE0 = 29.79
omega = 2*np.pi/365.25
e = 0.016722
e1 = np.array([.9940, .1095, .003116])
e2 = np.array([-.05173, .4945, -.8677])
e1 /= np.linalg.norm(e1)
e2 /= np.linalg.norm(e2)
lambdap = np.deg2rad(281.93)
t1 = 79.3266
