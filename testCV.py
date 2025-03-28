import time

import lfit
import matplotlib.pyplot as plt
import numpy as np
from trm import roche

phi = np.linspace(0.85, 1.2, 1000, endpoint=False)
width = np.mean(np.diff(phi)) * np.ones_like(phi) / 2.0

q = 0.1
inc = 85.9
xl1 = roche.xl1(q)
dphi = roche.findphi(q, inc)
rwd = 0.01
rdisc = 0.6
rexp = 0.2
az = 157.0
frac = 0.2
scale = 0.039
exp1 = 2.0
exp2 = 1.0
tilt = 120.0
yaw = 1.0


start = time.time()
# calculate individual components
w = lfit.PyWhiteDwarf(rwd, 0.4)
d = lfit.PyDisc(q, rwd, rdisc, rexp, 1000)

# simple bright spot
s = lfit.PySpot(q, rdisc, az, frac, scale, exp1)
s = lfit.PySpot(
    q, rdisc, az, frac, scale, exp1=exp1, exp2=exp2, tilt=tilt, yaw=yaw, complex=True
)
rs = lfit.PyDonor(q, 400)

# each component has a calcFlux method
ywd = w.calcFlux(q, inc, phi, width)
yd = d.calcFlux(q, inc, phi, width)
ys = s.calcFlux(q, inc, phi, width)
yrs = rs.calcFlux(q, inc, phi, width)
stop = time.time()
print("LFIT components took %f" % (stop - start))

# The usual way to do this is to use the CV class
# which will calculate the fluxes for all components
# and sum them up
pars = np.array(
    [
        0.333,  # flux fraction of white dwarf
        0.333,  # flux fraction of disc
        0.333,  # flux fraction of bright spot
        0.05,  # flux fraction of donor
        q,
        dphi,
        rdisc,
        0.4,  # limb darkening of white dwarf
        rwd,
        scale,
        az,
        frac,
        rexp,
        0.0,  # phase offset
        exp1,
        exp2,
        tilt,
        yaw,
    ]
)

start = time.time()
cv = lfit.CV(pars)
flux2 = cv.calcFlux(pars, phi)
stop = time.time()
print("LFIT CV took %f" % (stop - start))
# flux from previous component calcs
flux = 0.333 * (ywd + yd + ys) + 0.05 * yrs

fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})
ax1 = ax[0]
ax2 = ax[1]
ax1.plot(phi, 0.333 * ywd, "--b")
ax1.plot(phi, 0.333 * yd, "--r")
ax1.plot(phi, 0.333 * ys, "-g")
ax1.plot(phi, 0.05 * yrs, "--y")
ax1.plot(phi, flux2, "-g", label="components")
ax2.plot(phi, (flux - flux2) / flux2)
plt.show()
