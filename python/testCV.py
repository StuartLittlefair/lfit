from __future__ import absolute_import
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from trm import roche
import lfit
import time
import os


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
w = lfit.PyWhiteDwarf(rwd, 0.4)
d = lfit.PyDisc(q, rwd, rdisc, rexp, 1000)
# s = lfit.PySpot(q,rdisc,az,frac,scale)
s = lfit.PySpot(
    q, rdisc, az, frac, scale, exp1=exp1, exp2=exp2, tilt=tilt, yaw=yaw, complex=True
)
rs = lfit.PyDonor(q, 400)

ywd = w.calcFlux(q, inc, phi, width)
yd = d.calcFlux(q, inc, phi, width)
ys = s.calcFlux(q, inc, phi, width)
yrs = rs.calcFlux(q, inc, phi, width)
stop = time.time()
print("LFIT components took %f" % (stop - start))

"""
start = time.time()
os.system("../lfit_fake gfit.in 0.333 0.333 0.333 0.05 0.85 1.2 1000")
stop = time.time()
print("C++ version took %f" % (stop - start))
"""

pars = np.array(
    [
        0.333,
        0.333,
        0.333,
        0.05,
        q,
        dphi,
        rdisc,
        0.4,
        rwd,
        scale,
        az,
        frac,
        rexp,
        0.0,
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
flux = 0.333 * (ywd + yd + ys) + 0.05 * yrs

fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})
ax1 = ax[0]
ax2 = ax[1]
ax1.plot(phi, 0.333 * ywd, "--b")
ax1.plot(phi, 0.333 * yd, "--r")
ax1.plot(phi, 0.333 * ys, "-g")
ax1.plot(phi, 0.05 * yrs, "--y")
ax1.plot(phi, flux, "-g", label="Python")

ax2.plot(phi, (flux - flux2) / flux2)
# ax2.plot(phi-1,cv.ywd,label='Python')
# ax2.plot(phi-1,ywd,label='C++')
# ax2.legend()
plt.show()
