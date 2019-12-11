import lfit
import numpy as np
import matplotlib.pyplot as plt
import time

parlist = [0.05246,
 0.01703,
 1.0, # 0.02051,
 0.0061,
 0.05852,
 0.03861,
 0.55342,
 0.284,
 0.0261,
 0.06986,
 117.24748,
 0.36464,
 1.66211,
 0.00014]

N_data = 100
N_trials = 500


# phi = np.linspace(0.05971, 0.059716, N_data)
phi = np.linspace(-0.05, 0.2, N_data)

t_calc = []
for i in range(N_trials):
    t0 = time.perf_counter()
    cv = lfit.CV(parlist)
    flx = cv.calcFlux(parlist, phi)
    t1 = time.perf_counter()
    print("done {:>04d}".format(i), end='\r')
    t_calc.append(t1 - t0)


av_t = np.mean(t_calc)*1000.
std_t = np.std(t_calc)*1000.
print("Flux calculation, across {:d} data, took {:.3}+/-{:.3}ms (trialed {:d} times)".format(N_data, av_t, std_t, N_trials))

plt.scatter(phi, flx, marker='x')
plt.xlim((phi.min(), phi.max()))
plt.ylim((flx.min(), flx.max()))
plt.title("Bright Spot flux, only")
plt.show()
