import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

G = 6.67 * 10**-8 #

gamma = 4.0 / 3.0
gamma_mo = 1.0 / 3.0
K = 4.897 * 10**14 #
M = 1.37 * 10**33 # g

rl = 1.45 * 10**6 # cm
ru = 1.65 * 10**6 # cm
rc = 1.55 * 10**6 # cm

rhol = 3.2 * 10**13 # g / cm**3

def analyticf(r, Omega, alpha):
	therml = rhol**( gamma_mo )
	gravl  = (G * M) / rl
	grav   = (G * M) / r
	rotl   = (rl * rl * Omega * Omega) * (rl / rc)**(2 * alpha) \
			 / ( 2 * alpha + 2.0 )
	rot    = (r * r * Omega * Omega  ) * (r / rc )**(2 * alpha) \
             / ( 2 * alpha + 2.0 )
	return therml + (gamma_mo ) * (grav - gravl + rot - rotl) \
                   / (K * gamma)

def analyticrho(r, Omega, alpha):
	return analyticf(r, Omega, alpha)**( 1.0 / ( gamma_mo ) )

def q(r, f, Omega, alpha):
	grav = - G * M / (r * r)
	rot = r * Omega * Omega * (r / rc)**(2 * alpha)
	return (1.0 / (gamma * K)) * gamma_mo * (grav + rot)

def solve(r, Omega, alpha):

	f0 = [rhol**( gamma_mo )]

	sol = solve_ivp(fun=q, t_span=(r[0],r[-1]), y0=f0, method='DOP853',
    	            t_eval=r, atol=1e-10, rtol=1e-13, args = (Omega, alpha))

	f = sol.y[0,:]

	rho = f**(1.0 / gamma_mo)

	return f, rho

r = np.linspace( rl, ru, 10000 )

rho_norot = solve(r, 0, 0)[1]
f_norot   = solve(r, 0, 0)[0]

rho  = solve(r, 1900, -1.25)[1]
f    = solve(r, 1900, -1.25)[0]
arho = analyticrho(r, 1900, -1.25)
af   = analyticf  (r, 1900, -1.25)

plt.figure(1)
plt.plot(r, rho,  '-',  linewidth=3)
plt.plot(r, arho, '--', linewidth=3)
plt.title(r'Newtonian ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = -1.25)$')
plt.xlabel(r'$\varpi$')
plt.ylabel(r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)')
plt.legend([ 'Numerical', 'Analytic'])

# How to plot with two different y-axes
# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/two_scales.html

plt.figure(2)
plt.plot(r / 10**6, rho / 10**13, 'k--', linewidth=3)
plt.plot(r / 10**6, 1900 * (r / rc)**-1.25 / 10**3,'k-.', linewidth=3)
plt.title(r'Newtonian ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = -1.25)$')
plt.xlabel('r ($km$)')
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.plot(r / 10**6, - G * M / r / 10**19, 'k-', linewidth=3)
ax1.legend([r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)', r'$\Omega$ ($10^{3}$ $Hz$)'], loc='upper center', handlelength=3)
ax2.legend([r'$\Phi$ ($10^{19}$ $erg$ $g^{-1}$)'], loc='lower center')

plt.figure(3)
plt.plot(r, (f - af) / rhol**( gamma_mo ), linewidth=3)
plt.plot(r, (rho - arho) / rhol, linewidth=3)
plt.show()
