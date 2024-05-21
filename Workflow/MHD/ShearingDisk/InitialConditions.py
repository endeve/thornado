import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import h5py

G = 6.67 * 10**-8 #

gamma = 4.0 / 3.0
gamma_mo = 1.0 / 3.0
K = 4.897 * 10**14 # cgs
M = 1.37 * 10**33 # g
c =  2.99792458 * 10**10 # cm / s
rl = 1.45 * 10**6 # cm
ru = 1.65 * 10**6 # cm
rc = 1.55 * 10**6 # cm

rhol = 3.2 * 10**13 # g / cm**3

def psi(r):

	return (1.0 + X(r))

def dpsi(r):

	return -X(r) / r

def alpha(r):

	return (1.0 - X(r)) / (1.0 + X(r))

def dlnalpha(r):

	return (2.0 * X(r) / r) * (1.0 / (1.0 - X(r)**2))

def gamma33(r):

	return psi(r)**4 * r**2

def dgamma33(r):

	return 2.0 * r * psi(r)**4 + 4.0 * r**2 * psi(r)**3 * dpsi(r)

def Omega(r, Omega0, alphaO):

	return Omega0 * (r / rc)**alphaO

def W(r, Omega0, alphaO):

	return 1.0 / np.sqrt(1.0 - gamma33(r) * Omega(r, Omega0, alphaO)**2 / c**2)

def X(r):

	return (G * M) / (2.0 * r * c**2)

def Delta(r, Omega0):

	print(r)
	print(Omega0)
	print(r**2 * Omega0**2 * (1.0 + X(r))**6)
	print(c**2 * (1.0 - X(r))**2)
	return np.sqrt(-r**2 * Omega0**2 * (1.0 + X(r))**6 + c**2 * (1.0 - X(r))**2)

def analyticf_N(r, Omega0, alphaO):

	therml = rhol**gamma_mo
	gravl  = (G * M) / rl
	grav   = (G * M) / r
	rotl   = (rl**2 * Omega(rl, Omega0, alphaO)**2) \
			 / (2.0 * alphaO + 2.0)
	rot    = (r**2  * Omega(r, Omega0, alphaO)**2)   \
             / (2.0 * alphaO + 2.0)
	return therml + gamma_mo * (grav - gravl + rot - rotl) \
                   / (K * gamma)

def analyticf_GR_norot(r):

	const_GR = (1.0 / (1.0 + X(rl))) \
               * (rhol**gamma_mo * (1.0 - X(rl)) - (G * M / rl) * gamma_mo / (K * gamma))
	grav     = const_GR * (1.0 + X(r)) / (1.0 - X(r))
	therml   = (G * M / r) * (gamma_mo / (K * gamma)) * (1.0 / (1.0 - X(r)))
	return grav + therml

def analyticf_GR_constrot(r, Omega0):

	const_GR = 4.0 * c**4 * (rhol**gamma_mo + gamma_mo * c**2 / (K * gamma)) \
               * (Delta(rl, Omega0) / (1.0 + X(rl)))
	grav     = const_GR * (1.0 + X(r)) / (4.0 * c**4 * Delta(r, Omega0))
	therml   = gamma_mo * c**2 / (K * gamma)
	return grav - therml

def analyticrho_N(r, Omega0, alphaO):

	return analyticf_N(r, Omega0, alphaO)**(1.0 / gamma_mo)

def analyticrho_GR_norot(r):

	return analyticf_GR_norot(r)**(1.0 / gamma_mo)

def analyticrho_GR_constrot(r, Omega0):

	return analyticf_GR_constrot(r, Omega0)**(1.0 / gamma_mo)

def q_N(r, Omega0, alphaO):

	grav = - G * M / r**2
	rot = r * Omega(r, Omega0, alphaO)**2
	return (1.0 / (gamma * K)) * (grav + rot)

def diffEq_N(r, f, Omega0, alphaO):

	return gamma_mo * q_N(r, Omega0, alphaO)

def q_GR(r, Omega0, alphaO):

	grav = dlnalpha(r)
	rot = Omega(r, Omega0, alphaO)**2 \
          * dgamma33(r) / (2.0 * alpha(r)**2 * c**2)
	return (W(r, Omega0, alphaO)**2 / gamma_mo) * (grav - rot)

def r_GR(r, Omega0, alphaO):

	grav = dlnalpha(r)
	rot = Omega(r, Omega0, alphaO)**2 \
          * dgamma33(r) / (2.0 * alpha(r)**2 * c**2)
	return (c**2 * W(r, Omega0, alphaO)**2 / (K * gamma)) * (rot - grav)

def diffEq_GR(r, f, Omega0, alphaO):

	return gamma_mo * r_GR(r, Omega0, alphaO) \
           - gamma_mo * q_GR(r, Omega0, alphaO) * f

def solve_N(r, Omega0, alphaO):

	f0 = [rhol**( gamma_mo )]

	sol = solve_ivp(fun = diffEq_N, t_span = (r[0],r[-1]), y0 = f0, method = 'DOP853',
    	            t_eval = r, atol = 1e-13, rtol = 1e-13, args = (Omega0, alphaO))

	f = sol.y[0,:]

	rho = f**(1.0 / gamma_mo)

	return f, rho

def solve_GR(r, Omega0, alphaO):

	f0 = [rhol**( gamma_mo )]

	sol = solve_ivp(fun = diffEq_GR, t_span = (r[0],r[-1]), y0 = f0, method = 'DOP853',
    	            t_eval = r, atol = 1e-13, rtol = 1e-13, args = (Omega0, alphaO))

	f = sol.y[0,:]

	rho = f**(1.0 / gamma_mo)

	return f, rho

def output(r, rho, filename):

	file = h5py.File(filename + '.hdf5', 'w')

	r_dset   = file.create_dataset("r",        r.shape,   dtype = 'd')
	rho_dset = file.create_dataset("rho_dset", rho.shape, dtype = 'd')

	r_dset[:] = r
	rho_dset[:] = rho

r = np.linspace( rl, ru, 10000 )

arho_N_norot_NRparams     = analyticrho_N(r, 0, 0)
arho_N_constrot_NRparams  = analyticrho_N(r, 1900, 0)
arho_N_diffrot_NRparams   = analyticrho_N(r, 1900, -1.25)
arho_GR_norot_NRparams    = analyticrho_GR_norot(r)
arho_GR_constrot_NRparams = analyticrho_GR_constrot(r, 1900)

rho_N_norot_NRparams     = solve_N(r, 0, 0)[1]
rho_N_constrot_NRparams  = solve_N(r, 1900, 0)[1]
rho_N_diffrot_NRparams   = solve_N(r, 1900, -1.25)[1]
rho_GR_norot_NRparams    = solve_GR(r, 0, 0)[1]
rho_GR_constrot_NRparams = solve_GR(r, 1900, 0)[1]
rho_GR_diffrot_NRparams  = solve_GR(r, 1900, -1.25)[1]

output(r, rho_GR_norot_NRparams,    'GR_norot_NRparams'   )
output(r, rho_GR_constrot_NRparams, 'GR_constrot_NRparams')
output(r, rho_GR_diffrot_NRparams,  'GR_diffrot_NRparams' )

plt.figure(0)
plt.plot(r, rho_N_diffrot_NRparams,  '-',  linewidth=3)
plt.plot(r, arho_N_diffrot_NRparams, '--', linewidth=3)
plt.title(r'Newtonian ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = -1.25)$')
plt.xlabel(r'$\varpi$')
plt.ylabel(r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)')
plt.legend([ 'Numerical', 'Analytic'])

plt.figure(10)
plt.plot(r, rho_GR_norot_NRparams,  '-',  linewidth=3)
plt.plot(r, arho_GR_norot_NRparams, '--', linewidth=3)
plt.title(r'GR ($\Omega_{0} = 0 s^{-1}$, $\alpha_{\Omega} = 0)$')
plt.xlabel(r'$\varpi$')
plt.ylabel(r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)')
plt.legend([ 'Numerical', 'Analytic'])

plt.figure(20)
plt.plot(r, rho_GR_constrot_NRparams,  '-',  linewidth=3)
plt.plot(r, arho_GR_constrot_NRparams, '--', linewidth=3)
plt.title(r'GR ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = 0)$')
plt.xlabel(r'$\varpi$')
plt.ylabel(r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)')
plt.legend([ 'Numerical', 'Analytic'])

# How to plot with two different y-axes
# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/two_scales.html

plt.figure(2)
plt.plot(r / 10**5, rho_N_diffrot_NRparams / 10**13, 'k--', linewidth=3)
plt.plot(r / 10**5, 1900 * (r / rc)**-1.25 / 10**3,'k-.', linewidth=3)
plt.title(r'Newtonian ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = -1.25)$')
plt.xlabel('r ($km$)')
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.plot(r / 10**5, - G * M / r / 10**19, 'k-', linewidth=3)
ax1.legend([r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)', r'$\Omega$ ($10^{3}$ $Hz$)'],
           loc='upper center', handlelength=3)
ax2.legend([r'$\Phi$ ($10^{19}$ $erg$ $g^{-1}$)'], loc='lower center')
ax1.set_ylabel(r'$\rho$ and $\Omega$')
ax2.set_ylabel(r'$\Phi$')

plt.figure(30)
plt.plot(r / 10**5, arho_N_norot_NRparams,     linewidth=3, color='b')
plt.plot(r / 10**5, arho_N_constrot_NRparams,  linewidth=3, color='y')
plt.plot(r / 10**5, arho_N_diffrot_NRparams,   linewidth=3, color='r')
plt.plot(r / 10**5, rho_GR_norot_NRparams,     linewidth=3, color='b',  linestyle = '--')
plt.plot(r / 10**5, rho_GR_constrot_NRparams,  linewidth=3, color='y',  linestyle = '--')
plt.plot(r / 10**5, rho_GR_diffrot_NRparams,   linewidth=3, color='r',  linestyle = '--')
plt.xlabel('r ($km$)')
plt.ylabel(r'$\rho$ ($10^{13}$ $g$ $cm^{-3}$)')
plt.legend([ 'Analytic Newtonian (No Rotation)',
             r'Analytic Newtonian ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} =  0.00$)',
             r'Analytic Newtonian ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = -1.25$)',
             'Numerical GR (No Rotation)',
             r'Numerical GR ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} =  0.00$)',
             r'Numerical GR ($\Omega_{0} = 1900 s^{-1}$, $\alpha_{\Omega} = -1.25$)'])

plt.figure(1)
plt.plot(r, (rho_N_diffrot_NRparams - arho_N_diffrot_NRparams) / rhol, linewidth=3)

plt.figure(11)
plt.plot(r, (rho_GR_norot_NRparams - arho_GR_norot_NRparams) / arho_GR_norot_NRparams, linewidth=3)

plt.figure(21)
plt.plot(r, (rho_GR_constrot_NRparams - arho_GR_constrot_NRparams) / arho_GR_constrot_NRparams, linewidth=3)

plt.figure(31)
plt.plot(r / 10**5, psi(r),   linewidth=3)
plt.plot(r / 10**5, alpha(r), linewidth=3)
plt.plot(r / 10**5, 1 + (gamma / gamma_mo) * K * arho_GR_norot_NRparams**gamma_mo / c**2, linewidth=3)
plt.xlabel('r ($km$)')
plt.title('GR (No Rotation)')
plt.legend([r'$\psi$', r'$\alpha$', r'h ($c^{2}$)'])

#plt.show()
