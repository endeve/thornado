import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import h5py
from numpy import pi

G = 6.674 * 10**-8 #
c =  2.99792458 * 10**10 # cm / s

gamma = 4.0 / 3.0
gamma_mo = 1.0 / 3.0
SolarMass = 1.989e+33

Params =  [1.44 * 10**6, 1.66 * 10**6, 4.897 * 10**14, 0.68, 3.2 * 10**13, 1900.0, -1.25] #Less Relativistic
#Params = [1.44 * 10**6, 1.66 * 10**6, 9.1 * 10**14, 2.70, 3.2 * 10**13, 0, 0] # 1739 * pi, -1.33] #More Relativistic

Filename = 'GR_LR_diffrot'

rl = Params[0]
ru = Params[1]
rc = (rl + ru) / 2

K = Params[2]
M = Params[3] * SolarMass
rhol = Params[4]
Omega0 = Params[5]
alphaO = Params[6]

def X(r):

	return (G * M) / (2.0 * r * c**2)

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

def Omega(r):

	return Omega0 * (r / rc)**alphaO

def W(r):

	return 1.0 / np.sqrt(1.0 - gamma33(r) * Omega(r)**2 / (alpha(r)**2 * c**2))

def Delta(r):

	return np.sqrt(-r**2 * Omega0**2 * (1.0 + X(r))**6 + c**2 * (1.0 - X(r))**2)

def p(r, rho):

	return K * rho**gamma

def analyticf_N(r):

	therml = rhol**gamma_mo
	gravl  = (G * M) / rl
	grav   = (G * M) / r
	rotl   = (rl**2 * Omega(rl)**2) \
			 / (2.0 * alphaO + 2.0)
	rot    = (r**2  * Omega(r)**2)   \
             / (2.0 * alphaO + 2.0)
	print(rotl)
	print(rot)
	return therml + (gamma_mo / (K * gamma)) * (grav - gravl + rot - rotl)

def analyticf_GR_norot(r):

	const_GR = (1.0 / (1.0 + X(rl))) \
               * (rhol**gamma_mo * (1.0 - X(rl)) - (G * M / rl) * gamma_mo / (K * gamma))
	grav     = const_GR * (1.0 + X(r)) / (1.0 - X(r))
	therml   = (G * M / r) * (gamma_mo / (K * gamma)) * (1.0 / (1.0 - X(r)))
	return grav + therml

def analyticf_GR_constrot(r):

	const_GR = 4.0 * c**4 * (rhol**gamma_mo + gamma_mo * c**2 / (K * gamma)) \
               * (Delta(rl) / (1.0 + X(rl)))
	grav     = const_GR * (1.0 + X(r)) / (4.0 * c**4 * Delta(r))
	therml   = gamma_mo * c**2 / (K * gamma)
	return grav - therml

def analyticrho_N(r):

	return analyticf_N(r)**(1.0 / gamma_mo)

def analyticrho_GR_norot(r):

	return analyticf_GR_norot(r)**(1.0 / gamma_mo)

def analyticrho_GR_constrot(r):

	return analyticf_GR_constrot(r)**(1.0 / gamma_mo)

def q_N(r):

	grav = - G * M / r**2
	rot = r * Omega(r)**2
	return (1.0 / (gamma * K)) * (grav + rot)

def diffEq_N(r, f):

	return gamma_mo * q_N(r)

def q_GR(r):

	grav = dlnalpha(r)
	rot = Omega(r)**2 \
          * dgamma33(r) / (2.0 * alpha(r)**2 * c**2)
	return (W(r)**2 / gamma_mo) * (grav - rot)

def r_GR(r):

	grav = dlnalpha(r)
	rot = Omega(r)**2 \
          * dgamma33(r) / (2.0 * alpha(r)**2 * c**2)
	return (c**2 * W(r)**2 / (K * gamma)) * (rot - grav)

def diffEq_GR(r, f):

	return gamma_mo * r_GR(r) \
           - gamma_mo * q_GR(r) * f

def solve_N(r):

	f0 = [rhol**( gamma_mo )]

	sol = solve_ivp(fun = diffEq_N, t_span = (r[0],r[-1]), y0 = f0, method = 'DOP853',
    	            t_eval = r, atol = 1e-13, rtol = 1e-13)

	f = sol.y[0,:]

	rho = f**(1.0 / gamma_mo)

	return f, rho

def solve_GR(r):

	f0 = [rhol**( gamma_mo )]

	sol = solve_ivp(fun = diffEq_GR, t_span = (r[0],r[-1]), y0 = f0, method = 'DOP853',
    	            t_eval = r, atol = 1e-13, rtol = 1e-13)

	f = sol.y[0,:]

	rho = f**(1.0 / gamma_mo)

	return f, rho

def output(r, rho, filename):

	file = h5py.File(filename + '.h5', 'w')

	r_dset   = file.create_dataset("r",    r.shape, dtype = 'd')
	rho_dset = file.create_dataset("rho",  r.shape, dtype = 'd')
	p_dset   = file.create_dataset("pres", r.shape, dtype = 'd')
	V3_dset  = file.create_dataset("V3",   r.shape,  dtype = 'd')

	alpha_dset = file.create_dataset("alpha", r.shape, dtype = 'd')
	psi_dset   = file.create_dataset("psi",   r.shape, dtype = 'd')

	r_dset[:] = r
	rho_dset[:] = rho
	p_dset[:] = p(r, rho)
	V3_dset[:] = Omega(r) / alpha(r)

	alpha_dset[:] = alpha(r)
	psi_dset[:] = psi(r)

	file.close()

r = np.linspace( rl, ru, 10000 )

arho_N = analyticrho_N(r)

if (alphaO == 0.0 and Omega0 == 0.0):
	arho_GR = analyticrho_GR_norot(r)
	print('No Rotation')
elif(alphaO == 0.0 and Omega0 != 0.0):
	arho_GR = analyticrho_GR_constrot(r)
	print('Const. Rotation')
else:
	arho_GR = np.zeros(np.size(r))

rho_N  = solve_N(r)[1]
rho_GR = solve_GR(r)[1]

print(rho_N)
print(rho_GR)

plt.figure(1)
plt.plot(rho_N)
plt.plot(arho_N)

plt.figure(2)
plt.plot(abs(rho_N  - arho_N ) / arho_N )
#plt.plot(abs(rho_GR - arho_GR) / arho_GR )
plt.show()
print(max(abs(rho_N  - arho_N ) / arho_N ))
#print(max(abs(rho_GR - arho_GR) / arho_GR ))

output(r, rho_GR, Filename)
