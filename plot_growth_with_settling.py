# Script to plot growth rate of instability

import numpy as np
import math
import matplotlib.pyplot as plt

# Set parameters for run
step_size=1
n_steps=400
near_0=0.00001
w=0.002
R=10
N=1

# Set up domain coordinates
xlist = np.linspace(-n_steps*step_size,n_steps*step_size,n_steps)
ylist = np.linspace(-n_steps*step_size,n_steps*step_size,n_steps)
X, Y = np.meshgrid(xlist, ylist)

# Initialise output vectors
roots_out_real=np.zeros((n_steps,n_steps))
asin_out=np.zeros((n_steps,n_steps))
roots_out_imag=np.zeros((n_steps,n_steps))

# Calculate growth at each grid cell
for k_i,k in enumerate(xlist):
	for l_i,l in enumerate(ylist):
		p=[1,-1j*w*k,l**2/(k**2+l**2)*(1-1/R),-l**2/(k**2+l**2)*w*1j*k] # Expression for polynomial
		roots_out_real[k_i,l_i]=max(np.real(np.roots(p)))
		roots_out_imag[k_i,l_i]=max(np.imag(np.roots(p)))
		asin_out[k_i,l_i]=l/math.sqrt(k**2+l**2)

# Plot outputs
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(top=0.85)
plt.subplot(2, 1, 1)
plt.imshow(roots_out_real,origin='lower',extent=[-n_steps*step_size,n_steps*step_size,-n_steps*step_size,n_steps*step_size])
plt.hot()
plt.colorbar()
plt.title('max(Real(lambda)) W='+str(w) +' R='+str(R))
plt.xlabel('l')
plt.ylabel('k')

plt.subplot(2, 1, 2)
plt.imshow(roots_out_imag,origin='lower',extent=[-n_steps*step_size,n_steps*step_size,-n_steps*step_size,n_steps*step_size])
plt.hot()
plt.colorbar()
plt.title('max(Imag(lambda))')
plt.xlabel('l')
plt.ylabel('k')

# Save output figure
plt.savefig('roots_out w='+str(w)+'.pdf',format='pdf')


# https://docs.scipy.org/doc/numpy/reference/generated/numpy.roots.html#numpy.roots