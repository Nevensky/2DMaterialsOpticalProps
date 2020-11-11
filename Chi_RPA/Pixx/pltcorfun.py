#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt4Agg')
plt.style.use('ggplot')

plt.figure(num=None, figsize=(10, 6), dpi=100)

fl= './Corrfun_xx'

lines = []
freqs = []
imvals,revals = [],[]
with open(fl,'r') as f:
	lines = f.readlines()
	for i,ln in enumerate(lines):
		if i%2==0:
			freqs.append(float(ln.split()[1]))
		else:
			revals.append(float(ln.split()[0]))
			imvals.append(float(ln.split()[1]))

# data = np.array([freqs,revals])

plt.xlabel('frequency [Ha]')
plt.ylabel(r'$\frac{-\Im\left[ \chi (\omega,\vec{G_1},\vec{G_2}) \right]}{\pi}$')
plt.plot(freqs,revals,color='goldenrod',label='real part')
plt.legend(loc='upper left')
plt.twinx()
plt.tick_params(axis='y', labelcolor='darkseagreen')
plt.plot(freqs,imvals,color='darkseagreen',label='imaginary part')
plt.legend(loc='upper right')
plt.show()
plt.close()

fl2= './Pi_RPA_xx'
freqsEV,valsRe,valsIm = [],[],[]
with open(fl2,'r') as f:
	lines = f.readlines()
	for i,ln in enumerate(lines):
		freq = float(ln.split()[0])
		ln2 = ln.split()[1].lstrip('(').rstrip(')').split(',')
		re,im = float(ln2[0]),float(ln2[1])
		freqsEV.append(freq)
		valsRe.append(re)
		valsIm.append(im)
plt.figure(num=None, figsize=(10, 6), dpi=100)
plt.xlabel(r'frequency [eV]')
plt.ylabel(r'$\Re \left(\Pi_{\mu\nu}\right)$')
# plt.xlim([0,5])
plt.plot(freqsEV,valsRe,color='palevioletred')
plt.show()

plt.figure(num=None, figsize=(10, 6), dpi=100)
plt.xlabel(r'frequency [eV]')
plt.ylabel(r'$\Im \left(\Pi_{\mu\nu}\right)$')
# plt.xlim([0,5])
plt.plot(freqsEV,valsIm,color='cornflowerblue')
plt.show()