#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt4Agg')
plt.style.use('ggplot')

plt.figure(num=None, figsize=(10, 6), dpi=100)

# fl= './Corrfun_xx'
fl= '/Users/Nevensky/Downloads/Pixx/Corrfun_xx'


Nlf,Nlf_counted = 0,0
lines = []
freqs = []
imvals,revals = [],[]
with open(fl,'r') as f:
	lines = f.readlines()
	for i,ln in enumerate(lines):
		# if i%2==0:
		if "omega" in ln:
			freqs.append(float(ln.split()[1]))
			Nlf_counted += 1
			if Nlf_counted==2:
				Nlf = i-1
				print('Nlf: ',int(np.sqrt(Nlf))) 
		else:
			revals.append(float(ln.split()[0]))
			imvals.append(float(ln.split()[1]))
freqs = np.array(freqs)
revals = np.array(revals)
imvals = np.array(imvals)

print('revals shape',revals[::].shape)
print('revals assumed shape',(freqs.shape[0])*Nlf)
# data = np.array([freqs,revals])
# plt.figure(num=None, figsize=(10, 6), dpi=100)
plt.xlabel('frequency [Ha]')
plt.ylabel(r'$\frac{-\Im\left[ \chi (\omega,\vec{G_1},\vec{G_2}) \right]}{\pi}$')
plt.plot(freqs,revals[1:revals[::].shape[0]:Nlf],color='goldenrod',label='real part')
plt.legend(loc='upper left')
plt.twinx()
plt.tick_params(axis='y', labelcolor='darkseagreen')
plt.ylim([min(revals[1:revals[::].shape[0]:Nlf]),max(revals[1:revals[::].shape[0]:Nlf])])
plt.plot(freqs,imvals[1:imvals[::].shape[0]:Nlf],color='darkseagreen',label='imaginary part')
plt.legend(loc='upper right')
plt.show()
plt.close()

print("Plotting Corrfun for first 12 local field vectors.")
for j in range(1,12+1):
	plt.subplot(4,3,j)
	plt.plot(freqs,revals[j:revals.shape[0]:Nlf],color='goldenrod',label='real part')
	plt.twinx()
	plt.ylim([min(revals[j:revals[::].shape[0]:Nlf]),max(revals[j:revals[::].shape[0]:Nlf])])
	plt.tick_params(axis='y', labelcolor='darkseagreen')
	plt.plot(freqs,imvals[j:imvals.shape[0]:Nlf],color='darkseagreen',label='imaginary part')
			
plt.show()
plt.close()

# for j in range(1,Nlf):
# 	if j<=15:
# 		plt.subplot(3,5,j)
# 		if j==1:
# 			plt.xlabel('frequency [Ha]')
# 			plt.ylabel(r'$\frac{-\Im\left[ \chi (\omega,\vec{G_1},\vec{G_2}) \right]}{\pi}$')
# 			plt.legend(loc='upper left')
# 		plt.plot(freqs,revals[j:revals[::].shape[0]:Nlf],color='goldenrod',label='real part')
		
# 		plt.twinx()
# 		plt.tick_params(axis='y', labelcolor='darkseagreen')
# 		plt.plot(freqs,imvals[j:revals[::].shape[0]:Nlf],color='darkseagreen',label='imaginary part')
# 		if j==1:
# 			plt.legend(loc='upper right')
plt.show()
plt.close()


# fl2= '/Users/Nevensky/Downloads/Pixx/Pi_RPA_xx'
fl2= '/Users/Nevensky/Downloads/Pixx/fort.99'
# fl2= '/Users/Nevensky/Downloads/fort.99'
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
plt.plot(freqsEV,valsRe,color='palevioletred')
plt.show()

# zoom in 5eV
plt.figure(num=None, figsize=(10, 6), dpi=100)
plt.xlabel(r'frequency [eV]')
plt.ylabel(r'$\Re \left(\Pi_{\mu\nu}\right)$')
plt.plot(freqsEV,valsRe,color='palevioletred')
plt.xlim([0.,5.])
plt.show()

plt.figure(num=None, figsize=(10, 6), dpi=100)
plt.xlabel(r'frequency [eV]')
plt.ylabel(r'$\Im \left(\Pi_{\mu\nu}\right)$')
plt.plot(freqsEV,valsIm,color='cornflowerblue')
plt.show()

# zoom in 5eV
plt.figure(num=None, figsize=(10, 6), dpi=100)
plt.xlim([0.,5.])
plt.xlabel(r'frequency [eV]')
plt.ylabel(r'$\Im \left(\Pi_{\mu\nu}\right)$')
plt.xlim([0,5])
plt.plot(freqsEV,valsIm,color='cornflowerblue')
plt.show()