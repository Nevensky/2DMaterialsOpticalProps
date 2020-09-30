#!/usr/bin/env python3
import os
import subprocess
from sys import argv

try:
	configFile = argv[1]
except IndexError:
	configFile = './config.in'

rundir ='/'.join(configFile.split('/')[:-1])
print('rundir: ',rundir)

def findMin_NGd(savedir):
	paths = []
	igwxs = []
	for root, dirs, files in os.walk(savedir):
		for f in files:
			if (os.path.splitext(f)[1] == '.dat' and 'evc' in f):
				fullpath = os.path.join(root, f)
				paths.append(fullpath)
				igwx = subprocess.check_output("grep --line-buffered --text \"igwx\" {}".format(fullpath), shell=True)
				igwx = int(igwx.decode('ascii').split()[2].split("\"")[1])
				igwxs.append(igwx)

	NGd = min(igwxs)

	print('paths: ',paths)
	print('igwxs: ',igwxs)
	print("Minimal igwx --> NGd =",NGd)

	return NGd

default_config = """&directories
rundir    = ''
savedir   = ''
scf_file  = ''
band_file = ''
/
&system
lf        = 1           ! Crystal local field effect included in z for lf=1 or in x,y,z direction lf=3
loss      = 1           ! ??
jump      = 1           ! za 1 preskace trazenje wfn. u IBZ za sve bands m i n
omin      = 1.0D-5      ! [Hartree] frequency range, lower bound
omax      = 2.0D0       ! [Hartree] frequency range, upper bound
/
&config
 NGd      = 0           ! number of wave vectors in IBZ
 NkI      = 0           ! number of bands
 Nband    = 0           ! Number of electrons(unit cell)
 NelQE    = 0           ! number of coefficients CG shulod be less than minimum number of coefficients all over all evc.n 
 NG       = 0           ! total number of G vectors  
 no       = 2001        ! broj frekvencija
 nq       = 1           ! broj valnih vektora tu je 2 jer je rucno paralelizirano!
 Nlfd     = 0           ! dimenzija polja za local field zasto prozivoljno 50, ne moze se znati unaprijed
/
&parameters
Efermi   = 0.0          ! [eV] Fermi en. 
a0       = 0.0          ! [a.u.]  unit cell parameter in parallel direction 
c0       = 0.0          ! [a.u.]  unit cell parameter in perependicular direction 
eps      = 1.0D-4       ! 1.0D-4 threshold
T        = 0.01         ! [eV] temperature 
eta      = 0.05         ! damping i\\eta
Ecut     = 0.0          ! [Hartree] cutoff energy for crystal local field calculations , for Ecut=0 S matrix is a scalar ?
Vcell    = 0.0          ! [a.u.^3] unit-cell volume 
aBohr    = 0.0          ! [a.u.] unit cell parameter in perpendicular direction (z-separation between supercells)   
/
&parallel
Nthreads = 1            ! number of OpenMP threads
/
"""


NGd_comment ='       ! number of wave vectors in IBZ'

try:
	lns_new = []
	savedir = ''
	with open(configFile,'r') as fl:
		lns = fl.readlines()
		for ln in lns:
			if 'savedir' in ln:
				savedir = ln.split("\'")[1]
				# print(savedir)
				lns_new.append(ln)
			elif 'NGd' in ln:
				NGd = 0
				NGd = findMin_NGd(savedir)
				ln2 = ' NGd      = {} {}\n'.format(NGd,NGd_comment)
				lns_new.append(ln2)
			else:
				lns_new.append(ln)
	print('\n'+"".join(lns_new))
	
	if savedir != '':
		os.system('iotk convert {} {}'.format(savedir+'/gvectors.dat',rundir+'/gvectors.xml'))
	
	with open(configFile,'w') as fl:
		fl.writelines(lns_new)

except FileNotFoundError:
	with open(configFile,'w') as fl2:
		fl2.writelines(default_config)
	print('config.ini not found, generated default config.\n Edit the default values and rerun this script.')

