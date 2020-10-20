#!/usr/bin/env python3
import os
import subprocess
from sys import argv

try:
	configFile = argv[1]
except IndexError:
	configFile = './config.in'

execdir ='/'.join(configFile.split('/')[:-1])
print('execdir: ',execdir)

def findMin_NGd(savedir):
	""" Finds the smallest number of G vectors for each k point
	in QE nsc/bands calculation. """
	paths = []
	igwxs = []
	for root, dirs, files in os.walk(savedir):
		for f in files:
			if (os.path.splitext(f)[1] == '.dat' and 'evc' in f):
				fullpath = os.path.join(root, f)
				paths.append(fullpath)
				print(fullpath)
				igwx = subprocess.check_output("grep --line-buffered --text \"igwx\" {}".format(fullpath), shell=True)
				igwx = igwx.split('/>'.encode('ascii'))[0]#print(igwx)
				igwx = int(igwx.decode('ascii').split()[2].split("\"")[1])
				igwxs.append(igwx)
				print(igwx)

	NGd = min(igwxs)

	print('paths: ',paths)
	print('igwxs: ',igwxs)
	print("Minimal igwx --> NGd =",NGd)

	return NGd


def findNelQE(rundir,scf_file):
	""" Finds number of occupied bands in QE scf output. """
	NelQE = 0
	path = rundir+scf_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'number of electrons' in ln:
				NelQE = round(float(ln.split()[4]))
				break
	return NelQE

# def findNkI(rundir,nsc_file):
# 	""" Finds number of k-points. """
# 	NkI = 0
# 	path = rundir+nsc_file
# 	with open(path,'r') as f:
# 		lns = f.readlines()
# 		for idx,ln in enumerate(lns):
# 			if 'number of k points=' in ln:
# 				NkI = int(ln.split()[4])
# 				break
# 	return NkI

def findNband(rundir,band_file):
	""" Finds number of occupied bands in QE scf output. """
	Nband = 0
	path = rundir+band_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'nbnd=' in ln:
				Nband = round(float(ln.split()[2].rstrip(',')))
				break
	return Nband

def findNkI(rundir,band_file):
	""" Finds number of k-points.  """
	NkI = 0
	path = rundir+band_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'nks=' in ln:
				print(ln)
				NkI = int(ln.split()[4].rstrip(','))
				break
	return NkI

def findNG(rundir,scf_file):
	""" Finds total number of G-vectors in QE scf. """
	NG = 0
	path = rundir+scf_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'G-vecs:' in ln:
				NG = round(float(lns[idx+1].split()[4]))
				break
	return NG
	

def findVcell(rundir,scf_file):
	""" Finds number of occupied bands in QE scf output. """
	Vcell = 0.0
	path = rundir+scf_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'unit-cell volume' in ln:
				Vcell = round(float(ln.split()[3]),4)
				break
	return Vcell

def findaBohr(rundir,scf_file):
	""" Finds number of occupied bands in QE scf output. """
	alat = 0.0
	path = rundir+scf_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'lattice parameter (alat)' in ln:
				alat = round(float(ln.split()[4]),8)
				break
	return alat

def findCelldm(rundir,scf_file):
	""" Finds number of occupied bands in QE scf output. """
	celldm = 6*[0.0]
	path = rundir+scf_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'celldm(1)' in ln:
				celldm[0] = round(float(ln.split()[1]),6)
				celldm[1] = round(float(ln.split()[3]),6)
				celldm[2] = round(float(ln.split()[5]),6)
			elif 'celldm(4)' in ln:
				celldm[3] = round(float(ln.split()[1]),6)
				celldm[4] = round(float(ln.split()[3]),6)
				celldm[5] = round(float(ln.split()[5]),6)
				break
	print(celldm)
	return celldm

def findNocc(rundir,scf_file):
	""" Finds number of occupied bands in QE scf output. """
	Nocc = 0
	path = rundir+scf_file
	with open(path,'r') as f:
		lns = f.readlines()
		for idx,ln in enumerate(lns):
			if 'occupation numbers' in ln:
				j = idx + 1
				while (True):
					if lns[j] != '\n':
						ln2 = lns[j]
						Nocc += sum([round(float(i)) for i in ln2.split()])
						j += 1
					else:
						break
				break
	return Nocc

default_config = """&directories
 rundir    = ''
 savedir   = ''
 scf_file  = ''
 band_file = ''
/
&system
 lf        = 'z'         ! Crystal local field effect dirrection, available: 'z' or 'xyz'
 jump      = 1           ! jump = 1 skips reloading wfns in IBZ for all bands m (occ.) and n (virt.)
 omin      = 1.0D-5      ! [Hartree] frequency range, lower bound
 omax      = 2.0D0       ! [Hartree] frequency range, upper bound
/
&config
 NG       = 0           ! Total number of G vectors  
 NGd      = 0           ! Number of coefficients CG shulod be less than minimum number of coefficients all over all evc.n 
 NkI      = 0           ! Number of wave vectors in IBZ
 Nband    = 0           ! Number of bands (unit cell)
 Nocc     = 0           ! Number of occupied bands (unit cell)
 NelQE    = 0           ! Number of electrons(unit cell)
 no       = 2001        ! Number of frequencies
 nq       = 2           ! broj valnih vektora tu je 2 jer je rucno paralelizirano!
 Nlfd     = 0           ! dimenzija polja za local field zasto prozivoljno 50, ne moze se znati unaprijed
/
&parameters
 Efermi   = 0.0          ! [eV]      Fermi energy 
 a0       = 0.0          ! [a.u.]    unit cell parameter in parallel direction 
 c0       = 0.0          ! [a.u.]    unit cell parameter in perependicular direction 
 eps      = 1.0D-4       ! 1.0D-4    threshold
 T        = 0.01         ! [eV]      Temperature 
 eta      = 0.05         ! [eV]      Damping i\\eta
 Ecut     = 0.0          ! [Hartree] Cutoff energy for crystal local field calculations , for Ecut=0 S matrix is a scalar ?
 Vcell    = 0.0000       ! [a.u.^3]  Unit-cell volume 
/
&parallel
 Nthreads = 1            ! Number of OpenMP threads
/
"""

NG_comment    = '      ! Total number of G vectors'
NGd_comment   = '       ! Number of coefficients CG; should be less than minimum number of coefficients over all evc.n'
NelQE_comment = '         ! Number of electrons(unit cell)'
Nband_comment = '         ! Number of bands (unit cell)'
NkI_comment   = '           ! Number of wave vectors in IBZ'
Nocc_comment  = '          ! Number of occupied bands (unit cell)'
Vcell_comment = '    ! [a.u.^3]  Unit-cell volume '
Efermi_comment= '       ! [eV]      Fermi energy '
a0_comment    = '      ! [a.u.]    unit cell parameter in parallel direction '
c0_comment    = '     ! [a.u.]    unit cell parameter in perependicular direction '

try:
	lns_new = []
	savedir = ''
	with open(configFile,'r') as fl:
		lns = fl.readlines()
		for ln in lns:
			if 'savedir' in ln:
				savedir = ln.split("\'")[1]
				if savedir == '':
					print('WARNING QE savedir not defined.')
					exit()
				lns_new.append(ln)
			elif 'rundir' in ln:
				rundir = ln.split("\'")[1]
				if rundir == '':
					print('WARNING QE rundir not defined.')
					exit()
				lns_new.append(ln)
			elif 'scf_file' in ln:
				scf_file = ln.split("\'")[1]
				if scf_file=='':
					print('WARNING scf_file not defined.')
					exit()
				lns_new.append(ln)
			elif 'nsc_file' in ln:
				nsc_file = ln.split("\'")[1]
				if nsc_file=='':
					print('WARNING nsc_file not defined.')
					exit()
				lns_new.append(ln)
			elif 'band_file' in ln:
				band_file = ln.split("\'")[1]
				if band_file =='':
					print('WARNING band_file not defined.')
					exit()
				lns_new.append(ln)
			elif 'NelQE' in ln:
				NelQE = 0
				NelQE = findNelQE(rundir,scf_file)
				ln2 = ' NelQE    = {} {}\n'.format(NelQE,NelQE_comment)
				lns_new.append(ln2)
			elif 'Nband' in ln:
				Nband = 0
				Nband = findNband(rundir,band_file)
				ln2 = ' Nband    = {} {}\n'.format(Nband,Nband_comment)
				lns_new.append(ln2)
			elif 'Nocc' in ln:
				Nocc = 0
				Nocc = findNocc(rundir,scf_file)
				ln2 = ' Nocc     = {} {}\n'.format(Nocc,Nocc_comment)
				lns_new.append(ln2)
			elif 'NkI' in ln:
				NkI = 0
				NkI = findNkI(rundir,band_file)
				ln2 = ' NkI      = {} {}\n'.format(NkI,NkI_comment)
				lns_new.append(ln2)
			elif 'NG ' in ln:
				NG = 0
				NG = findNG(rundir,scf_file)
				ln2 = ' NG       = {} {}\n'.format(NG,NG_comment)
				lns_new.append(ln2)
			elif 'NGd ' in ln:
				NGd = 0
				NGd = findMin_NGd(savedir)
				ln2 = ' NGd      = {} {}\n'.format(NGd,NGd_comment)
				lns_new.append(ln2)
			elif 'Vcell' in ln:
				Vcell = 0.0
				Vcell = findVcell(rundir,scf_file)
				ln2 = ' Vcell    = {:.4f} {}\n'.format(Vcell,Vcell_comment)
				lns_new.append(ln2)
			elif 'a0' in ln:
				a0 = 0.0
				a0 = findCelldm(rundir,scf_file)[0]
				ln2 = ' a0       = {:.4f} {}\n'.format(a0,a0_comment)
				lns_new.append(ln2)
			elif 'c0' in ln:
				c0 = 0.0
				c0 = findCelldm(rundir,scf_file)[2]*a0
				ln2 = ' c0       = {:.4f} {}\n'.format(c0,c0_comment)
				lns_new.append(ln2)
			else:
				lns_new.append(ln)
	print('\n'+"".join(lns_new))
	
	if savedir != '':
		os.system('iotk convert {} {}'.format(savedir+'/gvectors.dat',execdir+'/gvectors.xml'))
	
	with open(configFile,'w') as fl:
		fl.writelines(lns_new)

except FileNotFoundError:
	with open(configFile,'w') as fl2:
		fl2.writelines(default_config)
	print('config.in not found, generated default config.\n Edit the default values and rerun this script.')

