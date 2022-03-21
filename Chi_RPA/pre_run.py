#!/usr/bin/env python3
import xmltodict as xml
import numpy as np
from sys import argv

# debug
# savedir ="/Users/Nevensky/Downloads"
# fname = "/".join([savedir,"MoS2.xml"])
# fname = "/".join(savedir,"data-file-schema.xml")
Ha2eV = 27.211386246

try:
	configFile = argv[1]
except IndexError:
	configFile = './config.in'

execdir ='/'.join(configFile.split('/')[:-1])
print('execdir: ',execdir)

def find_cell(cell):
	""" Returns lattice parameters for hexagonal 2D lattice."""
	a = np.linalg.norm(cell[0,:]) 
	c = np.linalg.norm(cell[2,:]) 
	Vcell = a**2 * c * np.sin(np.deg2rad(60)) 
	return a,c,Vcell

def find_NGd_and_Nocc(BS):
	NGd = BS[0]["npw"]
	Nocc = 0.0
	for idx,item in enumerate(BS):
		Npw_i = item["npw"]
		if NGd > Npw_i:
			NGd = Npw_i
		# Nocc_i = np.asarray(item["occupations"]["#text"],dtype=np.float64,delim="\n")
		Nocc_i = np.fromstring(item["occupations"]["#text"], dtype=np.float64, sep='\n' ).sum()
		if Nocc_i > Nocc:
			Nocc = Nocc_i
	return np.int32(NGd), np.int32(Nocc)	

def initalChecks(dat):
	calc_type = dat["qes:espresso"]["input"]["control_variables"]["calculation"]
	symm_noinv = dat["qes:espresso"]["input"]["symmetry_flags"]["noinv"]
	symm_symorph = dat["qes:espresso"]["input"]["symmetry_flags"]["force_symmorphic"]
	tmpdir = dat["qes:espresso"]["input"]["control_variables"]["outdir"]

	if calc_type != "nscf":
		print("calc_type: ",calc_type)
		raise Exception("Invalid calculation type (should be nscf!).")

	if symm_symorph=="false":
		print("WARNING: force_symporphic=.true. not specified in PW input! Symmetry errors might occur.")

	if symm_noinv=="false":
		print("WARNING: noinv=.true. not specified in PW input! Symmetry errors might occur.")

	if tmpdir!=savedir:
		print("WARNING: outdir specified in data-file-schema.xml and config.in do not match")
		print("tmpdir (data-file-schema.xml): ",tmpdir)
		print("tmpdir (config.in)           : ",savedir)	

	return None


def parseParameters(dat):
	global EFermi,Nbands,NelQE,NkI,Nsymm,NG,NGd,Nocc,a,c,Vcell
	Nk_x = dat["qes:espresso"]["input"]["k_points_IBZ"]["monkhorst_pack"]["@nk1"]
	Nk_y = dat["qes:espresso"]["input"]["k_points_IBZ"]["monkhorst_pack"]["@nk2"]
	Nk_z = dat["qes:espresso"]["input"]["k_points_IBZ"]["monkhorst_pack"]["@nk3"]
	EFermi = Ha2eV*np.float64(dat["qes:espresso"]["output"]["band_structure"]["fermi_energy"]) # converted from Ha to eV
	Nbands = np.int32(dat["qes:espresso"]["output"]["band_structure"]["nbnd"])
	NelQE = np.int32(np.float64(dat["qes:espresso"]["output"]["band_structure"]["nelec"]))
	NkI = np.int32(dat["qes:espresso"]["output"]["band_structure"]["nks"])
	Nsymm = np.int32(dat["qes:espresso"]["output"]["symmetries"]["nsym"])
	NG = np.int32(dat["qes:espresso"]["output"]["basis_set"]["ngm"])
	BS = dat["qes:espresso"]["output"]["band_structure"]["ks_energies"]
	cell_info = dat["qes:espresso"]["output"]["atomic_structure"]

	alat = np.float64(cell_info["@alat"])
	a1 = np.fromstring(cell_info["cell"]["a1"],dtype=np.float64,sep=" ")
	a2 = np.fromstring(cell_info["cell"]["a2"],dtype=np.float64,sep=" ")
	a3 = np.fromstring(cell_info["cell"]["a3"],dtype=np.float64,sep=" ")
	cell = np.array([a1,a2,a3]) # in units of [alat]

	a,c,Vcell = find_cell(cell)
	NGd, Nocc = find_NGd_and_Nocc(BS)
	
	print("kmesh: ",Nk_x," x ",Nk_y," x ",Nk_z)
	print("alat:",alat)
	print("a: ",a," [bohr] ","c: ",c," [bohr] ","\nVolume:",Vcell,"[bohr^3](hexagonal 2D)")
	print("NG:",NG,"\nNGd: ",NGd)
	print("Nsymm: ",Nsymm)
	print("Nocc: ",Nocc)
	print("NelQE: ",NelQE,"\nNbands: ",Nbands,"\nNkI: ",NkI)
	print("Fermi Energy: ",EFermi," [eV]")
	
	return None

default_config ="""&directories
 rundir    = ''
 savedir   = ''
 scf_file  = ''
 band_file = ''
/
&system
 lf        ='z'          ! Crystal local field effect included in z for lf=1 or in x,y,z direction lf=3
 jump      = 1           ! jump = 1 skips reloading wfns in IBZ for all bands m (occ.) and n (virt.)
 omin      = 1.0d-5      ! [Hartree] frequency range, lower bound
 omax      = 2.0D0       ! [Hartree] frequency range, upper bound
 qmin      = 2		     ! minimum transfer wave-vector
 qmax      = 2           ! maximum transfer wave-vector 
/
&config
 NG       = 0           ! Total number of G vectors
 NGd      = 0           ! Number of coefficients CG; should be less than minimum number of coefficients over all evc.n
 NkI      = 0           ! Number of wave vectors in IBZ
 Nband    = 0           ! Number of bands (unit cell)
 Nocc     = 0           ! Number of occupied bands (unit cell)
 NelQE    = 0           ! Number of electrons(unit cell)
 No       = 2001        ! Number of frequencies
 Nlfd     = 100         ! Number of loacal field vectors (eventually should be fully dynamically allocated)
/
&parameters
 Efermi   = 0.0000       ! [eV]      Fermi energy 
 eps      = 1.0d-4       ! 1.0D-4    threshold
 T        = 0.01         ! [eV]      Temperature 
 eta      = 0.05         ! [eV]      Damping i\eta
 Ecut     = 0.00         ! [Hartree] Cutoff energy for crystal local field calculations , for Ecut=0 S matrix is a scalar ?
 a0       = 0.00000      ! [a.u.]    unit cell parameter in parallel direction 
 c0       = 0.00000      ! [a.u.]    unit cell parameter in perependicular direction 
 Vcell    = 0.00000      ! [a.u.^3]  Unit-cell volume 
/
&parallel
 Nthreads = 1            ! Number of OpenMP threads
/
"""


if __name__ == '__main__':
	NG_comment    = '      ! Total number of G vectors'
	NGd_comment   = '       ! Number of coefficients CG; should be less than minimum number of coefficients over all evc.n'
	NelQE_comment = '         ! Number of electrons(unit cell)'
	Nband_comment = '         ! Number of bands (unit cell)'
	NkI_comment   = '           ! Number of wave vectors in IBZ'
	Nocc_comment  = '          ! Number of occupied bands (unit cell)'
	Vcell_comment = '    ! [a.u.^3]  Unit-cell volume '
	EFermi_comment= '      ! [eV]      Fermi energy '
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
					fname = "/".join([savedir,"data-file-schema.xml"])
					if savedir == '':
						print('WARNING: QE savedir not defined.')
						exit()
					lns_new.append(ln)
					try:
						with open(fname,"r") as xmlFile:
							dat = xml.parse(xmlFile.read())
						initalChecks(dat)
						parseParameters(dat)
					except FileNotFoundError:
						print("WARNING: data-file-schema.xml not found.")
						exit()
				elif 'rundir' in ln:
					rundir = ln.split("\'")[1]
					if rundir == '':
						print('WARNING: QE rundir not defined.')
						exit()
					lns_new.append(ln)
				elif 'scf_file' in ln:
					scf_file = ln.split("\'")[1]
					if scf_file=='':
						print('WARNING: scf_file not defined.')
						exit()
					lns_new.append(ln)
				elif 'nsc_file' in ln:
					nsc_file = ln.split("\'")[1]
					if nsc_file=='':
						print('WARNING: nsc_file not defined.')
						exit()
					lns_new.append(ln)
				elif 'band_file' in ln:
					band_file = ln.split("\'")[1]
					if band_file =='':
						print('WARNING: band_file not defined.')
						exit()
					lns_new.append(ln)
				elif 'NelQE' in ln:
					ln2 = ' NelQE    = {} {}\n'.format(NelQE,NelQE_comment)
					lns_new.append(ln2)
				elif 'Nband' in ln:
					ln2 = ' Nband    = {} {}\n'.format(Nbands,Nband_comment)
					lns_new.append(ln2)
				elif 'Nocc' in ln:
					ln2 = ' Nocc     = {} {}\n'.format(Nocc,Nocc_comment)
					lns_new.append(ln2)
				elif 'NkI' in ln:
					ln2 = ' NkI      = {} {}\n'.format(NkI,NkI_comment)
					lns_new.append(ln2)
				elif 'NG ' in ln:
					ln2 = ' NG       = {} {}\n'.format(NG,NG_comment)
					lns_new.append(ln2)
				elif 'NGd ' in ln:
					ln2 = ' NGd      = {} {}\n'.format(NGd,NGd_comment)
					lns_new.append(ln2)
				elif 'Efermi ' in ln:
					ln2 = ' Efermi   = {:.4f} {}\n'.format(EFermi,EFermi_comment)
					lns_new.append(ln2)
				elif 'Vcell' in ln:
					ln2 = ' Vcell    = {:.4f} {}\n'.format(Vcell,Vcell_comment)
					lns_new.append(ln2)
				elif 'a0' in ln:
					ln2 = ' a0       = {:.4f} {}\n'.format(a,a0_comment)
					lns_new.append(ln2)
				elif 'c0' in ln:
					ln2 = ' c0       = {:.4f} {}\n'.format(c,c0_comment)
					lns_new.append(ln2)
				else:
					lns_new.append(ln)
			print('\n'+"".join(lns_new))

		with open(configFile,'w') as fl:
			fl.writelines(lns_new)

	except FileNotFoundError:
		with open(configFile,'w') as fl2:
			fl2.writelines(default_config)
		print('config.in not found, generated default config.\n Edit the default values and rerun this script.')
	

