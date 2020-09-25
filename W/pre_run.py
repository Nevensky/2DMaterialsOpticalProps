#!/usr/bin/env python3
import os
import subprocess
from sys import argv

# config_file = sys.argv[1]
configFile = '/Users/nevensky/Repositories/2d-quasiparticle-optical-properties/W/config.in'
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

# confs= """
# &CONFIGURATION
#  NGd={NGd}
#  /
#  """.format(NGd=NGd)
# print(confs)

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
			ln2 = 'NGd = {}\n'.format(NGd)
			lns_new.append(ln2)
		else:
			lns_new.append(ln)
print('\n'+"".join(lns_new))

if savedir != '':
	os.system('iotk convert {} {}'.format(savedir+'/gvectors.dat',rundir+'/gvectors.xml'))

with open(configFile,'w') as fl:
	fl.writelines(lns_new)



