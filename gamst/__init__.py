'''
GALAMOST - GPU-Accelerated Large-Scale Molecular Simulation Toolkit
Version 5
COPYRIGHT
	GALAMOST Copyright (c) (2020) G.U.G.D.
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of GALAMOST do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with GALAMOST are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu et al.,
	J. Comput. Chem. 2013,34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu, 
	Email: ylzhu@galamost.com
'''

import sys
import ctypes
import os
import zipfile

curr_file = os.path.abspath(__file__)
poetry_path = curr_file.replace("gamst/__init__.py", "poetry")
import json
import random
if os.path.exists(poetry_path):
	tang = random.randint(0, 56)
	tang = tang * 1000
	filename = 'poet.tang.' + repr(tang) + '.json'
	f1 = zipfile.ZipFile(poetry_path+'/poetry.zip')
	f1.extract(filename)	
	with open(filename, 'r', encoding='utf8')as fp:
		json_data = json.load(fp)
	x = random.randint(0, 999)
	#json1 = json.dumps(json_data[x], indent=1, ensure_ascii=False)
	# print (type(json1))
	# print (json1)
	ser_dic = json_data[x]
	print("Read a Chinese Tang poem and wait a few seconds for the program initialization")
	print(ser_dic['title'] + '  作者:  ' + ser_dic['author'])
	for i in ser_dic['paragraphs']:
		print(i)
	os.remove(filename)

from gamst import application
from gamst import chare
from gamst import force
from gamst import integration
from gamst import plist
from gamst import snapshot
from gamst import tinker
from gamst import dump

from numba import cuda
from optparse import OptionParser



GALAMOST_VERSION="5"

print("GALAMOST ", GALAMOST_VERSION )
print("GALAMOST - GPU-Accelerated Large-Scale Molecular Simulation Toolkit")
print("COPYRIGHT" )
print("	GALAMOST Copyright (c) (2020) G.U.G.D." )
print("LICENSE" )
print("	This program is a free software: you can redistribute it and/or " )
print("	modify it under the terms of the GNU General Public License." ) 
print("	This program is distributed in the hope that it will be useful, " )
print("	but WITHOUT ANY WARRANTY; without even the implied warranty of" )
print("	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE." ) 
print("	See the General Public License v3 for more details." )
print("	You should have received a copy of the GNU General Public License" )
print("	along with this program. If not, see <http://www.gnu.org/licenses/>." )
print("DISCLAIMER" )
print("	The authors of GALAMOST do not guarantee that this program and its " )
print("	derivatives are free from error. In no event shall the copyright " )
print("	holder or contributors be liable for any indirect, incidental," ) 
print("	special, exemplary, or consequential loss or damage that results " )
print("	from its use. We also have no responsibility for providing the " )
print("	service of functional extension of this program to general users." )
print("USER OBLIGATION " )
print("	If any results obtained with GALAMOST are published in the scientific " )
print("	literature, the users have an obligation to distribute this program " )
print("	and acknowledge our efforts by citing the paper \"Y.-L. Zhu et al.," )
print("	J. Comput. Chem. 2013, 34, 2197-2211\" in their article." )
print("CORRESPONDENCE" )
print("	Dr. You-Liang Zhu," ) 
print("	Email: ylzhu@galamost.com" )

global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()


device_id = 0
if _options.gpu is not None:
	device_id = int(_options.gpu)

cuda.select_device(device_id)
	
d=cuda.cudadrv.driver.Device(devnum=device_id)
cc=str(d.compute_capability[0])+'.'+str(d.compute_capability[1])
id='['+str(d.id)+']'
print('GPU id', id, ' ', d.name, 'compute capability:', cc)	

if __name__ == '__main__':
    print ('main program run')
else:
    print ('gamst engine initialization')
	
	