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

from gamst import snapshots
import numpy as np
import numba as nb
from numba import cuda

# read a snapshot file with certain format
class read:
	def __init__(self, filename):
		file_array=filename.split(".")
		if file_array[len(file_array)-1]=="mst":
			self.data = snapshots.read_mst.read_mst(filename)
		self.compute_properties = {'temperature':False, 'pressure':False, 'momentum':False, 'potential':False, 'stress_tensor':False}
		self.variant = {'position':True, 'velocity':True, 'type':False, 'mass':False, 'image':True, 'box':False, 
						'force':True, 'potential':True, 'virial':True, 'bond':False, 'angle':False, 'dihedral':False}		
		# system information
		self.npa = self.data.num_particles
		self.pitch = (self.npa + (16 - (self.npa & 15)))
		self.dt=0.001
		self.timestep = self.data.timestep
		self.dimension = self.data.dimension
		
		self.particle_set = []
		self.comp_info = []
		self.plist = []
		self.typemap=[]
		self.bond = None
		self.angle = None
		self.dihedral = None		
	
		# host arrays		
		self.pos = np.zeros([self.npa, 4], dtype=np.float32)
		self.vel = np.zeros([self.npa, 4], dtype=np.float32)
		self.image = np.zeros([self.npa, 3], dtype=np.int32)
		self.box = np.asarray(self.data.box, dtype=np.float32)
		self.tag = np.zeros(self.npa, dtype=np.int32)		
		self.rtag = np.zeros(self.npa, dtype=np.int32)		
		
		
		for i in range(0, self.npa):
			self.pos[i][0] = self.data.position[i][0]
			self.pos[i][1] = self.data.position[i][1]
			self.pos[i][2] = self.data.position[i][2]
			type_i = self.data.type[i]
			type_id = self.add_name_to_id(type_i)
			self.pos[i][3] = np.float32(type_id)
			self.tag[i] = i 
			self.rtag[i] = i

			
		for i in range(0, self.npa):
			if len(self.data.velocity)==self.npa:
				self.vel[i][0] = self.data.velocity[i][0]
				self.vel[i][1] = self.data.velocity[i][1]
				self.vel[i][2] = self.data.velocity[i][2]
			if len(self.data.mass)==self.npa:
				self.vel[i][3] = self.data.mass[i]
			else:
				self.vel[i][3] = 1.0		

		if len(self.data.image)==self.npa:
			for i in range(0, self.npa):
				self.image[i][0] = self.data.image[i][0]
				self.image[i][1] = self.data.image[i][1]
				self.image[i][2] = self.data.image[i][2]

		self.force = np.zeros([self.npa, 3], dtype=np.float32)
		self.virial_potential = np.zeros([self.npa, 2], dtype=np.float32)
		self.ntypes=len(self.typemap)
		
		if len(self.data.bond)>0:
			self.bond = snapshots.bonded_data.bond_data(self, self.data.bond)
			
		if len(self.data.angle)>0:
			self.angle = snapshots.bonded_data.angle_data(self, self.data.angle)

		if len(self.data.dihedral)>0:
			self.dihedral = snapshots.bonded_data.dihedral_data(self, self.data.dihedral)			

		# device arrays		
		self.d_pos = cuda.to_device(self.pos)
		self.d_box = cuda.to_device(self.box)
		self.d_image = cuda.to_device(self.image)
		self.d_vel = cuda.to_device(self.vel)
		self.d_force = cuda.to_device(self.force)
		self.d_virial_potential = cuda.to_device(self.virial_potential)
		self.d_tag = cuda.to_device(self.tag)
		self.d_rtag = cuda.to_device(self.rtag)
		self.d_velo = None

	# add type name and convert it to type index		
	def add_name_to_id(self, typename):
		for i in range(0, len(self.typemap)):
			if self.typemap[i] == typename:
				return i
		self.typemap.append(typename)
		return len(self.typemap)-1

	# convert type name to type index		
	def convert_name_to_id(self, typename):
		for i in range(0, len(self.typemap)):
			if self.typemap[i] == typename:
				return i
		raise RuntimeError('Error! type '+typename+' is not existed')
		

	def find_particle_set(self, group):
		for i in self.particle_set:
			if i.group==group:
				return i
			
	# find an existed comp_info with same group
	def find_comp_info(self, group):
		for i in self.comp_info:
			if i.group==group:
				return i
			
	def find_plist(self, rcut):
		for i in self.plist:
			if i.rcut>=rcut:
				return i			
			
