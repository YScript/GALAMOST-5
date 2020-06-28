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

# read a snapshot file with mst format
class read_mst:
	# define init function
	def __init__(self, filename):
		self.num_particles=0
		self.timestep=0
		self.dimension=3
		self.mst_read =False
		self.invariant_data = False
		self.variant_data = False
		self.read_indicator = {'box':False, 'position':False, 'type':False, 'image':False, 'mass':False, 
								'velocity':False, 'force':False, 'virial':False, 'bond':False, 'angle':False, 'dihedral':False}		
		self.init_data()	
		self.read_file(filename)
	
	def read_file(self, filename):
		with open(filename) as file_object:
			self.reset_params()
			for line in file_object:
				lin=line.strip('\n')
				line_array = lin.split()
#				print(line_array)
				if len(line_array)==0:
					continue
					
				if line_array[0]=="mst_end":
					break

				if line_array[0]=="mst_version" and float(line_array[1])== 1.0:
					print("info : read mst file with version 1.0")
					self.mst_read = True
					continue
						
				if self.mst_read:
					if line_array[0]=="invariant_data":
						self.invariant_data = True
						self.variant_data = False						
						continue
						
					if line_array[0]=="variant_data":
						self.invariant_data = False
						self.variant_data = True	
						continue
						
					if line_array[0]=="frame":
						if self.variant_data:
							self.init_data()
						else:
							raise RuntimeError('Error! mst files with multiple frames without the label of "variant_data"')
						continue

					if line_array[0]=="num_particles":
						self.reset_params()
						self.num_particles_read=True
						continue
						
					if line_array[0]=="timestep":
						self.reset_params()
						self.timestep_read=True
						continue
						
					if line_array[0]=="dimension":
						self.reset_params()
						self.dimension_read=True
						continue						

					if line_array[0]=="position":
						self.reset_params()				
						self.position_read=True
						if self.invariant_data:
							self.read_indicator['position'] = True						
						continue
					
					if line_array[0]=="type":
						self.reset_params()				
						self.type_read=True
						if self.invariant_data:
							self.read_indicator['type'] = True						
						continue
						
					if line_array[0]=="image":
						self.reset_params()
						self.image_read=True
						if self.invariant_data:
							self.read_indicator['image'] = True
						continue
						
					if line_array[0]=="mass":
						self.reset_params()
						self.mass_read=True
						if self.invariant_data:
							self.read_indicator['mass'] = True						
						continue
						
					if line_array[0]=="velocity":
						self.reset_params()
						self.mass_read=True
						if self.invariant_data:
							self.read_indicator['velocity'] = True
						continue
						
					if line_array[0]=="bond":
						self.reset_params()
						self.bond_read=True	
						if self.invariant_data:
							self.read_indicator['bond'] = True
						continue
						
					if line_array[0]=="angle":
						self.reset_params()
						self.angle_read=True
						if self.invariant_data:
							self.read_indicator['angle'] = True
						continue
						
					if line_array[0]=="dihedral":
						self.reset_params()
						self.dihedral_read=True
						if self.invariant_data:
							self.read_indicator['dihedral'] = True
						continue						
				
					if line_array[0]=="box":
						self.reset_params()
						self.box_read=True
						if self.invariant_data:
							self.read_indicator['box'] = True
						continue						
						
					# read data
					if self.num_particles_read and len(line_array) >= 1:
						self.num_particles=int(line_array[0])
						print("info : number of particles", self.num_particles)
					
					if self.timestep_read and len(line_array) >= 1:
						self.timestep=int(line_array[0])
						print("info : timestep", self.timestep)
						
					if self.dimension_read and len(line_array) >= 1:
						self.dimension=int(line_array[0])		
						print("info : dimension", self.dimension)
						
					if self.box_read and len(line_array) >= 3:
						self.box.append(float(line_array[0]))
						self.box.append(float(line_array[1]))
						self.box.append(float(line_array[2]))
						print("info : box size", line_array[0], line_array[1], line_array[2])
						
					if self.position_read and len(line_array) >= 3:
						self.position.append([float(line_array[0]), float(line_array[1]), float(line_array[2])])
						
					if self.type_read and len(line_array) >= 1:
						self.type.append(line_array[0])
						
					if self.image_read and len(line_array) >= 3:
						self.image.append([int(line_array[0]), int(line_array[1]), int(line_array[2])])
						
					if self.mass_read and len(line_array) >= 1:
						self.mass.append(float(line_array[0]))
						
					if self.velocity_read and len(line_array) >= 3:
						self.velocity.append([float(line_array[0]), float(line_array[1]), float(line_array[2])])
						
					if self.bond_read and len(line_array) >= 3:
						self.bond.append([line_array[0], int(line_array[1]), int(line_array[2])])

					if self.angle_read and len(line_array) >= 4:
						self.angle.append([line_array[0], int(line_array[1]), int(line_array[2]), int(line_array[3])])

					if self.dihedral_read and len(line_array) >= 5:
						self.dihedral.append([line_array[0], int(line_array[1]), int(line_array[2]), int(line_array[3]), int(line_array[4])])	
						
		# check data
		if not self.mst_read:
			raise RuntimeError('Error! the input file is not a mst file with version 1.0!')
		
		if len(self.position) != self.num_particles:
			raise RuntimeError('Error! number of position ', len(self.position), ' is not equal to the number of particles ', self.num_particles)	
		print("info :", len(self.position), "positions")

		if len(self.type) != self.num_particles:
			raise RuntimeError('Error! number of type ', len(self.type), ' is not equal to the number of particles ', self.num_particles)
		print("info :", len(self.type), "types")
		
		if len(self.image) > 0:
			if len(self.image) != self.num_particles:
				raise RuntimeError('Error! number of image ', len(self.image), ' is not equal to the number of particles ', self.num_particles)
			print("info :", len(self.image), "images")
					
		if len(self.mass) > 0:
			if len(self.mass) != self.num_particles:
				raise RuntimeError('Error! number of mass ', len(self.mass), ' is not equal to the number of particles ', self.num_particles)
			print("info :", len(self.mass), "masses")
			
		if len(self.velocity) > 0:
			if len(self.velocity) != self.num_particles:
				raise RuntimeError('Error! number of velocity ', len(self.velocity), ' is not equal to the number of particles ', self.num_particles)
			print("info :", len(self.velocity), "velocities")
			
		if len(self.bond) > 0:
			print("info :", len(self.bond), "bonds")
			
		if len(self.angle) > 0:
			print("info :", len(self.angle), "angles")

		if len(self.dihedral) > 0:
			print("info :", len(self.dihedral), "dihedrals")			

	def init_data(self):
		# data
		if not self.read_indicator['box']:				
			self.box=[]
		if not self.read_indicator['position']:				
			self.position=[]
		if not self.read_indicator['type']:			
			self.type=[]
		if not self.read_indicator['image']:			
			self.image=[]
		if not self.read_indicator['mass']:			
			self.mass=[]
		if not self.read_indicator['velocity']:			
			self.velocity=[]
		if not self.read_indicator['force']:			
			self.force=[]
		if not self.read_indicator['virial']:			
			self.virial=[]
		if not self.read_indicator['bond']:
			self.bond=[]
		if not self.read_indicator['angle']:
			self.angle=[]
		if not self.read_indicator['dihedral']:
			self.dihedral=[]
	
	# reset parameters
	def reset_params(self):
		# indicators
		self.num_particles_read=False
		self.timestep_read=False
		self.dimension_read=False		
		self.box_read=False			
		self.position_read=False
		self.type_read=False
		self.image_read=False
		self.mass_read=False
		self.velocity_read=False
		self.bond_read=False
		self.angle_read=False
		self.dihedral_read=False
