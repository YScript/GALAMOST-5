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

import math
import gamst.snapshots.box as box_func
import numpy as np
import numba as nb
from numba import cuda
import gamst.plists.clist as clist
import time


@cuda.jit("void(int32, float32[:, :], int32[:], float32[:], float32[:], float32[:], int32[:, :], int32[:], float32[:, :, :], int32[:], int32[:, :], float32[:, :], float32, int32[:])")
def cu_neighbor_build(npa, pos, dim, box_low_boundary, inv_width, box, cell_adj, cell_size, cell_list, neighbor_size, neighbor_list, last_pos, cutoff_rsq, situation):
	i = cuda.grid(1)
	if i < npa:
		pix = pos[i][0]
		piy = pos[i][1]		
		piz = pos[i][2]
		
		dpix = pix - box_low_boundary[0]
		dpiy = piy - box_low_boundary[1]
		dpiz = piz - box_low_boundary[2]
		
		ix = nb.int32(dpix*inv_width[0])
		iy = nb.int32(dpiy*inv_width[1])
		iz = nb.int32(dpiz*inv_width[2])		
		if ix == dim[0]:
			ix = nb.int32(0)
		if iy == dim[1]:
			iy = nb.int32(0)			
		if iz == dim[2]:
			iz = nb.int32(0)

		cell_id = iz + dim[2] * (iy + ix * dim[1])
		nneigh=nb.int32(0)
		nneigh_needed=nb.int32(0)
		for ci in range(nb.int32(0), nb.int32(cell_adj.shape[1])):
			neigh_cell = cell_adj[cell_id][ci]
			for j in range(nb.int32(0), cell_size[neigh_cell]):
				pj = cell_list[neigh_cell][j]
				dp = cuda.local.array(3, dtype = nb.float32)
				dp[0] = pix - pj[0]
				dp[1] = piy - pj[1]
				dp[2] = piz - pj[2]
				box_func.cu_box_min_dis(dp, box)
				rsq = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]
				if rsq < cutoff_rsq:
					neigh_idx = nb.int32(pj[3])
					if neigh_idx != i:
						if nneigh < neighbor_list.shape[1]:
							neighbor_list[i][nneigh] = neigh_idx
						else:
							nneigh_needed = nneigh+nb.int32(1)
						nneigh += nb.int32(1)
		neighbor_size[i] = nneigh
		last_pos[i][0] = pix
		last_pos[i][1] = piy
		last_pos[i][2] = piz	
		if nneigh_needed > nb.int32(0):
			cuda.atomic.max(situation, nb.int32(1), nneigh_needed)
			
@cuda.jit("void(int32, float32[:, :], float32[:], float32[:, :], float32, int32[:])")
def cu_check_build(npa, pos, box, last_pos, rbuff_half_rsq, situation):
	i = cuda.grid(1)
	if i < npa:
		pi = pos[i]
		lpi = last_pos[i]
		dp = cuda.local.array(3, dtype = nb.float32)
		dp[0] = pi[0] - lpi[0]
		dp[1] = pi[1] - lpi[1]
		dp[2] = pi[2] - lpi[2]
		box_func.cu_box_min_dis(dp, box)
		rsq = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]
		if rsq > rbuff_half_rsq:
			situation[0] = nb.int32(1)
			
			
@cuda.jit("void(int32[:])")
def cu_zero(situation):
	i = cuda.grid(1)
	if i < 3:
		situation[i] = nb.int32(0)			
			
class nlist:
	#定义构造方法
	def __init__(self, info, rcut, rbuff):
		self.cutoff_rsq = rcut*rcut
		self.rbuff_half_rsq = (rbuff/2.0)*(rbuff/2.0)
#		print(self.rcut)
		self.info=info
		self.clist = clist.clist(info, rcut+rbuff)
		self.nmax = np.int32(8)
		self.neighbor_list = np.zeros([self.info.npa, self.nmax], dtype = np.int32)
		self.block_size = 64
		self.last_pos = np.zeros([self.info.npa, 3], dtype = np.float32)		
		self.neighbor_size = np.zeros([self.info.npa], dtype = np.int32)

		self.situation = np.zeros(3, dtype = np.int32)
		self.d_situation = cuda.to_device(self.situation)
		
#		device arrays
		self.d_neighbor_list = cuda.to_device(self.neighbor_list)
		self.d_last_pos = cuda.to_device(self.last_pos)
		self.d_neighbor_size = cuda.to_device(self.neighbor_size)


#		jump check parameters
		self.call = np.int32(0)
		self.p_call = np.int32(0)
		self.update_period = np.int32(1)
		self.p_update = np.int32(0)
		self.update_offset = np.int32(0)
		
		
	def calculate(self, build_list):
		if self.if_jump():
			return

		nblocks = math.ceil(self.info.npa / self.block_size)
		# start = time.time()
		# self.d_situation = cuda.to_device(self.situation_zero)
		cu_zero[1,32](self.d_situation)
		cu_check_build[nblocks, self.block_size](self.info.npa, self.info.d_pos, self.info.d_box, self.d_last_pos, self.rbuff_half_rsq, self.d_situation)
		self.situation = self.d_situation.copy_to_host()

		if self.situation[0]==np.int32(1):
			build_list = True
			self.update_offset = self.call
			self.p_update += np.int32(1)
		# end1 = time.time()			
		# cuda.synchronize()
		# print("check build",end1 - start)
		while build_list:
			# end1 = time.time()				
			self.clist.calculate()
			# self.d_situation = cuda.to_device(self.situation_zero)
			# cuda.synchronize()
			# end2 = time.time()	
			# print("cell list",end2 - end1)
			cu_zero[1,32](self.d_situation)
			cu_neighbor_build[nblocks, self.block_size](self.info.npa, self.info.d_pos, self.clist.d_dim, self.clist.d_box_low_boundary, self.clist.d_inv_width, self.info.d_box, self.clist.d_cell_adj, 
			                                            self.clist.d_cell_size, self.clist.d_cell_list, self.d_neighbor_size, self.d_neighbor_list, self.d_last_pos, self.cutoff_rsq, self.d_situation)
			self.situation = self.d_situation.copy_to_host()
			# cuda.synchronize()
			# end3 = time.time()
			# print("neigbor list",end3 - end2)
			if self.situation[1]>self.nmax:
				self.nmax = self.situation[1]
				self.nmax = self.nmax + np.int32(8) - (self.nmax & np.int32(7))
				self.neighbor_list = np.zeros([self.info.npa, self.nmax], dtype = np.int32)
				self.d_neighbor_list = cuda.to_device(self.neighbor_list)
			else:
				build_list = False
		# nlist = self.d_neighbor_list.copy_to_host()
		# nlist_size = self.d_neighbor_size.copy_to_host()
		# for i in range(0, nlist.shape[0]):
			# if nlist_size[i]>0:
				# print(i, nlist_size[i], nlist[i])
	def if_jump(self):
		self.call += np.int32(1)
		self.p_call += np.int32(1)
		if self.call%np.int32(1000)==np.int32(0):
			self.update_period = np.int32(self.p_call/self.p_update)
			self.p_update = np.int32(0)
			self.p_call = np.int32(0)
			# print(self.update_period)
		if self.call - self.update_offset < np.int32(np.float32(self.update_period)*np.float32(0.8)):
			return True
		return False
		

		
		
		
		
		
		
