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
from numba import cuda
import numba as nb

# device functions
def cu_nonbonded(force_name):
	@cuda.jit(device=True)
	def _lj(rsq, param, fp):
		lj1 = param[0]
		lj2 = param[1]
		rcutsq = param[2]
		if rsq<rcutsq:
			r2inv = nb.float32(1.0)/rsq
			r6inv = r2inv * r2inv * r2inv
			f = r2inv * r6inv * (nb.float32(12.0) * lj1  * r6inv - nb.float32(6.0) * lj2)	
			p = r6inv * (lj1 * r6inv - lj2)
			fp[0]=f
			fp[1]=p
			
	@cuda.jit(device=True)
	def _harmonic(rsq, param, fp):
		alpha = param[0]
		r_cutINV = param[1]
		rcutsq = param[2]
		if rsq<rcutsq:
			rinv = nb.float32(1.0)/math.sqrt(rsq)
			omega = rinv - r_cutINV
			f = alpha*omega
			p = nb.float32(0.5)*alpha*omega*omega*rsq
			fp[0]=f
			fp[1]=p

	if force_name=="lj":
		return _lj
	elif force_name=="harmonic":
		return _harmonic

# host functions
def nonbonded(force_name):
	def _lj(rsq, param, fp):
		lj1 = param[0]
		lj2 = param[1]
		rcutsq = param[2]
		if rsq<rcutsq:
			r2inv = 1.0/rsq;
			r6inv = r2inv * r2inv * r2inv;
			f = r2inv * r6inv * (12.0 * lj1  * r6inv - 6.0 * lj2)	
			p = r6inv * (lj1 * r6inv - lj2)
			fp[0]=f
			fp[1]=p
	if force_name=="lj":
		return _lj

		
