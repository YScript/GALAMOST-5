Non-bonded interactions
=======================
   
**Overview**

The gamst MD engine provides a few of functions for non-bonded interactions.
However, it supports well self-defined analytical functions via writting codes of 
device function in script.

=================   ===========================
:ref:`non-bonded`   :py:class:`force.nonbonded`
=================   ===========================


.. _non-bonded:

Non-bonded functions
--------------------

Description:

   The function of non-bonded interactions could be either one called from non-bonded interaction function libary, or a self-defined device function.
   Non-bonded interaction function libary contains Lennard-Jones function named as 'lj' and harmonic function named as 'harmonic'.
   
   Lennard-Jones function
    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)  = & 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                          \alpha \left( \frac{\sigma}{r} \right)^{6} \right] & r < r_{\mathrm{cut}} \\
                            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    The following coefficients must be set per unique pair of particle types:

    - :math:`\epsilon` - *epsilon* (in energy units)
    - :math:`\sigma` - *sigma* (in distance units)
    - :math:`\alpha` - *alpha* (unitless)
    - :math:`r_{\mathrm{cut}}` - *r_cut* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command   
   

.. py:class:: force.nonbonded(info, rcut, func)

   The constructor of non-bonded interaction calculation object.
	  
   :param info: system information.
   :param rcut: cut-off radius of interactions.
   :param func: function name.

   .. py:function:: setParams(type_i, type_j, param)
 
      specifies interaction parameters with type_i, type_j, a list of parameters.  
   
   Example::
   
      fn = gamst.force.nonbonded(info=mst, rcut=3.0, func='lj')
      fn.setParams(type_i="a", type_j="a", param=[1.0, 1.0, 1.0, 3.0])
      app.add(fn)



.. _self-defined-function:

Self-defined functions
----------------------

Description:

   The device function for non-bonded interactions could be written in script and conveyed 
   to kernal funcitons for calculation.
   
   Non-bonded interactions with potential form :math:`p(r)`

   * p = :math:`p(r)`
   * f = :math:`-(\triangle p(r)/\triangle r)(1/r)`

   Function code template::
   
   		@cuda.jit(device=True)
		def func(rsq, param, fp):
			rcut = param[0]
			p1 = param[1]
			p2 = param[2]
			p3 = param[3]
			...
			if rsq<rcut*rcut:
				calculation codes
				...
				fp[0]=f
				fp[1]=p
				
		fn = gamst.force.nonbonded(info, rcut, func)
		fn.setParams(type_i, type_j, param=[rcut, p1, p2, p3, ...])
		....
		app.add(fn)		
   
   Example::
   
		from numba import cuda
		import numba as nb
		
		@cuda.jit(device=True)
		def lj(rsq, param, fp):
			epsilon = param[0]
			sigma = param[1]
			alpha = param[2]
			rcut = param[3]
			if rsq<rcut*rcut:
				sigma2 = sigma*sigma
				r2inv = sigma2/rsq;
				r6inv = r2inv * r2inv * r2inv;
				f = nb.float32(4.0) * r2inv * r6inv * (nb.float32(12.0) 
				    * r6inv - nb.float32(6.0) * alpha)/sigma2	
				p = nb.float32(4.0) * r6inv * ( r6inv - nb.float32(1.0))
				fp[0]=f
				fp[1]=p
				
		fn = gamst.force.nonbonded(info=mst, rcut=3.0, func=lj)
		fn.setParams(type_i="a", type_j="a", param=[1.0, 1.0, 1.0, 3.0])
		app.add(fn)	
 
