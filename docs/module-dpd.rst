Dissipative particle dynamics
=============================

DPD force
---------

Description:

    The DPD force consists of pair‐wise conservative, dissipative and random terms.

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        \vec{F}_{ij}^{C}&=&\alpha\left(1-\frac{r_{ij}}{r_{cut}}\right)\vec{e}_{ij} \\
        \vec{F}_{ij}^{D}&=&-\gamma\omega^{D}(r_{ij})(\vec{e}_{ij} \cdot \vec{v}_{ij} )\vec{e}_{ij}  \\	
        \vec{F}_{ij}^{R}&=&T\sigma\omega^{R}(r_{ij})\xi_{ij}\vec{e}_{ij} \\			
       \end{eqnarray*}

	   
    - :math:`\gamma=\sigma^{2}/2k_{B}T`
    - :math:`\omega^{D}(r_{ij})=[\omega^{R}(r_{ij})]^2=(1-r_{ij}/r_{cut})^2`	
    - :math:`\xi_{ij}` - a random number with zero mean and unit variance
    - :math:`T` - `temperature`
      - *optional*: defaults to 1.0	
    - :math:`r_{cut}` - *r_cut* (in distance units)	
      - *optional*: defaults to 1.0

    The following coefficients must be set per unique pair of particle types:
	
    - :math:`\alpha` - *alpha* (in energy units)
    - :math:`\sigma` - *sigma* (unitless)


.. py:class:: DpdForce(all_info, nlist, r_cut, temperature, rand_num)

   The constructor of a DPD interaction object.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.
   :param float temperature: The temperature.
   :param int rand_num: The seed of random number generator.   
   
.. py:class:: DpdForce(all_info, nlist, r_cut, rand_num)

   The constructor of a DPD interaction object. The default temperature is 1.0.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.
   :param int rand_num: The seed of random number generator.
   
   .. py:function:: setParams(string typei, string typej, float alpha, float sigma)
   
      specifies the DPD interaction parameters per unique pair of particle types.
	  
   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time step.
	  
   .. py:function:: setDPDVV()
   
      calls the function to enable DPDVV method (the default is GWVV).
	  
   Example::
   
      dpd = galamost.DpdForce(all_info, neighbor_list, 1.0, 12345)
      dpd.setParams('A', 'A', 25.0, 3.0)
      app.add(dpd)
	  
GWVV integration
----------------

Description:

    Integration algorithm.


    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        &v_i^0&\leftarrow v_i+ \lambda\frac{1}{m}(F_i^c \Delta t + F_i^d \Delta t + F_i^r \sqrt{\Delta t})\\ 
        &v_i&\leftarrow v_i+ \frac{1}{2}\frac{1}{m}(F_i^c \Delta t + F_i^d \Delta t + F_i^r \sqrt{\Delta t})\\
        &r_i&\leftarrow r_i+ v_i \Delta t\\
        &&Calculate \quad F_i^c\{r_j\}, F_i^d\{r_j, v_j^0\}, F_i^r\{r_j\}\\
        &v_i&\leftarrow v_i+ \frac{1}{2}\frac{1}{m}(F_i^c \Delta t + F_i^d \Delta t + F_i^r \sqrt{\Delta t})		
       \end{eqnarray*}

    - :math:`\lambda` - *lambda* (unitless)
      - *optional*: defaults to 0.65
	  
.. py:class:: DpdGwvv(AllInfo all_info, ParticleSet group)

   The constructor of a GWVV NVT thermostat for a group of DPD particles.

   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	   

   .. py:function:: setLambda(float lambda)
   
      specifies lambda.	  
	  
   Example::

      gwvv = galamost.DpdGwvv(all_info, group)
      app.add(gwvv)
	  
