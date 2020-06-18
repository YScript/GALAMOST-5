Usage
=====

With a prepared script, you could run gamst MD engine for abstaining trajectory.

   Examples::
   
      python3 yourscript.gala --gpu=0 >a.log&
	  
Where you could specify the GPU id (default value is 0) with the ``--gpu=`` option and output the screen information into ``a.log`` file.

For version 4, you could use multiple GPU for parallel computation by following commands:	 

   Examples::
   
      mpirun -n 4 python yourscript.gala --gpu=0,1,2,3 >a.log& 

Here is an example of script for DPD simulation. 

Firstly, importing the gamst module installed as a package of python3 and reading system information by :py:class:`snapshot.read` from a mst file 

   Examples::

      import gamst
      mst = gamst.snapshot.read("AB.mst")
	  
After that, we need to build up an application by :py:class:`application.dynamics` which will call defined and added objects.

   Examples::
   
      app = gamst.application.dynamics(info=mst, dt=0.04)

Further, we should define objects by the classes of gamst and pass them to the application, such as the following example: DPD force :py:class:`force.dpd`
, NVT thermostat with GWVV algorithm :py:class:`integration.gwvv`, and th dump of system collective information :py:class:`dump.data`.
      
   Examples::
  	  
      fn = gamst.force.dpd(info=mst, rcut=1.0)
      fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
      fn.setParams(type_i="A", type_j="B", alpha=40.0, sigma=3.0)
      fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
      app.add(fn)
      
      
      gw = gamst.integration.gwvv(info=mst, group='all')
      app.add(gw)
      
      di = gamst.dump.data(info=mst, group='all', file='data.log', period=500)
      app.add(di)

	  
Finally, running the simulation with the number of time steps.
      
   Examples::
  	  
      app.run(10000)


