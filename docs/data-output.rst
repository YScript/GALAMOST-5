Data output
===========

Collective information
----------------------

.. py:class:: dump.data(info, group, file, period)

   Constructor of an information dump object for a group of particles.
   
   :param info: System information.	
   :param group: A group of particles.	   
   :param filename: Output file name.  
   :param period: The period to output data.
  
	  
   Example::
   
		dd = gamst.dump.data(info=mst, group=['a'], file='data.log', period=100)
		app.add(dd)

