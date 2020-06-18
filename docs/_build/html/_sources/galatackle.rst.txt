galaTackle
==========

Usage
-----

The **galaTackle** is a plugin of GALAMOST to analyze some important properties by reading the generated configuration files. 
The **galaTackle** can be compiled and installed by "sh compile.sh". You can use this tool to analyze one or more files at a time.

   Examples for running::
   
      galaTackle particle.mst
      galaTackle particle0.mst particle1.mst
      galaTackle *.mst
	  
After pressing enter, a menu of function options will be listed. You can choose one or more functions by 
the number indexes separated by blank. The parameters for a function can be input after ":" and seperated by "|". 

   Such as::
   
      4:npot=1000 3:gpu=0|rmax=2.0

The **galaTackle** also can be called and passed parameters in a Shell script.

   Such as::
   
      echo -e "18:gpu=1|rmax=3.0" | galaTackle particles.*.mst

The **galaTackle** plugin supports the configuration files with MST and DCD formats. 
The MST files can be tackled independently. However, a DCD trajectory file has to be tackled along with 
a MST file for particles attributes and topological information.

   Examples::

      galaTackle particle.mst trajectory.dcd

For the help about the parameters, you could input "function number:h" after the function list shown.
 	  
   Examples::

      14:h	  
	  
Functions
---------

   Function list::
   
      -------------------------------------------------------------
      1  Rg^2              2  Ed^2             3  RDF              
      4  bond_distri       5  angle_distri     6  dihedral_distri  
      7  stress tensor     8  density          9  unwrapping       
      10 MSD               11 RDF-CM           12 MSD-CM           
      13 ents              14 strfac           15 domain size      
      16 dynamic strfac    17 config check     18 RDF between types
      19 MST conversion 
      -------------------------------------------------------------

1  Rg^2:
^^^^^^^^

   Description:
      To calculate the mean square of gyration radius and output result to ``rg2.log``.

    .. math::
        :nowrap:

        \begin{eqnarray*}
		R{_{g}^{2}}&=&\frac{1}{N}\sum\limits_{i=1}^{N}{{{({{{\vec{R}}}_{i}}-{{{\vec{R}}}_{cm}})}^{2}}} \\
		{{\vec{R}}_{cm}}&=&\frac{1}{N}\sum\limits_{i=1}^{N}{{{{\vec{R}}}_{i}}}
        \end{eqnarray*}

    Coefficients:

    - :math:`{\vec{R}}_{i}` - monomer position vector 

2  Ed^2:	  
^^^^^^^^
   
   Description:
      To calculate the mean square of end to end distance and output result to ``ed2.log``.
	  
    .. math::
        :nowrap:

        \begin{eqnarray*}
		E{_{d}^{2}}={( {{{\vec{R}}}_{0}}-{{{\vec{R}}}_{N-1}} )}^{2}
        \end{eqnarray*}

    Coefficients:

    - :math:`{\vec{R}}_{i}` - monomer position vector 	  
	  
3  RDF:	  
^^^^^^^
   
   Description:
      To calculate the radial distribution function of all particles and output result to ``*.rdf`` and ``rdf.log``.
	  
   Parameters:
      :maxbin=100|gpu=0|rmax=Lx/2|bondex=false|angleex=false|molex=false
	  
4  bond_distri:	  
^^^^^^^^^^^^^^^

   Description:
      To calculate the distribution of bond lengths and output result to ``bond_distr.log``.

    .. math::
        :nowrap:

        \begin{eqnarray*}
		bond\_distri(i \cdot dr)=N(i)/(N \cdot dr)
        \end{eqnarray*}

    Coefficients:

    - :math:`dr` - the space of bond length `L/(2npot)`, where `L` is the box size
    - :math:`N(i)` - the number of bonds in the range of `idr < r < (i+1)dr`, where `i` is an integer
    - :math:`N` - the total number of bonds		
	  
   Parameters:
      :npot=2001

5  angle_distri:	  
^^^^^^^^^^^^^^^^
   
   Description:
      To calculate the distribution of angle degrees and output result to ``angle_distr.log``.
	  
    .. math::
        :nowrap:

        \begin{eqnarray*}
		angle\_distri(i \cdot da)=N(i)/(N \cdot da)
        \end{eqnarray*}

    Coefficients:

    - :math:`da` - the space of angle radian `pi/npot`
    - :math:`N(i)` - the number of angles in the range of `ida < angle < (i+1)da`, where `i` is an integer	  
    - :math:`N` - the total number of angles	
	
   Parameters:
      :npot=2001
	  
6  dihedral_distri:	  
^^^^^^^^^^^^^^^^^^^
   
   Description:
      To calculate the distribution of dihedral degrees and output result to ``dihedral_distr.log``.

    .. math::
        :nowrap:

        \begin{eqnarray*}
		dihedral\_distri(i \cdot da)=N(i)/(N \cdot da)
        \end{eqnarray*}

    Coefficients:

    - :math:`da` - the space of dihedral angle radian `2pi/npot`
    - :math:`N(i)` - the number of dihedrals in the range of `ida < dihedral angle < (i+1)da`, where `i` is an integer
    - :math:`N` - the total number of dihedrals		
	  
   Parameters:
      :npot=2001
	  
7  stress tensor:	  
^^^^^^^^^^^^^^^^^
   
   Description:
      To calculate the stress tensor by inputing the parameters of force calculation and output result to ``stress_tensor.log``.
	  
   Parameters:
      :bondex=true|bodyex=true|diameter=true 

8  density:	  
^^^^^^^^^^^
   
   Description:
      To calculate the real density (g/cm^3) with basic units [amu] and [nm] and output result to ``density.log``.
	  
9  unwrapping:	  
^^^^^^^^^^^^^^
   
   Description:
      To unwrap or shift molecules by changing the image information
	  
   Parameters:
      :unwrap_molecule=true|label_free_particle=particle type|molecule_center_in_box=false| shiftx=0.0|shifty=0.0|shiftz=0.0|remove_image=false|convert_constraint_to_bond=false| remove_bond_cross_box=false

	  
10 MSD:	  
^^^^^^^
   
   Description:
      To compute mean square displacement of all particles and output result to ``msd.log``.

   Parameters:
      :direction=XYZ (candidates are X,Y,Z,XY,YZ,XZ,XYZ) 
	  
11 RDF-CM:	  
^^^^^^^^^^
   
   Description:
      To calculate the radial distribution function of the mass center of molecules and output result to ``rdf_cm.log``.
	  
   Parameters:
      :maxbin=100|gpu=0|rmax=Lx/2	  
	  
12 MSD-CM:	  
^^^^^^^^^^
   
   Description:
      To compute mean square displacement of the mass center of molecules and output result to ``msd_cm.log``.
	  
   Parameters:
      :direction=XYZ (candidates are X,Y,Z,XY,YZ,XZ,XYZ)  
	  
13 ents:	  
^^^^^^^^
   
   Description:
      To analyze the entanglements of polymers and output result to ``ents.log``.
	  
14 strfac:	  
^^^^^^^^^^
   
   Description:
      To calculate the structure factor of particles and output result to ``*.strf`` and ``strf.log``.
	  
   Parameters:
      :qmax=160pi/Lmin|gpu=0|deltaq=2pi/Lmin|direction=XYZ|2D=false

15 domain size:	  
^^^^^^^^^^^^^^^
   
   Description:
      To calculate the domain size of components in mixtures and output result to ``domsize.log``.
	  
   Parameters:
      :kmax=20|qc=0.4690|gpu=0

16 dynamic strfac:	  
^^^^^^^^^^^^^^^^^^
   
   Description:
      Dynamic structure factor (incoherent intermediate) measures the decorrelation of the positions 
      of individual monomers with the time on length scale :math:`1/q`, where :math:`q=2\pi\sqrt{x^2+y^2+z^2}/L`, and :math:`L` is cubic box length. 
      :math:`\mbox{kmax}` limits the space in which the :math:`q` with possible combinations of x, y, z will be generated.

      The results will be output to ``dstrf.log``.

   Parameters:
      :kmax=int(L)|q=7.0
	  
   `Maintainer`: Shu-Jia Li
	  
17 config check:	  
^^^^^^^^^^^^^^^^
   
   Description:
      To check the configuration including the minimum distance of particles, and the Maximum and minimum length of bonds, etc. and output result to ``config_check.log``.

   Parameters:
      :bondex=true|angleex=true|dihedralex=true|bodyex=true|rcut=2.0


18 RDF between types:	  
^^^^^^^^^^^^^^^^^^^^^
   
   Description:
      To compute the radial distribution function between types and output result to ``*.type.rdf`` and ``rdf_by_type.log``.
	  
   Parameters:
      :maxbin=100|gpu=0|rmax=Lx/2|bondex=false|angleex=false|molex=false


19 MST conversion:	  
^^^^^^^^^^^^^^^^^^
   
   Description:
      Converting MST files into other formats
	  
   Parameters:
      :lammps=false|gromacs=false|xml=false
  
	  
 	  