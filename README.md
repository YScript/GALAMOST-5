# GALAMOST: GPU-Accelerated LArge-scale MOlecular Simulation Toolkit (GALAMOST) 
## Version 5, Python3 Numba version 
## Copyright: G.U.G.D.(Geek Union of GALAMOST Developers)

Coarse-grained molecular dynamics (CGMD) simulations are exceptionally important in the research field of soft matter. 
In general, CGMD targets problems typically at nano- to meso-scales that are not easily coped with using all-atom molecular 
dynamics simulations. `GALAMOST` was designed to provide a set of open-source and specific tools of employing GPUs, to accelerate CGMD simulations. 

These tools are, for example:

1. Modelling “polymerizations” in a stochastic way in CGMD simulations. The topological connections between beads could be changed according to pre-defined probability. 
This tool can be used in the studies of chain-growth polymerization, reversible reaction, exchange reaction, and so on.

2. Using anisotropic particle models to describe soft particles as well as rigid ellipsoidal particles with/without patches by combining harmonic repulsive 
potential/Gay-Berne potential with anisotropic “patchy” potentials.

3. Using hybrid particle-field technique to accelerate CGMD simulations (G. Milano and T. Kawakatsu, J. Chem. Phys. 130, 214106, 2009). 
This method could largely speed up some slowly evolving processes in CGMD simulations, such as microphase separation and self-assembly of polymeric systems.

4. Reading in numerical potentials derived from iterative Boltzmann inversion, inverse Monte Carlo, or other structure-based bottom-up coarse-graining methods 
(Mirzoev et al., Comput. Phys. Commun. 237, 263, 2019), and applying the potentials in CGMD simulations.

In GALAMOST, it is also possible to perform conventional CGMD, Brownian dynamics, and dissipative particle dynamics simulations with various potential forms. 
The trajectories obtained in GALAMOST could be visualized and analyzed by OVITO. 

## Installation
gamst:

    python3 setup.py install

molgen:

    python3 setup.py build

    python3 setup.py install

galaTackle:

    sh compile.sh

### Requrements:
1. Python3 including numba, numpy, cupy, and pybind11 packages
2. NVIDIA CUDA Toolkit >= 7.0

## Citation

```
To cite GALAMOST in publications use:
 
Y.-L. Zhu, H. Liu, Z.-W. Li, H.-J. Qian, G. Milano, Z.-Y. Lu, J. Comput. Chem., 34, 2197, 2013.

A BibTeX entry for LaTeX users is
@ARTICLE{Zhu2013,
   author = {Y.-L. Zhu and H. Liu and Z.-W. Li and H.-J. Qian and G. Milano and Z.-Y. Lu},
   title = {GALAMOST: GPU-accelerated large-scale molecular simulation toolkit},
   journal = {J. Comput. Chem.},
   volume = {34},
   pages = {2197},
   year = {2013}
}
``` 


## Documentation

Online manual could be read here: https://galamost.readthedocs.io/en/latest/.
Tutorials written by jupyter notebook are given here: https://nbviewer.jupyter.org/github/pigwarrior/GALAMOST-5/tree/master/tutorials/index.ipynb.
More examples could be found here: https://github.com/pigwarrior/GALAMOST-5/tree/master/examples.

## Example: DPD simulation of diblock copolymer

1 First step: generate configuration

```
import molgen

mol1=molgen.Molecule(10)#particle number
mol1.setParticleTypes("A,A,A,A,A,B,B,B,B,B")#type
mol1.setTopology("0-1,1-2,2-3,3-4,4-5,5-6,6-7,7-8,8-9")#topology
mol1.setBondLength(0.75)#bond length
mol1.setMass(1.0)#mass


gen=molgen.Generators(20,20,20) # box size in x, y, and z direction
gen.addMolecule(mol1,2400)#molecule, the number of molecules
gen.outPutMST("A5B5") #file name
```

2 Second step: run simulation

```
import gamst
	
mst = gamst.snapshot.read("A5B5.mst")
app = gamst.application.dynamics(info=mst, dt=0.04)

fn = gamst.force.dpd(info=mst, rcut=1.0)
fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
fn.setParams(type_i="A", type_j="B", alpha=40.0, sigma=3.0)
app.add(fn)

fb = gamst.force.bond(info=mst, func='harmonic')
fb.setParams(bond_type = 'A-A', param=[4.0, 0.0])# param=[k, r0]
fb.setParams(bond_type = 'A-B', param=[4.0, 0.0])# param=[k, r0]
fb.setParams(bond_type = 'B-B', param=[4.0, 0.0])# param=[k, r0]
app.add(fb)

inn = gamst.integration.gwvv(info=mst, group='all')
app.add(inn)

dm = gamst.dump.mst(info=mst, group=['A', 'B'], file='p.mst', period=10000)
app.add(dm)

di = gamst.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(di)

app.run(500000)
```

## Contributing

We welcome contributions to GALAMOST. Whether it is reporting a bug, starting a discussion by asking a question, or proposing/requesting a new feature, 
please go by creating a new issue here (https://github.com/pigwarrior/GALAMOST-5/issues/) or writing an email to the author Dr. You-Liang Zhu (Email: ylzhu@galamost.com) 
so that we can talk about it. Please note that this project is released with a Contributor Code of Conduct. 
By participating in this project you agree to abide by its terms.


