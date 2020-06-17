Introduction
============

Coarse-grained molecular dynamics (CGMD) simulations are exceptionally important in the research field of polymers and soft matters. 
In general, CGMD targets problems typically at nano- to meso-scales that are not easily coped with using all-atom molecular 
dynamics simulations. GALAMOST was designed to provide a set of open-source and specific tools of 
employing GPUs, to accelerate CGMD simulations. 

This package is the version 5 of GALAMOST which includes MD engine **gamst**, molecular configuration generator **molgen**, and data tackler **dataTackle** etc. 
The **gamst** is purely written by Python language based on Python3 Numba compiler. The plugins **molgen** and **dataTackle** that are
written by C++ and CUDA C, need to be compiled. 

With GALAMOST, it is also possible to perform conventional CGMD, Brownian dynamics, and Dissipative Particle Dynamics simulations with various potential forms. 
The trajectories obtained in GALAMOST could be visualized and analyzed by OVITO. 

