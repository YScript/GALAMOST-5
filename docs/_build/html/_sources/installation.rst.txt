Installation
============

Source code
-----------

The entire GALAMOST package is a Free Software under the GNU General Public License. 
The **gamst** engine, **molgen** plugin, and **galaTackle** plugin could be separately installed. 
The installation instructions for **gamst**, **molgen**, and **galaTackle** are as follows.

1. Installation for **gamst**::
	
	Requrements:
	
	1. Python3 including numba, numpy, cupy, and pybind11 packages
	2. NVIDIA CUDA Toolkit >= 7.0
	
	Installation:
	
	python3 setup.py install

2. Installation for **molgen**::

    python3 setup.py build
    python3 setup.py install
	
3. Installation for **galaTackle**::

    sh compile.sh	
