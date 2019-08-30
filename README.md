cpex: Crystal Plasticity Data Extraction
===============================================

[![Build Status](https://travis-ci.org/casimp/cpex.svg?branch=master)](https://travis-ci.org/casimp/cpex) 
[![Build status](https://ci.appveyor.com/api/projects/status/9cc2aej45li1pm97?svg=true)](https://ci.appveyor.com/project/casimp/cpex/branch/master)

What is cpex?
-------------


Example Usage
-------------


Requirements
------------

cpex is built on Python’s scientific stack (numpy, scipy, matplotlib). Additionally, the h5py package is required for the manipulation and management of the NeXus data files. Testing and development were carried out using the Anaconda (v 2019.03) package manager, which built with the following versions:

-	Python: version 3.7.3
-	numpy: version 1.16.2
-	scipy: version 1.2.1
-	matplotlib: version 3.0.3
-	h5py: version 2.8.0

Backward compatability to Python 3.5 is likely but not guaranteed. 

Installation
------------

Install from the distribution using the setup.py script. The source is stored in the GitHub repo, which can be browsed at:

https://github.com/casimp/cpex

Simply download and unpack, then navigate to the download directory and run the following from the command-line:

```
pip install .
```