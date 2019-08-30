cpex: Crystal Plasticity Data Extraction
===============================================

What is cpex?
-------------


Example Usage
-------------


Requirements
------------

PyXe is built on Python’s scientific stack (numpy, scipy, matplotlib). Additionally, the h5py package is required for the manipulation and management of the NeXus data files. Testing and development were carried out using the Anaconda (v 2019.03) package manager, which built with the following versions:

-	Python: version 3.7.3
-	numpy: version 1.16.2
-	scipy: version 1.2.1
-	matplotlib: version 3.0.3
-	h5py: version 2.8.0

Backward compatability to Python 3.5 is likely but not guaranteed. Monochromatic XRD caking/azimuthal integration within pyXe relies on pyFAI (and fabIO), which is a software package developed at the ESRF, designed to reduce SAXS, WAXS and XRPD images recorded by area detectors to 1D plots or 2D patterns. This caking functionality is not currently under development within pyXe and recent developments within pyFAI may have broken this functionality. While this may be fixed in the future we currently advise that azimuthal integration be carried out as a pre-processing step at the beamline (using pyFAI at ESRF or DAWN at DLS); the pyXe monochromatic post-processing analysis platform should be flexible enough to deal with most data inputs (although interface modules will likely be required outside of the Diamond Light Source).

Installation
------------

Install from the distribution using the setup.py script. The source is stored in the GitHub repo, which can be browsed at:

https://github.com/casimp/cpex

Simply download and unpack, then navigate to the download directory and run the following from the command-line:

```
pip install .
```