****** README file for sxdm package of native-SAD phasing, data selection, merging, and offline-mode processing of Serial Synchrotron Crystallography (SSX) data collected at the SLS, Switzerland *****************

@Author: Shibom Basu, PSI

SXDM project aims at providing scaling and merging utilities to users at SLS beamlines (PXI,PXII,PXIII)
These codes are well-suited for Ra-cluster at PSI because some parts (SHELXD-grid search) require slurm queuing system.
All users at PSI/SLS/SwissFEL will have access to Ra-cluster as long as they have psi-account.

codes under 'src' folder are for scaling and merging purposes.
codes under 'guis' folder are for native-SAD phasing - SHELXD (brute-force grid search), model building
codes under 'xds' folder are for offline-mode processing of SSX data collected at the SLS. These codes eventually calls merging/scaling codes from 'src' folder to provide automated merging after processing of individual "minisets"

For general users, I am making wrapper scripts, which can easily called via command line to excute certain task

Potentially, they can be used in our beamline computing nodes as well, except the grid-search.
All scripts are written based on python 2.7

source /opt/gfa/python 2.7 for our MX beamlines computing nodes
Merging utility under 'src' folder is more independent (Ra, beamline nodes, even on personal computer)
and can be run locally even from one's home computer.

You can potentially add the folder into your PYTHONPATH and import as individual modules.
Then, you can use all methods/functions in your python code.

Dependencies for the package:

Python 2.7, with scipy, numpy, h5py, matplotlib.pyplot, pyQt4, and cctbx library.

-------------------------------------------------INSTALLATION-------------------------------------------------
I recommend to use conda (Anaconda to install all the required packages/libs to run the program)
usual conda command found: https://conda.io/docs/_downloads/conda-cheatsheet.pdf
conda search PACKAGENAME
conda install PACKAGENAME
conda list <Check what is available or installed in your computer>

cctbx installation:
cd local/installation/folder
git clone sxdm-package
mkdir cctbx
cd cctbx
copy bootstrap.py script from the link: https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
python bootstrap.py --builder=cctbx --nproc=4 --with-python=/path/to/your/conda-environment/bin/python

>>>bootstrap.py will git pull all the necessary libs and install dependencies. If it fails at certain dependant package, use your conda to install it.
>>>E.g. conda install six; conda install curl 

following environment has to be set/added in shell (.bashrc) before importing cctbx libs

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'/path/to/cctbx/build/lib:/path/to/cctbx/base/lib:/path/to/anaconda-environment/lib:/path/to/anaconda/lib/python2.7/site-packages'
export LIBTBX_BUILD='/path/to/cctbx/build'
export PYTHONPATH=$PYTHONPATH:'/path/to/cctbx/build/lib:/path/to/cctbx/modules/cctbx_project:/path/to/cctbx/modules/cctbx_project/boost_adaptbx:/path/to/conda-environment/lib'
export PYTHONPATH=$PYTHONPATH:'/path/to/sxdm/src:/path/to/sxdm/xds:/path/to/sxdm/guis'

cd sxdm
make install PREFIX=/path/to/installation/folder

Nothing needs to be installed if you are using beamline nodes or Ra cluster. you can call scripts directly.
Above installation procedure is to be followed only if you want to use locally (DONOT process SSX data collected at SLS on your local computer!!!)

========================== To run offline-mode processing of SSX data ======================

You MUST be either on Ra-cluster or on our beamline computing node.
copy the bash script called "run_proc.sh" located under merging-utility/xds folder to your working directory.
Then, edit run_proc.sh with input_paths (where your minisets data are), reference dataset, unit_cell, space_group_number, resolution cutoff, isa-cutoff. Then run it.

On Ra cluster:
<allocate computing nodes>
salloc -p day -N 2 -I2
<Run the script>
./run_proc.sh
<upon being finished, type>
scancel -u $user
