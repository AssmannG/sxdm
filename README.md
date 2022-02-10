****** README file for sxdm package  data selection, merging, and offline-mode processing on Ra-cluster of Serial Synchrotron Crystallography (SSX) data collected at the SLS, Switzerland *****************

@Author: Shibom Basu, PSI
@Coauthor: Greta Assmann, PSI

The SXDM project provides scaling and merging utilities to users at SLS beamlines (PXI,PXII,PXIII) for Serial Synchrotron Crystallography (SSX). All users at PSI/SLS/SwissFEL will have access to Ra-cluster as long as they have psi-account.

The actual code is in the 'src' folder  for scaling and merging purposes.
The script sxdm.merge in the 'bin' folder is used for offline processing on Ra-cluster. run_merge.py is also necessary for the offline processing on Ra-cluster. 

The code is written in Python 3.8, relevant environments and modules neccessary for offline processingg are described in the MX software documentation from MX group at PSI (https://www.psi.ch/de/sls/pxii/mx-software-documentation).Packages used by the sxdm pipeline are cctbx, scipy, numpy, matplotlib, h5puy(etc.).  Find the list of the conda environment in conda_list.txt. 



========================== To run offline-mode processing of SSX data ======================

You MUST be either on Ra-cluster or on our beamline computing node. Follow the instructions in the MX software documentation. 
