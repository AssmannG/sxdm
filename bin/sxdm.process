#!/bin/bash

input="/sls/X06DA/data/e11206/Data10/shibom/20170528_TD1/TB2/TB2_4"

#output="/sls/X06DA/data/e11206/Data10/shibom/20170528_TD1/TB2"
output=`echo $PWD`

spgroup="4"
#unit_cell="85 398 150 90 90 90"
unit_cell="73.82  91.2  82.99  90.000  97.32  90.000"

friedel=FALSE
res_cutoff="2.4"

#refs="/sls/X06SA/data/e10019/Data10/20170313_sercat/proc2/IMISX_PepT_1/minisets_20170313_163840/xtal_0/XDS_ASCII.HKL"
refs=" "               #uncomment this option if you do not want to use reference dataset
strong_pix="3.0"
min_pix_spot="1"
tot_degree="360"
bl="PXIII"
isa="4.5"
method="native-sad"

merge_paths="$output/proc"
#"/sls/X06SA/data/e10010/Data10/PXI/20170325/Kathrin/KJ0765_proc1/KJ0765 \
#/sls/X06SA/data/e10010/Data10/PXI/20170325/Kathrin/KJ0765_proc1/KJ0766 \
#/sls/X06SA/data/e10010/Data10/PXI/20170325/Kathrin/KJ0765_proc1/KJ0767 \
#/sls/X06SA/data/e10010/Data10/PXI/20170325/Kathrin/KJ0765_proc1/KJ0768 \
#/sls/X06SA/data/e10010/Data10/PXI/20170325/Kathrin/KJ0765_proc1/KJ0769 \
#/sls/X06SA/data/e10010/Data10/PXI/20170325/Kathrin/KJ0765_proc1/KJ0770"

#salloc -N4

#=========== do not edit below part ============
if [[ -n "${BEAMLINE_XNAME}" ]]; then
  source /opt/gfa/python 2.7
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'/exchange/mx/daplus/imisx_processing/cctbx/build/lib:/exchange/mx/daplus/imisx_processing/cctbx/base/lib:/exchange/mx/daplus/conda_envs/adp_copy/lib:/exchange/mx/daplus/conda_envs/adp_copy/lib/python2.7/site-packages'
  export LIBTBX_BUILD='/exchange/mx/daplus/imisx_processing/cctbx/build'
  export PYTHONPATH='/exchange/mx/daplus/imisx_processing/cctbx/build/lib:/exchange/mx/daplus/imisx_processing/cctbx/modules/cctbx_project:/exchange/mx/daplus/imisx_processing/cctbx/modules/cctbx_project/boost_adaptbx:/exchange/mx/daplus/conda_envs/adp_copy/lib'
  export PYTHONPATH=$PYTHONPATH:'/exchange/mx/daplus/git/sxdm/src'
  python /exchange/mx/daplus/git/sxdm/xds/sxdp_v1.py --image_path $input --output_dir $output \
           --BeamID $bl --total_degree $tot_degree --method $method \
           --SG_num $spgroup --friedel $friedel --highres $res_cutoff \
           --cell "$unit_cell" --refs "$refs" --strong_pixel $strong_pix \
           --min_pix_spot $min_pix_spot  --ISa_cutoff $isa --merge_paths "$merge_paths"
else
    echo "using ra-cluster"
    echo "cctbx and sxdm have to be added to LD_LIBRARY_PATH and PYTHONPATH"
    salloc -N2 -p day -I2
    python $script/sxdp_v1.py --image_path $input --output_dir $output \
             --BeamID $bl --total_degree $tot_degree --method $method \
             --SG_num $spgroup --friedel $friedel --highres $res_cutoff \
             --cell "$unit_cell" --refs "$refs" --strong_pixel $strong_pix \
             --min_pix_spot $min_pix_spot  --ISa_cutoff $isa --merge_paths "$merge_paths"
fi
