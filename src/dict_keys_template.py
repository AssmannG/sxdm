import os,sys
import logging

logger = logging.getLogger('Scale&Merge')

output = {
'status': False,
'merge_output_path': None, #string type
'xtals_expected': None, # integer type showing totla number of xtals from hklpaths.
'xtals_found': None, #integer type based on found ASCII.HKL files
'xtals_after_idx_check': None, #integer type
'xtals_after_isa': None, #integer type
'xtals_after_cell': None, #integer type
'xtals_after_pCC': None, #integer type
'xtals_rejected': None, #integer type
'rejected_hkls': [], #list of hkl paths
'space_group': None, #space group number
'unit-cell': {'a': None, 'b': None, 'c': None, 'al': None, 'be': None, 'ga': None},
'nSAD_xscale_stats': [], #list of dictionaries, where each row of stat table is a dictionary.
'no_selection': [], #list of dictionaries, only for SX type expt.
#'Rmeas_selection': [], #only for SX type expt
'ISa_selection': [], #only for SX and nSAD type expt
'cell_selection': [], #only for SX type expt
'pCC_selection':[], #only for SX type expt
'iso-cluster': [], # stats from cluster 1 for nSAD and SX type expts - both
'clusters': None, # no. of Best clusters from iso-clustering
'anisotropicity': " ", #tuple of ('a:2.0', b:1.0, c:3.0)
'multiplicity':[], #list of lists to plot multiplicity of each reflection as scatter
'cell_array': [], #numpy array of shape(datasets,6)
'histogram': {}, #cell_array converted to dictionary for highcharts
'hclust_matrix': [], #numpy array of shape(datasets-1,4), storing linkage matrix from hierarchical clustering
'cell_dendrogram': {}, #stores dendrogram dictionary, use 'ivl' key-values for labels when make plot
'cc_dendrogram':{},
'cell_n_clusters':None,
'cc_n_clusters': None,
'info': [], #list of string types data to store all error messages.
}

def push_val2dict(key, value):
    if key in output.keys():
        output[key] = value
        output['status'] = True
    else:
        err = 'KeyError: the key does not exist'
        logger.info('Error:{}'.format(err))
        output['status'] = False
        return output
    return output

def add_new_key(key, value):
    try:
        output[key].append(value)
    except (KeyError, ValueError, TypeError) as err:
        output[key] = value
    return output
