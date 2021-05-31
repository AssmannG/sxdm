
__author__ = ["S. Basu"]
__license__ = "M.I.T"
__date__ = "02/02/2018"
__refactordate__ = "25/05/2021"

import sys, os
import numpy as np
import logging
from ascii import ASCII
from cellprobe import Cell
import multiprocessing as mp
import scipy.cluster.hierarchy as sch
from abstract import Abstract

logger = logging.getLogger('sxdm')

try:
    from cctbx.array_family import flex
    from cctbx import crystal
    from cctbx import miller
except (ImportError, RuntimeError) as err:
    logger.info('Error:{}'.format(err))

def mp_corr(x1, x2, q):
    corr = flex.linear_correlation(x1.data(),x2.data())
    if corr.is_well_defined():
        q.put("%5.4f" %corr.coefficient())
    else:
        q.put(0.0)

def CC_calc(file1,file2, **kwargs):
    indata1 = {"xds_ascii": file1}
    indata2 = {"xds_ascii": file2}

    hkl1 = ASCII(indata1)
    hkl2 = ASCII(indata2)
    highres = kwargs.get('highres', '4.0')
    lowres = kwargs.get('lowres', '8.0')
    hkl1.get_data(indata1)
    hkl2.get_data(indata2)
    data1 = hkl1.i_obs(); data2 = hkl2.i_obs()
    #assert (data1.is_similar_symmetry(data2)) #BUG: cause AssertionError for no reason, not always
    try:
        data1 = data1.resolution_filter(d_min=float(highres), d_max=float(lowres))
        data2 = data2.resolution_filter(d_min=float(highres), d_max=float(lowres))
        common_ref_data1, common_ref_data2 = data1.common_sets(data2, assert_is_similar_symmetry=False) #FIXME hack for AssertionError
        corr = flex.linear_correlation(common_ref_data1.data(), common_ref_data2.data())
        if corr.is_well_defined():
            return "%5.4f" %corr.coefficient()
        else:
            return 0.0
    except Exception as e:
        logger.info("no correlation: {}".format(e))
        return 0.0


def mp_cc_calc(file1, file2, **kwargs):
    hkl1 = ASCII(file1); hkl2 = ASCII(file2)
    highres = kwargs.get('highres', '4.0')
    lowres = kwargs.get('lowres', '8.0')
    data1 = hkl1.i_obs(); data2 = hkl2.i_obs()
    assert (data1.is_similar_symmetry(data2))
    data1 = data1.resolution_filter(d_min=float(highres), d_max=float(lowres))
    data2 = data2.resolution_filter(d_min=float(highres), d_max=float(lowres))
    common_ref_data1, common_ref_data2 = data1.common_sets(data2)
    return common_ref_data1, common_ref_data2

class CCestimator(Abstract):
    """docstring for ."""
    @property
    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "xds_ascii": {"type": "string"}
            }
        }

    def __init__(self, inData):
        super().__init__(inData)
        self.results['xscalefile'] = inData['xds_ascii']

        xscale = ASCII(inData)
        xscale.get_data(inData)
        self.results['setlist'] = []
        for k in xscale.results['input_files'].iterkeys():
            self.results['setlist'].append(xscale.results['input_files'][k][0])

        self.results['cc_dataset_list'] = {}
        self.results['cc_cluster_list'] = []
        self.results['cc_dendo'] = []
        self.results['n_clusters_cc'] = 0
        self.results['pcc_arr'] = np.zeros((len(self.results['setlist']), len(self.results['setlist'])))
        self.results['pcc_arr_symm'] = np.zeros((len(self.results['setlist']), len(self.results['setlist'])))
        return

    @staticmethod
    def add_to_queue(func,que,args=None, kwargs=None):
        """
        This method is to put returned values from a function into mp.Queue()
        Queue controls parallel jobs running via mp.Process. This method is set as staticmethod
        so that we can access it without class instance as well as with class instance in other codes
        """
        args=args if args is not None else []
        kwargs=kwargs if kwargs is not None else {}
        que.put(func(*args,**kwargs))

    def ccd_sorter(self, cutoff=0.8):

        for ii in range(len(self.results['setlist'])):
            self.results['cc_dataset_list'][self.results['setlist'][ii]] = CC_calc(self.results['xscalefile'],
                                                                                   self.results['setlist'][ii])

        for k, v in self.results['cc_dataset_list'].items():
            if v < cutoff:
                try:
                    del self.results['cc_dataset_list'][k]
                except (KeyError, ValueError):
                    pass

        return

    def pcc_matrix_mpi(self):
        datasize = len(self.results['setlist'])
        for j in range(datasize):
            start =(j*10)
            end = (j+1)*10
            proc = []
            queue = mp.Queue()
            for ii in range(start, end):
                try:

                    jobid = mp.Process(target=self.add_to_queue, args=[CC_calc, queue, [self.results['setlist'][j],
                                                                                        self.results['setlist'][ii]]],)
                    #jobid.start()
                    proc.append(jobid)
                    #self.pcc_arr[i][j] = CC_calc(self.setlist[i], self.setlist[j])
                    #self.pcc_arr[j][i] = self.pcc_arr[i][j]
                except IndexError:
                    pass
            print(len(proc))
            for p in proc:
                p.start()
            for p in proc:
                p.join()

            for ii in range(start,end):
                if not queue.empty():
                    self.results['pcc_arr'][j][ii] = queue.get()
                    self.results['pcc_arr'][ii][j] = self.results['pcc_arr'][j][ii]

        for k in range(datasize):
            self.results['pcc_arr'][k][k] = 1.0
        return

    def pool_matrix(self):
        datasize = len(self.results['setlist'])
        for j in range(datasize):
            chunk = []
            for i in range(j*10, (j+1)*10):
                try:
                    tmp = [self.results['setlist'][i], self.results['setlist'][j]]
                    chunk.append(tmp)
                except IndexError:
                    pass
            proc = mp.Pool(processes=10)
            output = proc.map(CC_calc, chunk, 1)
            self.results['pcc_arr'][j, j*10:(j+1)*10] = output
        return

    def pcc_matrix(self):
        datasize = len(self.results['setlist'])

        for j in range(datasize):
            for i in range(j+1):

                if i == j:
                    self.results['pcc_arr'][j,i] = 1.0
                else:
                    self.results['pcc_arr'][j,i] = CC_calc(self.results['setlist'][j], self.results['setlist'][i])

        #fast trick to generate symmetric array
        self.results['pcc_arr_symm'] = self.results['pcc_arr'] + self.results['pcc_arr'].T - \
                                       np.diag(self.results['pcc_arr'].diagonal())
        return

    def cc_cluster(self):
        self.results['data_points'] = []
        for item in self.results['setlist']:
            tmp = os.path.basename(item)
            self.results['data_points'].append(tmp.strip('.HKL'))

        try:
            Y = sch.linkage(self.results['pcc_arr_symm'], metric='correlation', method='average')
        except ValueError as e:
            logger.info('Error-CC linkage:{}'.format(e))
            logger.info('setting NaN or Inf values to zero and rerunning CC clustering')
            clean_array = np.nan_to_num(self.results['pcc_arr_symm'])
            Y = sch.linkage(clean_array, metric='correlation', method='average')
        try:
            assign = sch.fcluster(Y, 0.8, 'distance')
            MSG = "# clusters: %d" %assign.max()
            self.results['n_clusters_cc'] = int(assign.max())
            logger.info('MSG: {}'.format(MSG))

            idx = Cell.cluster_indices(assign)
            max_size = []; self.results['cc_cluster_list'] = []
            for k, i in enumerate(idx):
                max_size.append(len(i))
                logger.info("cluster: %d; size: %d" %((k+1),len(i)))
            msg = "most populated cluster: %d" %(max_size.index(max(max_size))+1)
            logger.info('MSG: {}'.format(msg))

            best_cluster_name = max_size.index(max(max_size))
            for i in range(len(idx[best_cluster_name])):

                item = idx[best_cluster_name][i]
                cc_cluster_hkl = self.results['setlist'][item]
                self.results['cc_cluster_list'].append(cc_cluster_hkl)

            logger.info('data_points_cc:{}'.format(self.results['data_points']))
            self.results['cc_dendo'] = sch.dendrogram(Y, labels=self.results['data_points'], no_plot=True)
            #sch.dendrogram(Y, truncate_mode='level', show_contracted=True, leaf_rotation=90)
        except Exception as e:
            logger.info('From_cc_cluster_func:{}')

    def cc_select(self, fom='ccd'):
        if fom == 'ccd':
            self.ccd_sorter()
            fh = open('cc-dataset.INP', 'w')
            for fname in self.results['cc_dataset_list'].keys():
                fh.write('INPUT_FILE=%s\n' %fname)
                fh.write('SNRC=0.0\n')
            fh.close()

        elif fom == 'pcc':
            logger.info('pcc_matrix called')
            self.pcc_matrix()

            logger.info('pcc_array:{}'.format(self.results['pcc_arr']))

            logger.info('pairCC class called')
            self.cc_cluster()
            fh = open('cc-sorted.INP', 'w')

            for fname in self.results['cc_cluster_list']:
                fh.write('INPUT_FILE=%s\n' %fname)
                fh.write('SNRC=0.0\n')
            fh.close()
        else:
            pass
        return

if '__name__'=='__main__':
    logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%a, %d %b %Y %H:%M:%S',
    filename='correlate.log',
    filemode='w')
    inData = {"xds_ascii": sys.argv[1]}
    CC = CCestimator(inData)
    #print CC.ccd_sorter()
    CC.cc_select(fom='pcc')
