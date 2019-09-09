'''
created by S.Basu
Date: 20 Feb, 2018
'''
import sys, os, time
import errno, numpy as np
import re, logging
from ascii import ASCII
import data_picker as dp
import multiprocessing as mp

logger = logging.getLogger('Scale&Merge')

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
	#hkl1 = ASCII(lst[0]); hkl2 = ASCII(lst[1])
	hkl1 = ASCII(file1); hkl2 = ASCII(file2)
	highres = kwargs.get('highres', '4.0')
	lowres = kwargs.get('lowres', '8.0')
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
	except Exception:
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

class CC_estimator(object):
	"""docstring for ."""
	def __init__(self, xscalefile):
		self.xscalefile = xscalefile
		self.xscale = ASCII(xscalefile)
		self.setlist = [];
		for k in self.xscale.input_files.iterkeys():
			self.setlist.append(self.xscale.input_files[k][0])

		self.cc_dataset_list = {}

	@staticmethod
	def add_to_queue(func,que,args=None, kwargs=None):
		'''
		This method is to put returned values from a function into mp.Queue()
		Queue controls parallel jobs running via mp.Process. This method is set as staticmethod
		so that we can access it without class instance as well as with class instance in other codes
		'''
		args=args if args is not None else []
		kwargs=kwargs if kwargs is not None else {}
		que.put(func(*args,**kwargs))

	def ccd_sorter(self, cutoff=0.8):
		for ii in range(len(self.setlist)):
			self.cc_dataset_list[self.setlist[ii]] = CC_calc(self.xscalefile, self.setlist[ii])

		for k, v in self.cc_dataset_list.items():
			if v < cutoff:
				try:
					del self.cc_dataset_list[k]
				except (KeyError, ValueError) as err:
					pass

		return self.cc_dataset_list.keys()

	def pcc_matrix_mpi(self):
		datasize = len(self.setlist)
		self.pcc_arr = np.zeros((datasize, datasize))
		#self.pcc_arr[0][1] = CC_calc(self.setlist[0], self.setlist[1])
		#self.pcc_arr[1][0] = self.pcc_arr[0][1]
		for j in range(datasize):
			start =(j*10);
			end = (j+1)*10;
			proc = [];
			queue = mp.Queue()
			for ii in range(start, end):
				try:

					jobid = mp.Process(target=self.add_to_queue, args=[CC_calc, queue, [self.setlist[j],self.setlist[ii]]],)
					#jobid.start()
					proc.append(jobid)
					#self.pcc_arr[i][j] = CC_calc(self.setlist[i], self.setlist[j])
					#self.pcc_arr[j][i] = self.pcc_arr[i][j]
				except IndexError:
					pass
			print len(proc)
			for p in proc:
				p.start()
			for p in proc:
				p.join()

			for ii in range(start,end):
				if not queue.empty():
					self.pcc_arr[j][ii] = queue.get()
					self.pcc_arr[ii][j] = self.pcc_arr[j][ii]

		for k in range(datasize):
			self.pcc_arr[k][k] = 1.0
		return

	def pool_matrix(self):
		datasize = len(self.setlist)
		self.pcc_arr = np.zeros((datasize, datasize))
		for j in range(datasize):
			chunk = [];
			for i in range(j*10, (j+1)*10):
				try:
					tmp = [self.setlist[i], self.setlist[j]]
					chunk.append(tmp)
				except IndexError:
					pass
			proc = mp.Pool(processes=10)
			output = proc.map(CC_calc, chunk, 1)
			self.pcc_arr[j, j*10:(j+1)*10] = output
		return

	def pcc_matrix(self):
		datasize = len(self.setlist)
		self.pcc_arr = np.zeros((datasize, datasize))
		for j in range(datasize):
			for i in range(j+1):

				if i == j:
					self.pcc_arr[j,i] = 1.0
				else:
					self.pcc_arr[j,i] = CC_calc(self.setlist[j], self.setlist[i])

		#fast trick to generate symmetric array
		self.pcc_arr_symm = self.pcc_arr + self.pcc_arr.T - np.diag(self.pcc_arr.diagonal())
		return

	def cc_select(self, LPfile, fom='ccd'):
		if fom == 'ccd':
			self.ccd_sorter()
			fh = open('cc-dataset.INP', 'w')
			for fname in self.cc_dataset_list.keys():
				fh.write('INPUT_FILE=%s\n' %fname)
				fh.write('MINIMUM_I/SIGMA=0.0\n')
			fh.close()

		elif fom == 'pcc':
			logger.info('pcc_matrix called')
			self.pcc_matrix()
			#pcc = dp.pairCC(LPfile)
			logger.info('pcc_array:{}'.format(self.pcc_arr))
			pcc = dp.pairCC(LPfile)
			logger.info('pairCC class called')
			pcc.cc_cluster(self.pcc_arr_symm, self.setlist)
			fh = open('cc-sorted.INP', 'w')
			self.cc_cluster_list = pcc.cc_cluster_list
			self.cc_dendo = pcc.cc_dendo
			self.n_clusters_cc = pcc.n_clusters_cc
			for fname in pcc.cc_cluster_list:
				fh.write('INPUT_FILE=%s\n' %fname)
				fh.write('MINIMUM_I/SIGMA=0.0\n')
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
	CC = CC_estimator(sys.argv[1])
	#print CC.ccd_sorter()
	CC.cc_select(sys.argv[2], fom='pcc')
