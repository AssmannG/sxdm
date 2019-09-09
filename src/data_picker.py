'''
Created by S.Basu
On 1-Aug-2017
'''
import os, sys
import logging
import numpy as np
from ascii import ASCII
from cellprobe import Cell

logger = logging.getLogger('Scale&Merge')

class pairCC(object):

	def __init__(self, fname):
		self.filename = fname
		self.cc_dict = {};
		self.error_b_dict = {};
		self.datasets = [];
		self.hkl_b_sorted = [];
		self.hkl_cc_sorted = [];

	def get_cc_error_b(self):

		if self.filename.endswith('.LP'):
			fh = open(self.filename, 'r')
			_all = fh.readlines()
			fh.close()

			try:
				cc_start = _all.index("  #i   #j     REFLECTIONS     BETWEEN i,j  INTENSITIES (i/j)  BETWEEN i,j\n")
				cc_end = _all.index(" K*EXP(B*SS) = Factor applied to intensities\n");
				cc_chunk = _all[cc_start+1:cc_end-5]
				b_start = _all.index('     a        b          ISa    ISa0   INPUT DATA SET\n')
				b_end = _all.index(' (ASSUMING A PROTEIN WITH 50% SOLVENT)\n')
				error_b_chunk = _all[b_start+1:b_end-3]
			except (ValueError, IndexError) as err:
				error = "check if XSCALE ran properly or XSCALE.INP screwed up \n"
				logger.info('Error: {}'.format(error))
				return
			try:
				for lines in cc_chunk:
					line = lines.split()
					try:
						data_key = tuple((int(line[0]),int(line[1])))
						self.cc_dict[data_key] = float(line[3])

					except IndexError:
						pass

				for lines in error_b_chunk:
					line = lines.split()
					try:
						hkl_key = line[4]
						self.datasets.append(hkl_key)
						self.error_b_dict[hkl_key] = float(line[1])
					except IndexError:
						pass

			except Exception as err:
				logger.info('Error:{}'.format(err))

		else:
			err = TypeError('wrong file type, XSCALE.LP needed')
			logging.info('ERROR:{}'.format(err))

		return

	def outliers_iqr(self, ys): #Tukey's outlier rejection
		quartile_1, quartile_3 = np.percentile(ys, [25, 75])
		iqr = quartile_3 - quartile_1
		lower_bound = quartile_1 - (iqr * 1.5)
		upper_bound = quartile_3 + (iqr * 1.5)
		return np.where((ys < upper_bound) | (ys > lower_bound))

	def error_b_sorter(self):

		b_vals = self.error_b_dict.values()
		b_array = np.asarray(b_vals)
		good_b_indx = self.outliers_iqr(b_array)  #index of good ones but as tuple
		clean_b = good_b_indx[0].tolist()

		for b in clean_b:
			hkl = self.error_b_dict.keys()[clean_b.index(b)]
			self.hkl_b_sorted.append(hkl)
		return

	def pair_corr_sorter(self, cutoff=0.9):
		mean_ccs = []; self.hkl_cc_sorted = self.datasets[0:2];
		for j in range(2,len(self.datasets)+1):
			cc_j = []
			for i in range(1,j):
				try:
					cc_j.append(self.cc_dict[(i,j)])
				except KeyError:
					pass
			mean_ccs.append(sum(cc_j)/len(cc_j))

		for ii in range(len(mean_ccs)):
			if mean_ccs[ii] > cutoff:
				try:
					self.hkl_cc_sorted.append(self.datasets[ii+2])
				except (IndexError, ValueError) as err:
					pass
		return

	def pair_cc_matrix(self):
		datasize = len(self.datasets)
		self.cc_arr = np.empty((datasize, datasize))
		for j in range(2, datasize):
			for i in range(1,j):
				try:
					self.cc_arr[i][j] = self.cc_dict[(i,j)]
					self.cc_arr[j][i] = self.cc_dict[i][j]
				except KeyError:
					pass
		for k in range(datasize):
			self.cc_arr[k][k] = 1.0
		return

	def cc_cluster(self, array, datasets):
		import scipy.cluster.hierarchy as sch
		import matplotlib.pyplot as plt
		self.data_points = [];
		for item in datasets:
			tmp = os.path.basename(item)
			self.data_points.append(tmp.strip('.HKL'))

		try:
			Y = sch.linkage(array, metric='correlation', method='average')
		except ValueError as e:
			logger.info('Error-CC linkage:{}'.format(e))
			logger.info('setting NaN or Inf values to zero and rerunning CC clustering')
			clean_array = np.nan_to_num(array)
			Y = sch.linkage(clean_array, metric='correlation', method='average')
		try:
			assign = sch.fcluster(Y, 0.8, 'distance')
			MSG = "# clusters: %d" %assign.max()
			self.n_clusters_cc = int(assign.max())
			logger.info('MSG: {}'.format(MSG))
			cell = Cell(datasets)
			idx = cell.cluster_indices(assign)
			max_size = []; self.cc_cluster_list = [];
			for k, i in enumerate(idx):
				max_size.append(len(i))
				logger.info("cluster: %d; size: %d" %((k+1),len(i)))
			msg = "most populated cluster: %d" %(max_size.index(max(max_size))+1)
			logger.info('MSG: {}'.format(msg))

			best_cluster_name = max_size.index(max(max_size))
			for i in range(len(idx[best_cluster_name])):

				item = idx[best_cluster_name][i]
				cc_cluster_hkl = datasets[item]
				self.cc_cluster_list.append(cc_cluster_hkl)

			logger.info('data_points_cc:{}'.format(self.data_points))
			self.cc_dendo = sch.dendrogram(Y, labels=self.data_points, no_plot=True)
			#sch.dendrogram(Y, truncate_mode='level', show_contracted=True, leaf_rotation=90)
		except Exception as e:
			logger.info('From_cc_cluster_func:{}'.format(e))
		return

	def make_xscale(self, List1):
		if len(List1) > 0:
			fh = open('XSCALE.INP', 'w')
			fh.write('OUTPUT_FILE=XSCALE.HKL\n')
			fh.write("FRIEDEL'S_LAW=FALSE\n")
			fh.write('SAVE_CORRECTION_IMAGES=FALSE\n')
			for name in List1:
				fh.write('INPUT_FILE=%s\n' %name)
			fh.close()
		else:
			logger.info('IndexError: cc-sorted hkl list was empty\n')

	def calc_CCset(self, hkl1, hkl2):
		f1 = ASCII(hkl1)
		f2 = ASCII(hkl2); cclst = []

		for index in set(f1.data_dict).intersection(set(f2.data_dict)):
			I1 = f1.data_dict[index][0]; sigI1 = f1.data_dict[index][1];
			I2 = f2.data_dict[index][0]; sigI2 = f2.data_dict[index][1];
			cc = (I1*I2)/(sigI1*sigI2)
			cclst.append(cc)
		return cclst

def main():
	pcc = pairCC(sys.argv[1])
	pcc.get_cc_error_b()
	'''
	hkl_b_sorted = error_b_sorter(b_dict)
	hkl_cc_sorted = pair_corr_sorter(cc_dict, datasets)
	fh = open('b-sorted.INP', 'w')
	for fname in hkl_b_sorted:
		fh.write('INPUT_FILE=%s\n' %fname)
	fh.close()

	fh = open('cc-sorted.INP', 'w')

	for fname in hkl_cc_sorted:
		fh.write('INPUT_FILE=%s\n' %fname)
	fh.close()
	print len(hkl_b_sorted), len(hkl_cc_sorted)
	'''
	pcc.pair_cc_matrix()
	D_arr = np.empty(pcc.cc_arr.shape)
	for j in range(pcc.cc_arr.shape[0]):
		for i in range(pcc.cc_arr.shape[1]):
			D_arr[i][j] = np.sqrt(1-(pcc.cc_arr[i][j]**2))

	pcc.cc_cluster(pcc.cc_arr, pcc.datasets)
	fh = open('cc-sorted.INP', 'w')

	for fname in pcc.cc_cluster_list:
		fh.write('INPUT_FILE=%s\n' %fname)
		fh.write('MINIMUM_I/SIGMA=0.0\n')
	fh.close()

	return


if __name__=='__main__':
	main()
