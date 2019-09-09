'''
Created by S.Basu
On 28-June-2017
'''

import sys, os
import numpy as np
import subprocess as sub
import matplotlib.pyplot as plt
from scipy import stats
import scipy.cluster.hierarchy as sch
import logging

logger = logging.getLogger('Scale&Merge')

class Cell(object):

	def __init__(self, files):
		self.status = False
		if type(files) is list:
			self.hkllist = files
			self.status = True
		elif type(files) is str and files.endswith(".LP"):
			self.hkllist = self.get_filenames(files)

		else:
			err = TypeError("Need correct type of input, either list of HKL files or XSCALE.LP file")
			logger.info('TypeError: {}'.format(err))
			self.status = False
			return self.status


		if len(self.hkllist) > 0:
			self.cell_ar = np.empty((len(self.hkllist),6))
			self.cell_vector = np.empty((len(self.hkllist),3))
			self.cell_select = []; self.dendo = {};
			self.hclust_matrix = np.empty((len(self.hkllist)-1,4))
			self.status = True
		else:
			err = ValueError("Reprocess with a lower ISa_cutoff, no XDS_ASCII.HKL file has enough ISa\n")
			logger.info('ValueError: {}'.format(err))
			self.status = False

	def get_filenames(self, filename):
		fLst1 = []; fLst2 = [];
		if os.path.isfile(filename):
			fh = open(filename, 'r')
			all_lines = fh.readlines()
			fh.close();

			for lines in all_lines:
				if "INPUT_FILE" in lines:
					line = lines.split('=')
					fLst1.append(line[1])
				#else:
				#	pass
			for name in fLst1:
				tmp = name.strip('\n')
				fLst2.append(tmp)
		else:
			err = OSError("%s does not exist" %filename)
			logger.info('IOError: {}'.format(err))
			self.status = False
			return fLst2
		return fLst2

	def get_cells(self, filename):
		try:
			fh = open(filename, 'r')
			all_lines = fh.readlines()
			fh.close()
			for lines in all_lines:
				if "!UNIT_CELL_CONSTANTS" in lines:
					line = lines.split()
					self.a = float(line[1]); self.b = float(line[2]); self.c = float(line[3]);
					self.al = float(line[4]); self.be = float(line[5]); self.ga = float(line[6]);
				elif "!SPACE_GROUP_NUMBER" in lines:
					line = lines.split()
					self.spg = int(line[1])
		except Exception as e:
			logger.info('Error: {}'.format(e))
			self.status = False
			return self.status
	def calc_vector(self):
		Vab = np.sqrt(self.a**2 + self.b**2 - 2*self.a*self.b*np.cos(np.pi-np.deg2rad(self.ga)))
		Vbc = np.sqrt(self.b**2 + self.c**2 - 2*self.b*self.c*np.cos(np.pi-np.deg2rad(self.al)))
		Vca = np.sqrt(self.c**2 + self.a**2 - 2*self.c*self.a*np.cos(np.pi-np.deg2rad(self.be)))

		return Vab, Vbc, Vca

	def reject_outlier(self, array):
		ohne_outlier = array[abs(array - np.mean(array)) <= 1.5*np.std(array)]
		ohne_indx = np.where(abs(array - np.mean(array)) <= 1.5*np.std(array))
		ohne_indx = ohne_indx[0].tolist()

		return ohne_outlier, ohne_indx

	def lcv_(self):

		self.cell_vector = np.empty((len(self.hkllist),3))

		for i in range(len(self.hkllist)):
			self.get_cells(self.hkllist[i])
			Vab, Vbc, Vca = self.calc_vector()
			if self.status == False:
				return self.status

			self.cell_vector[i,0] = Vab; self.cell_vector[i,1] = Vbc;
			self.cell_vector[i,2] = Vca;
		return self.status

	def cluster_indices(self, cluster_assignments):
		n = cluster_assignments.max();
		indices = []
		for cluster_number in range(1, n + 1):
			indices.append(np.where(cluster_assignments == cluster_number)[0])

		return indices

	def clustering(self):
		self.lcv_()
		self.data_points = [];
		for item in self.hkllist:
			tmp = os.path.basename(item)
			self.data_points.append(tmp.strip('.HKL'))

		if self.status == False:
			return self.status
		#hierarchical clustering using resultant vectors from unit cells
		Y = sch.linkage(self.cell_vector, method='ward')
		self.hclust_matrix = Y
		assign = sch.fcluster(Y, 2, 'distance')
		MSG = "# clusters: %d" %assign.max()
		self.n_clusters_cell = int(assign.max())
		logger.info('MSG: {}'.format(MSG))
		idx = self.cluster_indices(assign)
		max_size = [];
		for k, i in enumerate(idx):
			max_size.append(len(i))
			logger.info("cluster: %d; size: %d" %((k+1),len(i)))
		msg = "most populated cluster: %d" %(max_size.index(max(max_size))+1)
		logger.info('MSG: {}'.format(msg))
		best_cluster_name = max_size.index(max(max_size))
		self.cell_ar_best_cluster = np.empty((len(idx[best_cluster_name]), 6))
		for i in range(len(idx[best_cluster_name])):
			item = idx[best_cluster_name][i]
			cell_selected_file = self.hkllist[item]
			self.cell_select.append(cell_selected_file)
			self.get_cells(cell_selected_file)
			self.cell_ar_best_cluster[i,0] = self.a; self.cell_ar_best_cluster[i,1] = self.b;
			self.cell_ar_best_cluster[i,2] = self.c; self.cell_ar_best_cluster[i,3] = self.al;
			self.cell_ar_best_cluster[i,4] = self.be; self.cell_ar_best_cluster[i,5] = self.ga;
		#calculate median values for unit cell parameters from most populated cluster..
		self.a_mode = stats.mode(self.cell_ar_best_cluster[:,0])[0][0]
		self.b_mode = stats.mode(self.cell_ar_best_cluster[:,1])[0][0]
		self.c_mode = stats.mode(self.cell_ar_best_cluster[:,2])[0][0]
		self.al_mode = stats.mode(self.cell_ar_best_cluster[:,3])[0][0]
		self.be_mode = stats.mode(self.cell_ar_best_cluster[:,4])[0][0]
		self.ga_mode = stats.mode(self.cell_ar_best_cluster[:,5])[0][0]

		logger.info('Cell selection # HKLs: %d' %len(self.cell_select))

		self.dendo = sch.dendrogram(Y, p=10, labels=self.data_points, no_plot=True)
		#logger.info('dendrogram switched off')
		self.status = True
		return self.status

	def cell_analysis(self):


		tmp_fileList = []; #Local variable - a list of lists

		for i in range(len(self.hkllist)):
			self.get_cells(self.hkllist[i])
			if self.status == False:
				return self.status

			self.cell_ar[i,0] = self.a; self.cell_ar[i,1] = self.b; self.cell_ar[i,2] = self.c;
			self.cell_ar[i,3] = self.al; self.cell_ar[i,4] = self.be; self.cell_ar[i,5] = self.ga;

		self.a_ohne, a_indx = self.reject_outlier(self.cell_ar[:,0]);
		self.b_ohne, b_indx = self.reject_outlier(self.cell_ar[:,1]);
		self.c_ohne, c_indx = self.reject_outlier(self.cell_ar[:,2]);
		self.al_ohne, al_indx = self.reject_outlier(self.cell_ar[:,3]);
		self.be_ohne, be_indx = self.reject_outlier(self.cell_ar[:,4]);
		self.ga_ohne, ga_indx = self.reject_outlier(self.cell_ar[:,5])
		#calculate median values for unit cell parameters..
		self.a_mode = stats.mode(self.a_ohne)[0][0]
		self.b_mode = stats.mode(self.b_ohne)[0][0]
		self.c_mode = stats.mode(self.c_ohne)[0][0]
		self.al_mode = stats.mode(self.al_ohne)[0][0]
		self.be_mode = stats.mode(self.be_ohne)[0][0]
		self.ga_mode = stats.mode(self.ga_ohne)[0][0]

		for index in (a_indx, b_indx, c_indx, al_indx, be_indx, ga_indx):
			if index: # no outlier may lead to an empty array, avoiding it..
				tmp_fileList.append([self.hkllist[i] for i in index])

		self.cell_select = set(tmp_fileList[0]) #use set method to find common HKL files after outlier rejection
		for s in tmp_fileList[1:2]:
			self.cell_select.intersection_update(s)
			self.cell_select = list(self.cell_select)

		return self.status

	def dict_for_histogram(self):
		self.hist_dict = {}
		if self.cell_ar.size == 0:
			msg = "Cell selection did not work, so cell_array is empty"
			logger.info('MSG:{}'.format(msg))
			self.status = False
			return self.hist_dict
		self.hist_dict['a'] = list(self.cell_ar[:,0])
		self.hist_dict['b'] = list(self.cell_ar[:,1])
		self.hist_dict['c'] = list(self.cell_ar[:,2])
		self.hist_dict['al'] = list(self.cell_ar[:,3])
		self.hist_dict['be'] = list(self.cell_ar[:,4])
		self.hist_dict['ga'] = list(self.cell_ar[:,5])
		self.status = True
		return self.hist_dict


	def cell_histogram(self):
		if len(self.hkllist) == 0:
			err = ValueError("xscaling with ISa-selection did not work, reduce isa-cutoff")
			logger.info('ValueError: {}'.format(err))
			self.status = False
			return self.status

		fig, axs = plt.subplots(2,3, squeeze=True, facecolor='w', edgecolor='k')
		fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.4, wspace=0.3)
		axs = axs.flatten()
		self.cell_analysis()
		if self.status == False:
			return self.status


		axs[0].hist(self.cell_ar[:,0], 50, histtype='step')
		axs[0].set_xlabel('a-axis(Ang)', fontweight='bold')
		axs[0].set_ylabel('frequency', fontweight='bold')
		axs[0].set_title("a=%5.2f A" %self.a_mode)
		try:
			axs[0].hist(self.a_ohne)
		except ValueError, IndexError:
			pass

		axs[1].hist(self.cell_ar[:,1], 50, histtype='step')
		axs[1].set_xlabel('b-axis(Ang)', fontweight='bold')
		axs[1].set_ylabel('frequency', fontweight='bold')
		axs[1].set_title("b=%5.2f A" %self.b_mode)
		try:
			axs[1].hist(self.b_ohne)
		except ValueError, IndexError:
			pass

		axs[2].hist(self.cell_ar[:,2], 50, histtype='step')
		axs[2].set_xlabel('c-axis(Ang)', fontweight='bold')
		axs[2].set_ylabel('frequency', fontweight='bold')
		axs[2].set_title("c=%5.2f A" %self.c_mode)
		try:
			axs[2].hist(self.c_ohne)
		except ValueError, IndexError:
			pass

		axs[3].hist(self.cell_ar[:,3],50, histtype='step')
		axs[3].set_xlabel('alpha(deg)', fontweight='bold')
		axs[3].set_ylabel('frequency', fontweight='bold')
		axs[3].set_title("alpha=%5.2f A" %self.al_mode)
		try:
			axs[3].hist(self.al_ohne)
		except ValueError, IndexError:
			pass
		axs[4].hist(self.cell_ar[:,4], 50, histtype='step')
		axs[4].set_xlabel('beta(deg)', fontweight='bold')
		axs[4].set_ylabel('frequency', fontweight='bold')
		axs[4].set_title("beta=%5.2f A" %self.be_mode)
		try:
			axs[4].hist(self.be_ohne)
		except ValueError, IndexError:
			pass

		axs[5].hist(self.cell_ar[:,5], 50, histtype='step')
		axs[5].set_xlabel('gamma[deg]', fontweight='bold')
		axs[5].set_ylabel('frequency', fontweight='bold')
		axs[5].set_title("gamma=%5.2f A" %self.ga_mode)
		try:
			axs[5].hist(self.ga_ohne)
		except ValueError, IndexError:
			pass

		plt.savefig("cell-histogram.png")
		'''
		# run_command will not work because this is a python command and not bash shell command.
		try:
			run_command("Scale&Merge", os.getcwd(), os.environ['USER'], cmd, 'merge.log')
		except Exception:
			sub.call(cmd, shell=True)
		return self.status
		'''


def finder():
	root = sys.argv[1]
	import glob
	paths = glob.glob(os.path.join(root, 'minisets_*/xtal_*/XDS_ASCII.HKL'))
	return paths

import glob
paths = glob.glob('/Users/shibom/work/CY_IMISX/all_pepT_sad_v2/xtal_*.HKL')
cell = Cell(sorted(paths[0:100]))
cell.clustering()
