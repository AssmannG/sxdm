'''
Date Sept-2017
S.Basu
'''

import sys, os, time
import errno, numpy as np
import re, logging
from run_command import *
import subprocess as sub



logger = logging.getLogger("Scale&Merge")
try:
	from cctbx.array_family import flex
	from cctbx import crystal
	from cctbx import miller
except (ImportError, RuntimeError) as err:
	logger.info('Error:{}'.format(err))

class ASCII(object):
	def __init__(self, xds_ascii):
		self.fname = xds_ascii
		self.header = {}; self.input_files = {};
		self.unit_cell = {};
		self.anom = "TRUE"
		self.data_dict = {} #{(hkl): I, sigI}
		self.indices = []; self.iobs = [];
		self.sigI = [];
		if os.path.isfile(self.fname):
			self.readfile()
			self._extract()
		else:
			err = "XDS_ASCII file cannot be read\n"
			logger.info('Error:{}'.format(err))
		return


	def readfile(self):
		fh = open(self.fname, 'r')
		_all = fh.readlines()
		fh.close()
		start = _all.index('!END_OF_HEADER\n')
		try:
			end = _all.index('!END_OF_DATA\n')
		except ValueError:
			end = _all.index('!END_OF_DATA \n')

		chunk = _all[start+1:end]
		kwd = re.compile("([^ =]+)= *((?:(?! [^ =]+=).)*)")

		for lines in _all:
			if lines.startswith("! ISET"):
				params = dict(kwd.findall(lines))
				iset = int(params['ISET'])
				if iset not in self.input_files:
					self.input_files[iset] = [None, None]
				if "INPUT_FILE" in params:
					self.input_files[iset][0] = params["INPUT_FILE"]
				elif "X-RAY_WAVELENGTH" in params:
					self.input_files[iset][1] = params["X-RAY_WAVELENGTH"]
				else:
					pass
			elif "!" in lines:
				lst = kwd.findall(lines)
				try:
					for item in lst:
						self.header[item[0]] = item[1]
				except IndexError:
					pass
			else:
				pass

		for field in chunk:
			line = field.split()
			index = [int(line[0]), int(line[1]), int(line[2])]
			self.indices.append(index)
			self.iobs.append(float(line[3]))
			self.sigI.append(float(line[4]))
			#self.data_dict[index] = (float(line[3]), float(line[4]))
		self.indices = flex.miller_index(self.indices)
		self.iobs, self.sigI = map(lambda x:flex.double(x), (self.iobs, self.sigI))
		return

	def _extract(self):
		try:
			self.anom = self.header["FRIEDEL'S_LAW"]
			self.spg = self.header["!SPACE_GROUP_NUMBER"]
			self.cell = self.header["!UNIT_CELL_CONSTANTS"]
			self.cell = self.cell.split()
			self.unit_cell['a'] = self.cell[0]
			self.unit_cell['b'] = self.cell[1]
			self.unit_cell['c'] = self.cell[2]
			self.unit_cell['al'] = self.cell[3]
			self.unit_cell['be'] = self.cell[4]
			self.unit_cell['ga'] = self.cell[5]
			#cctbx crystal symmetry into flex array type
			a,b,c,al,be,ga = map(lambda x:float(x), self.cell)
			self.symm = crystal.symmetry(unit_cell=(a,b,c,al,be,ga), space_group=int(self.spg))

			self.wave = self.header["!X-RAY_WAVELENGTH"]
		except (KeyError, IndexError) as e:
			err = "ASCII-class: could be XSCALE file or check file for keyErrors\n"
			logger.info('INFO:{}'.format(err))

	def create_miller_set(self):
		if self.anom == "TRUE":
			anom_flag = True
		elif self.anom == "FALSE":
			anom_flag = False
		else:
			anom_flag = None

		return miller.set(crystal_symmetry=self.symm, indices=self.indices, anomalous_flag=anom_flag)

	def i_obs(self):
		array_info = miller.array_info(source_type="xds_ascii")
		return miller.array(self.create_miller_set(),
							data=self.iobs, sigmas=self.sigI).set_info(array_info).set_observation_type_xray_intensity()

	def multiplicity_check(self):
		data = self.i_obs()
		merge = data.merge_equivalents(use_internal_variance=False)
		self.mult_list = []
		redundancy = merge.redundancies().data()
		try:
			for x in sorted(set(redundancy)):
				each_data_point = [x, redundancy.count(x)]
				self.mult_list.append(each_data_point)
			#import matplotlib.pyplot as plt
			#plt.hist(self.mult_dict.values())
			#plt.savefig('multiplicity-distribution.png')
		except Exception as err:
			logger.info('multiplicity_check_error:{}'.format(err))
		return

	def xdsconv(self, form, res=None, dirname=None, user=None):

		if os.path.isfile(self.fname):
			try:
				linkname = "data.hkl"
				os.symlink(self.fname, linkname)
			except OSError, e:
				if e.errno == errno.EEXIST:
					os.remove(linkname)
					os.symlink(self.fname, linkname)
				else:
					logger.info('xdsconv-error:{}'.format(e))
		else:
			logger.info('xdsconv_error:{}'.format('file does not exist'))
			return
		namelist = os.path.basename(self.fname).split('.')
		self.rootname = namelist[0];
		if res is None:
			res = 0.0
		else:
			pass

		xdsconv_cmd = 'xdsconv > /dev/null'
		f2mtz_cmd = 'f2mtz HKLOUT %s < F2MTZ.INP > /dev/null' %(self.rootname+'_'+form+'.mtz')

		fh = open('XDSCONV.INP', 'w')
		fh.write("OUTPUT_FILE=temp.hkl  %s\n" %(form))
		fh.write("INPUT_FILE=%s\n" %(linkname))
		fh.write("FRIEDEL'S_LAW=%s\n" %(self.anom))
		fh.write("INCLUDE_RESOLUTION_RANGE=50 %s\n" %(res))
		fh.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS= 0.05\n")
		fh.close()
		logger.info("running xdsconv")

		try:
			run_command("Scale&Merge", dirname, user, xdsconv_cmd, 'xdsconv.log')
		except Exception:
			sub.call(xdsconv_cmd, shell=True)
		time.sleep(2)
		if os.path.isfile('F2MTZ.INP'):
			try:
				run_command("Scale&Merge", dirname, user, f2mtz_cmd, 'xdsconv.log')
			except Exception:
				sub.call(f2mtz_cmd, shell=True)
			time.sleep(2)
			os.remove('temp.hkl')
		else:
			logger.info('f2mtz_error:{}'.format('xdsconv did not run\n'))
		return
