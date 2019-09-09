'''
created by S.Basu
On 26-June-2017
'''
import sys, os, errno
import glob, shutil
import subprocess as sub
import logging
from cellprobe import Cell
import index_check
from xscale_output import *
from run_command import *
from scale_utl import *
from data_picker import *


logger = logging.getLogger('Scale&Merge')


class Merge_utls(object):

	def __init__(self, pathlist, experiment, **kwargs):
		self.status = False
		self.input_paths = pathlist
		self.hklpaths = []
		self.filelinks = [];
		self.selected_files = []
		self.scale_results = {'status': False, 'info': []}
		self.expt = experiment
		self.user = kwargs.get('user', os.environ['USER'])
		self.xscale_cmd = "xscale_par > /dev/null"
		self.support_expt = ['native-sad', 'inverse-beam', 'interleave-and-inverse-first', 'serial-xtal']
		
		self.isa_cut = kwargs.get('isa_cut', '4.0')
		self.reference = kwargs.get('reference', None)
		self.running_folder = kwargs.get('running_folder', None)

		self.scale_results['method'] = self.expt

	def outdir_make(self):
		outname = 'adm-'+self.expt
		if self.running_folder != None:
			self.running_folder = os.path.abspath(self.running_folder)
			if os.path.isdir(self.running_folder):
				self.output = os.path.join(self.running_folder, outname)

			else:
				err = "Error: wrong running folder name given. check it\n"
				logger.info({}.format(err))
				self.status = False

			if not os.path.exists(self.output):
				os.makedirs(self.output, 0755)
			else:
				pass
		else:
			cwd = os.getcwd()
			self.output = os.path.join(cwd, outname)
		if not os.path.exists(self.output):
			os.makedirs(self.output, 0755)
		else:
			pass
		try:
			os.chdir(self.output)
			self.subadm = 'adm-'+str(len(self.hklpaths))
			self.subadm = os.path.join(self.output, self.subadm)
			if not os.path.exists(self.subadm):
				os.makedirs(self.subadm, 0755)
			else:
				pass

		except OSError:
			err = OSError("adm folder either not there or not accessible\n")
			logger.info('OSError:{}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results

		self.scale_results['merge_output_path'] = self.subadm

	def sort_xtal(self, lists):
		sortlist = [];
		if lists:
			for val in lists:
				num_sep = val.split('_')
				num = num_sep[1].split('.')
				sortlist.append(int(num[0]))
				tmp = sorted(sortlist)
			sortlist = [];
			for val in tmp:
				element = 'xtal_'+str(val)+'.HKL'
				sortlist.append(element)
		return sortlist

	def find_HKLs(self):
		if not self.input_paths:
			err = IndexError("no hkl filename provided\n")
			logger.info('IndexError: {}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results
		else:
			msg = "# of xtals found: %d\n" %len(self.input_paths)
			logger.info('MSG: {}'.format(msg))
			self.status = True; self.scale_results['status'] = self.status
			self.scale_results['xtals_found'] = len(self.input_paths)
		for path in self.input_paths:
			if os.path.isdir(path):
				filepath = os.path.join(path, "XDS_ASCII.HKL")
				if os.path.isfile(filepath):
					self.hklpaths.append(os.path.join(path, "XDS_ASCII.HKL"))
					self.scale_results['status'] = True
				else:
					self.scale_results['status'] = False
					logger.info('Error: %s does not exist' %filepath)
			elif os.path.isfile(path) and path.endswith('.HKL'):
				self.hklpaths.append(path)
				self.scale_results['status'] = True
			else:
				err = "Incorrect format of files, couldn't find XDS_ASCII.HKLs"
				logger.info('Error:{}'.format(err))
				self.scale_results['status'] = False
				return self.scale_results
		return self.scale_results

	
	def create_file_links(self):
		if len(self.hklpaths) == 0:
			err = IndexError("no hkl filenames exist")
			logger.info('IndexError: {}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status;
			return self.scale_results
		for ii in range(len(self.hklpaths)):
			try:
				linkname = "xtal_"+str(ii)+".HKL"

				os.symlink(self.hklpaths[ii], linkname)
				self.filelinks.append(linkname)
				self.status = True; #self.scale_results['status'] = self.status
			except OSError, e:
				if e.errno == errno.EEXIST:
					os.remove(linkname)
					os.symlink(self.hklpaths[ii], linkname)
					self.filelinks.append(linkname)
					self.status = True
				else:
					logger.info('Error:{}'.format(e))
					self.status = False; self.scale_results['status'] = self.status
					self.scale_results['info'].append(e)
		return self.scale_results

	def indexing_(self):
		if  self.reference is None:
			self.reference = self.hklpaths[0]

		try:
			for ii in range(1,len(self.hklpaths)):
				_correct = index_check.is_correct(self.hklpaths[ii], self.reference, self.output, self.user)
				if _correct == True:
					logger.info("indexing is consistent\n")
				else:
					logger.info("wrong indexing\n")
					self.hklpaths.pop(ii)
				self.status = True
				self.scale_results['status'] = self.status


		except IndexError, ValueError:
			pass
		logger.info('MSG: # of cleaned xtals %s' %len(self.hklpaths))
		self.scale_results['xtals_after_idx_check'] = len(self.hklpaths)
		return self.status


	def create_inp(self, fList, **kwargs):
		a = kwargs.get('a', 79.0)
		b = kwargs.get('b', 79.0)
		c = kwargs.get('c', 39.0)
		al = kwargs.get('al', 90.0); be = kwargs.get('be', 90.0); ga = kwargs.get('ga', 90.0)
		reso = kwargs.get('reso',2.2); refs = kwargs.get('refs', None)
		cellkeys = ['a', 'b', 'c', 'al', 'be', 'ga']
		try:
			bl = os.environ['BEAMLINE_XNAME']
			if bl =='X06SA':
				self.nproc = 30
			elif bl == 'X06DA' or bl == 'X10SA':
				self.nproc = 20
		except KeyError:
			MSG = "XSCALE running locally, not from beamline\n"
			logger.info('MSG:{}'.format(MSG))
			self.nproc = 4

		fh = open("XSCALE.INP",'w')
		fh.write("OUTPUT_FILE=XSCALE.HKL\n")
		fh.write("MAXIMUM_NUMBER_OF_PROCESSORS=%d\n" %self.nproc)
		fh.write("SAVE_CORRECTION_IMAGES=FALSE\n")
		fh.write("FRIEDEL'S_LAW=FALSE\n\n")
		if all(k in kwargs.keys() for k in cellkeys):
			fh.write('UNIT_CELL_CONSTANTS= %f  %f  %f  %f  %f %f\n' %(a, b, c, al, be, ga))
		else:
			pass
		if 'refs' in kwargs.keys():
			try:
				os.symlink(refs, 'reference.HKL')
			except OSError, e:
				if e.errno == errno.EEXIST:
					os.remove('reference.HKL')
					os.symlink(refs, 'reference.HKL')
				else:
					logger.info('Error:{}'.format(e))
					self.scale_results['info'].append(e)
					self.scale_results['status'] = False
			fh.write('REFERENCE_DATA_SET= reference.HKL\n')

		if 'reso' in kwargs.keys():
			for f in fList:
				fh.write("INPUT_FILE=%s\n" %f)
				fh.write("MINIMUM_I/SIGMA=0.0\n")
				fh.write("INCLUDE_RESOLUTION_RANGE= 50 %f\n" %reso)
		else:
			for f in fList:
				fh.write("INPUT_FILE=%s\n" %f)
				fh.write("MINIMUM_I/SIGMA=0.0\n")
		fh.close()

	def Isa_select(self, ISa_th):
		#method to perform ISa based selection of 'good' HKL files for next round of XSCALEing
		if os.path.isfile("XSCALE.LP"):
			fh = open("XSCALE.LP", 'r')
			all_lines = fh.readlines()
			fh.close()
			try:
				start = all_lines.index('     a        b          ISa    ISa0   INPUT DATA SET\n')
				end = all_lines.index(' (ASSUMING A PROTEIN WITH 50% SOLVENT)\n')
				# select the table with ISa values from the file..
				Isa_list = all_lines[start+1:end-3]
				self.status = True
			except (ValueError, IndexError) as err:
				error = err("check if XSCALE ran properly or XSCALE.INP screwed up \n")
				logger.info('Error: {}'.format(error))
				self.scale_results['info'].append(error)
				self.status = False; self.scale_results['status']=self.status
				return self.scale_results

			try:
				for lines in Isa_list:
					line = lines.split()
					try:
						if float(line[2]) > float(ISa_th):
							self.selected_files.append(line[4])
					except IndexError:
						pass
				self.selected_files = self.sort_xtal(self.selected_files)
				self.create_inp(self.selected_files, reso=3.0)
				msg = "%d xtals selected by ISa-cut out of %d xtals\n" %(len(self.selected_files), len(self.hklpaths))
				self.scale_results['xtals_after_isa'] = len(self.selected_files)

				logger.info('MSG:{}'.format(msg))
				try:
					run_command("Scale&Merge", self.subadm, self.user, self.xscale_cmd, 'merge.log')
				except Exception:
					sub.call(self.xscale_cmd, shell=True)
					self.status = True
			except Exception as e:
				logger.info('Error:{}'.format(e))
				self.status = False; self.scale_results['status'] = self.status
				self.scale_results['info'].append(e)
				return self.scale_results
		else:
			err = "XSCALE.LP file does not exist"
			logger.info('IOError: {}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results
		return self.scale_results

	def xscale_for_sad(self):
		self.find_HKLs()
		try:
			self.outdir_make()
			os.chdir(self.subadm)
		except OSError:
			err = OSError("multi folder either not there or not accessible\n")
			logger.info('OSError:{}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
		if len(self.hklpaths) > 0:
			if len(self.hklpaths) < 50:
				self.indexing_()
			else:
				logger.info('Too many hkls, so skipping the indexing check')
				pass
			self.create_file_links()
			state, self.reference = ref_choice(self.hklpaths, fom='bfac')
			if state == False and self.reference == None:
				self.status = False; self.scale_results['status'] = self.status;
				logging.info('Error: bfactor-based reference choice did not work')
				
			elif state == True and self.reference != None:
				logging.info('lowest-Bfactor file: %s' %self.reference)
				self.create_inp(self.filelinks, refs=self.reference, reso=3.0)
			else:
				self.create_inp(self.filelinks)

			msg = "xscale-ing of native SAD data\n"
			logger.info('MSG:{}'.format(msg))
			try:
				run_command("Scale&Merge", self.subadm, self.user, self.xscale_cmd, 'merge.log')
			except Exception:
				sub.call(self.xscale_cmd, shell=True)

			try:
				state, stats = parse_xscale_output("XSCALE.LP")
				if state == False:
					self.status = False; self.scale_results['status'] = self.status
					return self.scale_results
				else:

					logger.info('stat_dict:{}'.format(stats))
					self.scale_results['nSAD_xscale_stats'] = stats
					shutil.copyfile('XSCALE.HKL', 'nSAD.HKL')
					shutil.copyfile('XSCALE.LP', 'nSAD.LP')
					self.status = True
			except Exception as err:
				logger.info('Error: {}'.format(err))
				self.scale_results['info'].append(err)
				self.status = False; self.scale_results['status'] = self.status
			return self.scale_results

	def xscale_for_sx(self):
		self.find_HKLs()
		try:
			self.outdir_make()
			os.chdir(self.subadm)
		except OSError:
			err = OSError("adm folder either not there or not accessible\n")
			logger.info('OSError:{}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results

		try:
			if len(self.hklpaths) < 100:
				self.indexing_()
			else:
				MSG = "Too many hkls, so skipping indexing check"
				logger.info('MSG:{}'.format(MSG))
				pass
			self.create_file_links()
			state, self.reference = ref_choice(self.hklpaths, fom='bfac')
			if state == False and self.reference == None:
				self.status = False; self.scale_results['status'] = self.status;
				logging.info('Error: bfactor-based reference choice did not work')
				
			elif state == True and self.reference != None:
				logging.info('lowest-Bfactor file: %s' %self.reference)
				self.create_inp(self.filelinks, refs=self.reference, reso=3.0)
			else:
				self.create_inp(self.filelinks)

			#self.create_inp(self.filelinks)

			msg = "Running 1st round of xscale-ing\n"
			logger.info('MSG:{}'.format(msg))
			try:
				run_command("Scale&Merge", self.subadm, self.user, self.xscale_cmd, 'merge.log')
			except Exception:
				sub.call(self.xscale_cmd, shell=True)
			try:
				state, stats = parse_xscale_output("XSCALE.LP")
				if state == False:
					self.status = False; self.scale_results['status'] = self.status
					return self.scale_results

				logger.info('stat_dict:{}'.format(stats))
				self.scale_results['no_selection'] = stats
				shutil.copyfile("XSCALE.LP", "noSelect.LP")
				shutil.copyfile("XSCALE.HKL", "noSelect.HKL")
				self.status = True; self.scale_results['status'] = self.status
			except OSError:
				err = OSError("xscaling may not have run\n")
				logger.info('OSError:{}'.format(err))
				self.scale_results['info'].append(err)
				self.status = False; self.scale_results['status'] = self.status
				return self.scale_results

			celler = Cell('noSelect.LP')
			celler.cell_histogram()
			if celler.status == False:
				self.status = False
				return self.scale_results

			mode_cells = {'a': celler.a_mode, 'b': celler.b_mode, 'c': celler.c_mode,
						'al':celler.al_mode, 'be':celler.be_mode, 'ga': celler.ga_mode}

			self.scale_results['unit-cell'] = mode_cells;
			self.scale_results['cell_array'] = celler.cell_ar
			self.scale_results['cell_array_clean'] = [celler.a_ohne, celler.b_ohne, celler.c_ohne, celler.al_ohne, celler.be_ohne, celler.ga_ohne]
			try:
				#self.create_inp(self.sort_xtal(celler.cell_select), **mode_cells)
				self.create_inp(celler.cell_select, reso=3.0, **mode_cells)
				msg = "xscaling after cell-analysis and rejecting outliers\n"
				logger.info('MSG:{}'.format(msg))
				try:
					run_command("Scale&Merge", self.subadm, self.user, self.xscale_cmd, 'merge.log')
				except (OSError, TypeError, Exception) as e :
					sub.call(self.xscale_cmd, shell=True)

				msg = "%d xtals got selected after cell analysis out of %d xtals" %(len(celler.cell_select), len(self.hklpaths))
				logger.info('MSG:{}'.format(msg))
				self.status = True
				self.scale_results['xtals_after_cell'] = len(celler.cell_select)
				self.scale_results['status'] = self.status
			except Exception as e:
				logger.info('Error:{}'.format(e))
				self.scale_results['info'].append(e)
				self.status = False; self.scale_results['status'] = self.status
				return self.scale_results
			try:
				state, stats = parse_xscale_output("XSCALE.LP")
				if state == False:
					self.status = False; self.scale_results['status'] = self.status
					return self.scale_results

				logger.info('stat_dict:{}'.format(stats))
				self.scale_results['cell_selection'] = stats
				shutil.copyfile("XSCALE.LP", "Cell_Select.LP")
				shutil.copyfile("XSCALE.HKL", "Cell_Select.HKL")
				self.status = True; self.scale_results['status'] = self.status
			except OSError:
				err = OSError("xscaling after Cell selection may not work\n")
				logger.info('OSError:{}'.format(err))
				self.scale_results['info'].append(err)
				self.status = False; self.scale_results['status'] = self.status
				return self.scale_results

			msg = "running xscale after ISa selection\n"
			logger.info('MSG:{}'.format(msg))
			self.Isa_select(self.isa_cut)

			try:
				if self.status == False:
					self.scale_results['status'] = self.status
					return self.scale_results
				state, stats = parse_xscale_output("XSCALE.LP")

				logger.info('stat_dict:{}'.format(stats))
				self.scale_results['ISa_selection'] = stats
				shutil.copyfile("XSCALE.LP", "ISa_Select.LP")
				shutil.copyfile("XSCALE.HKL", "ISa_Select.HKL")
				self.status = True; self.scale_results['status'] = self.status
			except OSError:
				err = OSError("xscaling after ISa selection may not work\n")
				logger.info('OSError:{}'.format(err))
				self.scale_results['info'].append(err)
				self.status = False; self.scale_results['status'] = self.status
				return self.scale_results

			cc_dict, b_dict, datasets = get_cc_error_b('ISa_Select.LP')
			self.hkl_cc_sorted = pair_corr_sorter(cc_dict, datasets)

			msg = "pair-correlation sorting produced: %d hkls" %len(self.hkl_cc_sorted)
			logger.info('MSG:{}'.format(msg))
			
			self.create_inp(self.hkl_cc_sorted, reso=3.0, **mode_cells)
			try:
				run_command("Scale&Merge", self.subadm, self.user, self.xscale_cmd, 'merge.log')
			except (OSError, TypeError, Exception) as e :
				sub.call(self.xscale_cmd, shell=True)

			self.status = True
			self.scale_results['xtals_after_pCC'] = len(self.hkl_cc_sorted)
			self.scale_results['status'] = self.status
			try:
				state, stats = parse_xscale_output("XSCALE.LP")
				if state == False:
					self.status = False; self.scale_results['status'] = self.status
					return self.scale_results

				logger.info('stat_dict:{}'.format(stats))
				self.scale_results['pCC_selection'] = stats
				shutil.copyfile("XSCALE.LP", "pCC_Select.LP")
				shutil.copyfile("XSCALE.HKL", "pCC_Select.HKL")
				self.status = True; self.scale_results['status'] = self.status
			except OSError:
				err = OSError("xscaling after Cell selection may not work\n")
				logger.info('OSError:{}'.format(err))
				self.scale_results['info'].append(err)
				self.status = False; self.scale_results['status'] = self.status
				return self.scale_results
		except Exception as e:
			logger.info('Error: {}'.format(e))
			self.scale_results['info'].append(e)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results
		
		return self.scale_results

	def isocluster(self, xscalefile):

		if os.path.isfile(xscalefile):
			comd = "xscale_isocluster %s > isocluster.log" %xscalefile
			msg = "iso-clustering based on Correlation\n"
			logger.info('MSG:{}'.format(msg))
			try:
				run_command("Scale&Merge", self.subadm, self.user, comd, 'merge.log')
			except Exception:
				sub.call(comd, shell=True)

				self.status = True; self.scale_results['status'] = self.status
		else:
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results

		if os.path.isfile("isocluster.log"):
			fkey = "best number of clusters"
			fh = open("isocluster.log", "r")
			_all = fh.readlines(); fh.close()
			for lines in _all:
				if "Error" in lines:
					self.status = False
					self.scale_results['status'] = self.status
					return self.scale_results
				elif fkey in lines:
					line = lines.split(':')
					val = line[-1].strip('\n')
					self.scale_results["clusters"] = val.strip(' ')
					self.status = True;
					self.scale_results['status'] = self.status
				else:
					pass
		if self.status == True and os.path.isfile("XSCALE.1.INP"):
			shutil.copyfile("XSCALE.1.INP", "XSCALE.INP")
			msg = "xscale-ing over 1st cluster only..\n"
			logger.info('MSG:{}'.format(msg))
			try:
				run_command("Scale&Merge", self.subadm, self.user, self.xscale_cmd, 'merge.log')
			except Exception:
				sub.call(self.xscale_cmd, shell=True)
			try:
				state, stats = parse_xscale_output("XSCALE.LP")
				logger.info('stat_dict:{}'.format(stats))
				shutil.copyfile("XSCALE.LP", "cluster1.LP")
				shutil.copyfile("XSCALE.HKL", "cluster1.HKL")
				self.scale_results['iso-cluster'] = stats
				self.status = True; self.scale_results['status'] = self.status
			except OSError:
				err = OSError("xscaling after iso-clustering may not work\n")
				logger.info('OSError:{}'.format(err))
				self.scale_results['info'].append(err)
				self.status = False; self.scale_results['status'] = self.status
				return self.scale_results
		return self.scale_results

	def aniso_check(self):
		point_cmd = "pointless -xdsin noSelect.HKL hklout noSelect_point.mtz > /dev/null"
		try:
			run_command("Scale&Merge", self.subadm, self.user, point_cmd, 'merge.log')
		except Exception:
			sub.call(point_cmd, shell=True)


		if os.path.isfile("noSelect_point.mtz"):
			batch_cmd = "phenix.python /exchange/mx/daplus/imisx_processing/mergetools/get_batch.py noSelect_point.mtz"
			try:
				batches = run_command('Scale&Merge', self.subadm, self.user, 'merge.log')

				batches = batches.strip('[').strip(']').split(",")
			except Exception:
				batches = sub.check_output(batch_cmd, shell=True)
				batches = batches.strip('[').strip(']').split(",")

			fh = open("onlymerge.inp", 'w')
			fh.write("ONLYMERGE\n")
			fh.write("run 1 batch %s to %s\n" %(batches[0], batches[100]))
			fh.close()
			aim_cmd = "aimless hklin noSelect_point.mtz hklout noSelect_aim.mtz < onlymerge.inp > aimless.log"
			try:
				run_command("Scale&Merge", self.subadm, self.user, aim_cmd, 'merge.log')
			except (OSError, TypeError, Exception) as e:
				sub.call(aim_cmd, shell=True); self.status = True
				self.scale_results['status'] = self.status
		else:
			err = OSError("pointless has not run properly, check if HKL file is correct \n")
			logger.info('OSError:{}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results
		try:
			fh = open("aimless.log", 'r')
			_all = fh.readlines()
			fh.close(); self.status = True;
			self.scale_results['status'] = self.status
		except OSError, IOError:
			err = OSError("aimless.log file doesn't exist")
			logger.info('OSError:{}'.format(err))
			self.scale_results['info'].append(err)
			self.status = False; self.scale_results['status'] = self.status
			return self.scale_results

		keyword = "Estimated maximum resolution limits"

		for lines in _all:

			if "Error" in lines:
				self.status = False
				self.scale_results['status'] = self.status
				return self.scale_results
			elif keyword in lines:
				line = lines.split(',')
				a = line[1]; b = line[2]; c = line[3].strip('\n')
				self.scale_results['anisotropicity'] = (a, b, c)
				logger.info('Anisotropy:{}'.format((a,b,c)))
				self.status = True; self.scale_results['status'] = self.status
			else:
				pass
		return self.scale_results

	def run_(self):
		
		try:
			if self.expt == 'native-sad':
				self.xscale_for_sad()
				if len(self.hklpaths) > 10:
					self.isocluster('XSCALE.HKL')
				return self.scale_results
			elif self.expt == 'serial-xtal':
				self.xscale_for_sx()
				self.isocluster('XSCALE.HKL')
				self.aniso_check()
				return self.scale_results
			elif self.expt == 'inverse-beam' or self.expt == 'interleave-and-inverse-first':
				pass
			else:
				logger.info('Error: experiment type not supported %s' %self.expt)
				self.scale_results['info'].append('Expt. type not supported')
				self.scale_results['status'] = False
		except Exception as e:
			self.scale_results['status'] = False
			logger.info('Error:{}'.format(e))
			self.scale_results['info'].append(e)
			return self.scale_results
		return self.scale_results


def finder():
	root = sys.argv[1]; expt = sys.argv[2];
	if expt == 'serial-xtal':
		paths = glob.glob(os.path.join(root, 'minisets_*/xtal_*/XDS_ASCII.HKL'))
	elif expt == 'native-sad':
		paths = glob.glob(os.path.join(root, 'set*'))
	return paths

def get_paths_xscale():
	xscale_file = sys.argv[1]; #expt = sys.argv[2];
	paths = [];
	if not xscale_file.endswith('LP'):
		msg = "TypeError: .LP file needed"
		logger.info(msg)
	else:
		fh = open(xscale_file, 'r')
		_all = fh.readlines()
		fh.close()
		for lines in _all:
			if 'INPUT_FILE' in lines:
				line = lines.split('=')
				abs_path = os.path.abspath(line[1].strip('\n'))
				paths.append(abs_path)
			else:
				pass
	return paths

def _run_(hklpath, expt_type, username=None):
	success = False; merge_result = {};
	try:
		scale = Merge_utls(hklpath, expt_type, user=username)
		if scale.expt == 'native-sad':
			xscale_out = scale.xscale_for_sad()
			success = True; merge_result['nSAD_xscale'] = xscale_out
			cluster = scale.isocluster('XSCALE.HKL')
			merge_result['isocluster'] = cluster

		elif scale.expt == 'serial-xtal':
			xscale_out = scale.xscale_for_sx()
			merge_result['imisx_xscale'] = xscale_out
			cluster = scale.isocluster("ISa_Select.HKL")
			merge_result['isocluster'] = cluster
			aniso = scale.aniso_check()
			merge_result['aniso_result'] = aniso
			success = True
		elif scale.expt == 'inverse-beam' or scale.expt == 'interleave-and-inverse-first':
			pass
		else:
			logger.info('Error: experiment type not supported %s' %scale.expt)
			success = False

	except Exception as e:
		success = False
		logger.info('Error:{}'.format(e))
		return success, merge_result
	return success, merge_result


if __name__ == '__main__':
	logging.basicConfig(level=logging.DEBUG,
	format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
	datefmt='%y-%m-%d %H:%M',
	filename='merge.log',
	filemode='a+')
	#hklpath = finder()
	hklpath = get_paths_xscale()
	
	scale = Merge_utls(hklpath, sys.argv[2])
	results = scale.run_()
	with open('example.dict', 'w') as ofh:
		for k, v in results.iteritems():
			ofh.write("%s: \t %s\n" %(k,v))
	ofh.close()

	#_run_('SX', 'e14365')
