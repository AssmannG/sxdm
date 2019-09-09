'''
Date: Aug, 2017
S.Basu
'''


import os, sys, errno
import glob, time
import subprocess as sub
import h5py
import numpy as np
import multiprocessing as mp
import argparse, logging
import xds_input

import Merge_utls as merge
from cellprobe import Cell
from xscale_output import *
import matplotlib.pyplot as plt


class Xtal(object):

	def __init__(self,xtalImgPath,xtalProcessPath,xtalNum, BL, tot_angle=360, **kwargs):


		self.xtalimgpath = xtalImgPath
		self.xtalprocesspath = xtalProcessPath
		self.xtalnum = xtalNum
		self.xtalname = None
		self.osci = None
		self.osci_range = tot_angle
		self.datrange_str = None
		self.bgrange_str = None
		self.sprange_str = None
		self.beamline = BL
		self.SG = kwargs.get('SG', '0')
		self.cell = kwargs.get('cell', "70 70 30 90 90 90")
		self.res_cut = kwargs.get('res_cut', "2.5")
		self.idx_res = kwargs.get('idx_res', "5.0")
		self.friedel = kwargs.get('friedel', "FALSE")
		self.refdata = kwargs.get('refdata', " ")
		self.strong_pixel = kwargs.get('strong_pixel', '6.0')
		self.min_pix_spot = kwargs.get('min_pix_spot', '3')
		self.content = {}

	def read_masterh5(self, master_file):

		#read master file headers for xds.inp preparation
		header = h5py.File(master_file, 'r')
		#beam center x and y
		beamx = header['/entry/instrument/detector/beam_center_x']
		beamx = np.array(beamx);
		self.content['beamX'] = beamx
		beamy = header['/entry/instrument/detector/beam_center_y']
		beamy = np.array(beamy)
		self.content['beamY'] = beamy

		#wavelength and detector distance
		wave = header['/entry/instrument/beam/incident_wavelength']
		wave = np.array(wave)
		self.content['wavelength'] = wave
		detZ = header['/entry/instrument/detector/detector_distance']
		detZ = np.array(detZ)
		self.content['detectorDistance'] = round(detZ*1e3) # convert distance into millimeter

		#omega and oscillation
		try:
			omega = header['/entry/sample/goniometer/omega_increment']
			omega = np.array(omega)
			self.osci = omega
			self.content['oscillationAngle'] = omega
		except KeyError:
			self.content['oscillationAngle'] = 0.1
		#number of images
		if self.osci_range is None:
		   nimages = header['/entry/instrument/detector/detectorSpecific/nimages']
		   nimages = np.array(nimages)
		else:
		   nimages = round(float(self.osci_range)/self.osci)
		   nimages = int(nimages)

		data_start = 1; data_end = nimages
		#data_start = self.xtalnum*nimages + 1; data_end = (self.xtalnum+1)*nimages
		self.datrange_str = (str(data_start), str(data_end))
		self.sprange_str = self.datrange_str
		self.bgrange_str = (str(data_start), str(data_start+10))
		self.content['firstIndex'] = data_start; self.content['numFrames'] = nimages;
		self.content['lastIndex'] = data_start+10

		# xtal filename template with path
		name = os.path.basename(master_file)
		name = name.split( 'master' )
		#self.xtalname = self.xtalimgpath + name[0]+"??????.h5"
		img = name[0]+"??????.h5"
		self.xtalname = os.path.join(self.xtalimgpath, img)
		self.content['xtalname'] = self.xtalname;

	def read_cbf(self, headerfile):
		#read cbf file header and store the info in dictionary and later prepare xds.inp

		cmd = "head -35 "+headerfile+" > filehead.txt"
		sub.call(cmd, shell=True)
		self.xtalname = headerfile[:-9]+"?????.cbf"
		self.content['xtalname'] = self.xtalname
		keys_head = ["Beam_xy","Wavelength", "Detector_distance","Angle_increment"]
		fh = open('filehead.txt', 'r')
		all_lines = fh.readlines()
		fh.close()

		for lines in all_lines:
			if any(key in lines for key in keys_head):
				line = lines.split()
				if line[1] == 'Beam_xy':
					self.content['beamX'] = str(line[2].strip('(').strip(','))

					self.content['beamY'] = str(line[3].strip(')'))
				elif line[1] == 'Wavelength':
					self.content['wavelength'] = line[2]
				elif line[1] == 'Detector_distance':
					self.content["detectorDistance"] = str(float(line[2])*1e3)
				else:
					self.content["oscillationAngle"] = line[2]

					self.content['numFrames'] = int(int(self.osci_range)/float(line[2]))
		self.content['firstIndex'] = 1; self.content['lastIndex'] = 11;


	def locatextalpath(self):
		#Locate and sanity check if the data exists and then read the headers/master files

		if not os.path.exists(self.xtalimgpath):
			print 'Error: path does not exist\n'
			sys.exit()
		if self.beamline == "PXI":
			master_file = glob.glob(os.path.join(self.xtalimgpath,"*_master.h5"))[0]
			self.content['lib'] = '/exchange/mx/xds/library/dectris-neggia-centos6.so'
			try:
				self.read_masterh5(master_file)
			except OSError:
				raise OSError("master file may not exist\n")
		elif self.beamline == "PXII" or self.beamline == "PXIII":
			cbf_header = glob.glob(os.path.join(self.xtalimgpath,"*_00001.cbf"))[0]
			try:
				self.read_cbf(cbf_header)
			except OSError:
				raise OSError("cbf may not have collected yet\n")

	def create_idx_inp(self):
		try:

			os.chdir(self.xtalprocesspath)
		except OSError:
			raise OSError("xtal process folder have not been created yet")

		self.content['jobs']='XYCORR INIT COLSPOT IDXREF'
		self.content['njobs'] = 4
		try:
			if os.environ['BEAMLINE_XNAME'] == 'X06SA':
				self.content['nproc'] = 18
				self.content['nodes'] = "x06sa-cn-117 x06sa-cn-118 x06sa-cn-119 x06sa-cn-120 x06sa-cn-121 x06sa-cn-122 x06sa-cn-123 x06sa-cn-124"
			elif os.environ['BEAMLINE_XNAME'] == 'X06DA':
				self.content['nproc'] = 12
				self.content['nodes'] = "x06da-cn-1 x06da-cn-2"
			elif os.environ['BEAMLINE_XNAME'] == 'X10SA':
				self.content['nproc'] = 12
				self.content['nodes'] = "x10sa-cn-1 x10sa-cn-2"

		except KeyError:
			if 'SLURM_NODELIST' in os.environ:
				node_string = os.environ['SLURM_NODELIST']
				node_num_list = node_string.strip('ra-c-[').strip(']').split('-')
				self.content['nproc'] = 12
				self.content['njobs'] = 2
				self.content['nodes'] = 'ra-c-%s ra-c-%s' %(node_num_list[0], node_num_list[1])
				if self.beamline == 'PXI':
					sub.call(['module load dectris-neggia/17.09'], shell=True)
					self.content['lib'] = os.path.join(os.environ['DECTRIS_NEGGIA_LIBRARY_DIR'], 'dectris-neggia.so')
			else:
				print 'On Ra-cluster, salloc was not done so using logging node, will be slow\n'
				self.content['nproc'] = 12
				self.content['njobs'] = 1
				self.content['nodes'] = 'ra-c-001 ra-c-002 ra-c-003 ra-c-0004'
				if self.beamline == 'PXI':
					sub.call(['module load dectris-neggia/17.09'], shell=True)
					self.content['lib'] = os.environ['DECTRIS_NEGGIA_LIBRARY_DIR']

		if self.refdata != " " and os.path.isfile(self.refdata):
			ref_link = 'reference.HKL'
			os.symlink(self.refdata, ref_link)
			self.content['referenceData'] = ref_link
		else:
			self.content['referenceData'] = " "

		self.content['SG'] = self.SG; self.content['unit_cell'] = self.cell;
		self.content['friedel'] = self.friedel; self.content['highres'] = 5.0
		self.content['strong_pixel'] = self.strong_pixel; self.content['min_pix_spot'] = self.min_pix_spot;

		inp_string = xds_input.INP[self.beamline]
		if not os.path.isfile("XDS.INP"):
			fh = open("XDS.INP", 'w')
			fh.write(inp_string.format(**self.content))
			fh.close()
		else:
			pass

	def create_integrate_inp(self):
		try:
			os.chdir(self.xtalprocesspath)
		except OSError:
			raise OSError('xtal process folder may not be created yet')
		self.content['jobs'] = 'DEFPIX INTEGRATE CORRECT'
		self.content['highres'] = self.res_cut

		if os.path.isfile("XDS.INP"):
			sub.call(["cp XDS.INP indexing.INP"], shell=True)
			inp_string = xds_input.INP[self.beamline]
			fh = open("XDS.INP",'w')
			fh.write(inp_string.format(**self.content))
			fh.close()
		else:
			pass

	def check_allfiles(self):
		if self.beamline == 'PXII' or self.beamline == 'PXIII':
			if len(str(self.content['numFrames'])) == 3:
				tmp_str = '00'+str(self.content['numFrames'])
			else:
				tmp_str = '0'+str(self.content['numFrames'])
			lastImage = self.content['xtalname'].strip('?????.cbf')+tmp_str+'.cbf'
			wait_max = self.content['numFrames']*10
			wait = 0;
			while not os.path.exists(lastImage):
				time.sleep(10)
				print("waiting for the last image: %s" %lastImage)
				wait += 5
				if wait > wait_max:
					print "all images were not saved, so processing timed out\n"
					break

		else:
			pass


	def runxds(self):
		self.check_allfiles()
		os.chdir(self.xtalprocesspath)
		sub.call(['xds_par > /dev/null'], shell=True)
		self.create_integrate_inp()
		sub.call(['xds_par > /dev/null'], shell=True)
		if not os.path.isfile("XDS_ASCII.HKL"):
			print "xtal: %d failed from %s\n" %(self.xtalnum, self.xtalprocesspath),
		else:
			print "xtal: %d processed\n" %self.xtalnum,


class Process(object):
	"""docstring for ClassName"""
	def __init__(self, data_dir, output, BL, tot_angle=360):
		self.data_folder = data_dir
		self.process_folder = None
		self.output = os.path.join(output, "proc")
		if not os.path.exists(self.output):
			os.makedirs(self.output, 0755)
		else:
			pass
		self.nxtals = None
		self.setname = None
		self.setfolder = [];
		self.dataName = None
		self.process_folder = None
		self.process_data = None
		self.beamline = BL
		self.ISa_th = 4.0
		self.total_range = tot_angle
		self.xscale_file_list = [];
		self.xtals_lists = []


	def Lookup(self):
		if len(self.data_folder) == 0:
			print "No image folder found \n"
			return
		for ii in range(len(self.data_folder)):
			for dirs in sorted(glob.glob(os.path.join(self.data_folder[ii], "*set*"))):
				if os.path.isdir(dirs) and len(os.listdir(dirs)) > 2:
					self.setfolder.append(dirs)
		return


	def get_xtals(self, **kwargs):
		for ii in range(len(self.data_folder)):
			if self.data_folder[ii].endswith('*'):
				parent_dir = self.data_folder[ii][:-1]
				if not os.path.exists(parent_dir):
					print "Error: data directory does not exist!\n"
					sys.exit()
			else:
				if not os.path.exists(self.data_folder[ii]):
					print "Error: data directory does not exist!\n"
					sys.exit()
		try:
			os.chdir(self.output)
			self.Lookup()
		except OSError:
			raise IOError("output path is not accessible\n")
			sys.exit()
		self.nxtals = len(self.setfolder)
		if self.nxtals > 0:
			for k in range(self.nxtals):
				self.setname = os.path.basename(self.setfolder[k])
				dir_von_sets = os.path.dirname(self.setfolder[k])
				self.dataName = os.path.basename(dir_von_sets)
				self.process_folder = os.path.join(self.output, self.dataName, self.setname)
				self.process_data = os.path.join(self.output, self.dataName)

				if not os.path.isdir(self.process_folder):
					print "creating processing directory %s\n" %(self.process_folder),
					os.makedirs(self.process_folder, 0755)
					os.chdir(self.process_folder)
					xtalobj = Xtal(self.setfolder[k], self.process_folder,k, self.beamline, self.total_range, **kwargs)
					xtalobj.locatextalpath()
					xtalobj.create_idx_inp()
					self.xtals_lists.append(xtalobj)
				else:
					print "folder may exist, skipping %s\n" %(self.process_folder),
				os.chdir(self.output)

		print "%d xtals have been found \n" %len(self.xtals_lists)
		return

	def get_serial_eiger_xtals(self, **kwargs):
		for ii in range(len(self.data_folder)):
			if self.data_folder[ii].endswith('*'):
				parent_dir = self.data_folder[ii][:-1]
				if not os.path.exists(parent_dir):
					print "Error: data directory does not exist!\n"
					sys.exit()
			else:
				if not os.path.exists(self.data_folder[ii]):
					print "Error: data directory does not exist!\n"
					sys.exit()
		try:
			os.chdir(self.output)
			self.Lookup()
		except OSError:
			raise IOError("output path is not accessible\n")
			sys.exit()
		if len(self.setfolder) > 0:
			for i in range(len(self.setfolder)):
				self.xtal_each_miniset = sorted(glob.glob(os.path.join(self.setfolder[i], "*data*.h5")))
				self.nxtals = len(self.xtal_each_miniset) #num of xtals in each minisets
				print "%d xtals in %s miniset \n" %(self.nxtals, self.setfolder[i]),
				self.mininame = os.path.basename(self.setfolder[i])
				dir_von_miniset = os.path.dirname(self.setfolder[i])
				self.dataName = os.path.basename(dir_von_miniset)
				puckname = os.path.basename(os.path.dirname(self.dataName))
				self.process_data = os.path.join(self.output, puckname, self.dataName)
				self.process_folder = os.path.join(self.output, puckname, self.dataName, self.mininame)
				if not os.path.exists(self.process_folder):
					print "creating processing directory %s" %(self.process_folder)
					os.makedirs(self.process_folder, 0755)
				else:
					pass
				os.chdir(self.process_folder)

				for k in range(self.nxtals):
					xtal_process_path = self.process_folder + '/xtal_' + str(k)
					if not os.path.exists(xtal_process_path):
						os.makedirs(xtal_process_path, 0755)
						os.chdir(xtal_process_path)
						xtalobj = Xtal(self.setfolder[i], xtal_process_path, k, self.beamline, self.total_range, **kwargs)
						xtalobj.locatextalpath()
						image_block = xtalobj.content['numFrames']
						start_image = xtalobj.xtalnum*image_block+1
						end_image = (xtalobj.xtalnum+1)*image_block
						xtalobj.content['firstIndex'] = start_image
						xtalobj.content['numFrames'] = end_image
						xtalobj.content['lastIndex'] = start_image+10
						xtalobj.create_idx_inp()
						self.xtals_lists.append(xtalobj)
					else:
						print "folder may exist, skipping it %s\n" %xtal_process_path,
				os.chdir(self.output)
		print "\n %d xtals have been gathered\n" %len(self.xtals_lists)
		return



	def find_HKLs(self, **kwargs):

		mergepaths = kwargs.get('mergepaths',[self.output])
		try:
			os.chdir(self.output)
		except OSError:
			print "check if the output folder exists\n"

		for path in mergepaths:
			for parent, dirs, files in os.walk(path):
				for fh in files:
					if fh == "XDS_ASCII.HKL":
						HKLpath = os.path.join(parent,fh)
						self.xscale_file_list.append(HKLpath)
					else:
						pass
		return

	def runeiger(self):
		job_cnt = 0
		if len(self.xtals_lists) > 0:
			for j in range(len(self.xtals_lists)):

				proc = [];
				for i in range (0,4):
					try:
						jobid = mp.Process(target=self.xtals_lists[(j*4)+i].runxds)
						proc.append(jobid)

					except IndexError:
						pass
				for p in proc:
					p.start()

				for p in proc:
					p.join()


			print "%d crystals have been attempted\n" %((j+1))

	def runpilatus(self, expt):
		job_cnt = 0
		if expt == 'native-sad':
			try:
				for j in range(len(self.xtals_lists)):
					self.xtals_lists[j].runxds()
					job_cnt += 1
				print "%d crystals have been attempted\n" %job_cnt
			except Exception:
				raise Exception("no xtals found to run xds or other error, check\n")
		elif expt == 'serial-xtal':
			try:
				for j in range(len(self.xtals_lists)):

					proc = [];
					for i in range (0,10):
						try:
							jobid = mp.Process(target=self.xtals_lists[(j*10)+i].runxds)
							proc.append(jobid)

						except IndexError:
							pass
					for p in proc:
						p.start()

					for p in proc:
						p.join()
				print "%d crystals have been attempted\n" %(j+1)

			except Exception:
				raise


def options():
	parser = argparse.ArgumentParser()
	parser.add_argument("--image_path", type=str, nargs='+', \
						help="provide path for each well, containing minisets folder or provide path for parent folder, containing all wells, e.g. your/path/to/parent/<well-id> or your/path/to/parent")
	parser.add_argument("--output_dir", type=str, \
						help="provide path where processing stuffs will be dumped using identical directory tree")
	parser.add_argument("--BeamID", type=str, help="Beamline ID needs to be specified, eg. PXI or PXII\n")
	parser.add_argument("--method", type=str, help="Mention either native-sad or serial-xtal\n")

	parser.add_argument("--total_degree", type=str, help="provide angular range to process.It's mutually exclusive with start/end_omega keywords")


	parser.add_argument("--SG_num", type=str, default="0", \
						help="optionally, Space-group number can be specified, default is 0")
	parser.add_argument("--cell", type=str, default="70 70 30 90 90 90", \
						help="optionally, unit-cell can be specified as 70 70 30 90 90 90; otherwise it will try to determine by itself")

	parser.add_argument("--highres", type=str, default="2.5", \
						help="optionally high-resolution limit can be given, default: 2.5")
	parser.add_argument("--friedel", type=str, default="FALSE", help="optionally, it can be changed to true..")

	parser.add_argument("--refs", type=str, help='optionally, reference data set for indexing can be provided..')

	parser.add_argument("--strong_pixel", type=str)
	parser.add_argument("--min_pix_spot", type=str)
	parser.add_argument("--ISa_cutoff", type=str, default= "3.0")
	parser.add_argument("--merge_paths", type=str, nargs='+')

	args = parser.parse_args()
	return args

def main():

	args = options()
	if args.image_path is None:
		sys.exit("you didn't tell me where data is\n")
	elif args.output_dir is None:
		print "no output path provided, dumping everything in current directory\n"
		args.output_dir = os.getcwd()

	elif args.BeamID is None:
		sys.exit("Beamline has to be mentioned, e.g. PXI, PXII, or PXIII\n")
	elif args.total_degree is None:
		args.total_degree = 360
	elif args.method is None:
		sys.exit("Please specify the method, either native-sad or serial-xtal\n")

	keywords = {}
	if args.SG_num != None:
		keywords['SG'] = args.SG_num
	if args.cell != None:
		keywords['cell'] = args.cell
	if args.highres != None:
		keywords['res_cut'] = args.highres
	if args.friedel != None:
		keywords['friedel'] = args.friedel
	if args.refs != None:
		keywords['refdata'] = args.refs
	if args.strong_pixel != None:
		keywords['strong_pixel'] = args.strong_pixel
	if args.min_pix_spot != None:
		keywords['min_pix_spot'] = args.min_pix_spot

	merge_keys = {}
	if args.ISa_cutoff != None:
		merge_keys['isa_cut'] = args.ISa_cutoff
	if args.highres != None:
		merge_keys['res_cut'] = args.highres

	merge_hkls = {}
	if args.merge_paths != None:
		merge_hkls['mergepaths'] = args.merge_paths[0].split()

	proc = Process(args.image_path, args.output_dir, args.BeamID, args.total_degree)
	if proc.beamline == "PXI" and args.method == 'serial-xtal':
		proc.get_serial_eiger_xtals(**keywords)
	else:
		proc.get_xtals(**keywords)

	if proc.beamline == "PXI":
		proc.runeiger()
	else:
		proc.runpilatus(args.method)

	#Merging with Merge_utls..
	proc.find_HKLs(**merge_hkls)


	mm = merge.Merge_utls(sorted(proc.xscale_file_list), args.method, **merge_keys)
	results = mm.run_()
	if args.method == 'native-sad':
		try:
			print "No selection table with xtals: %d\n" %results['xtals_found']
			print_table(results['nSAD_xscale_stats'])
			print'\nISa selection table with xtals: %d\n' %results['xtals_after_isa']
			print_table(results['ISa_selection'])
		except KeyError:
			print "Either ISa-select is not set or XSCALE stats not found\n"

	else:
		try:
			print "No selection table with xtals: %d\n" %results['xtals_found']
			print_table(results['no_selection'])
			print '\nISa selection table with xtals: %d\n' %results['xtals_after_isa']
			print_table(results['ISa_selection'])
			print '\nCell selection table with xtals: %d\n' %results['xtals_after_cell']
			print_table(results['cell_selection'])
			print '\npair-CC selection table with xtals: %d\n' %results['xtals_after_pCC']
			print_table(results['pCC_selection'])
			print '\nxscale_isocluster table from most populated cluster\n'
			print_table(results['iso-cluster'])
			print '\n\n'


			if len(proc.xscale_file_list) > 200:
				hkl_file_cell = os.path.join(mm.subadm, 'Cell_Select.LP')
				if os.path.isfile(hkl_file_cell):
					cell_hist = Cell(hkl_file_cell)
					cell_hist.cell_histogram()
			else:
				import scipy.cluster.hierarchy as sch
				fig = plt.figure()
				dn = sch.dendrogram(results['hclust_matrix'], p=10, labels=results['dendro_labels'], truncate_mode='level')
				fig.savefig('cell-dendrogram.png', dpi=300)
		except (KeyError, TypeError) as e:
			print "xscaleing had error, check \n"

	return

if __name__ == '__main__':
	main()
