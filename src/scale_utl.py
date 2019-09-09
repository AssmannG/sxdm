'''
created by S.Basu
On 28-July-2017
'''
import os, sys,glob
import logging
import subprocess as sub
import correct

logger = logging.getLogger('Scale&Merge')

def find_corrects(hklpaths):
	corr_paths = []
	for fname in hklpaths:
		folder = os.path.dirname(fname)
		path = os.path.join(folder, 'CORRECT.LP')
		if os.path.isfile(path):
			corr_paths.append(path)
			status = True
		else:
			logger.info('CORRECT.LP could not be found in %s' %folder)
			status = False
			return status, corr_paths
	return status, corr_paths

def check_bfactor(hklpaths):
	status, corr_paths = find_corrects(hklpaths)
	bfac_dicts = {}; sort_bfac_xasci = []
	if len(corr_paths) == 0:
		err = 'ValueError: no CORRECT.LP found'
		logger.info('ValueError: {}'.format(err))
		status = False;
		return status, sort_bfac_xasci
	else:
		for fname, cor_name in zip(hklpaths, corr_paths):
			fh = open(cor_name, 'r')
			_all = fh.readlines()
			fh.close()
			xasci = fname
			for lines in _all:
				if "WILSON LINE" in lines:
					line = lines.split()
					try:
						bfac_dicts[xasci] = float(line[9])
					except Exception:
						logger.info('B-factor might be negative, not considered')
				else:
					pass
			status = True

		sort_bfac_xasci = sorted(bfac_dicts.items(), key=lambda x : x[1])
	return status, sort_bfac_xasci

def rank_rmeas(hklpaths):
	status, corr_paths = find_corrects(hklpaths)
	rmeas_dict = {}; sort_rmeas_xasci = []
	if len(corr_paths) == 0:
		err = 'ValueError: no CORRECT.LP found'
		logger.info('ValueError: {}'.format(err))
		status = False;
		return status, sort_rmeas_xasci
	else:
		for fname, cor_name in zip(hklpaths, corr_paths):
			state, correct_stat = correct.parse_xds_stats(cor_name)
			if state == True and len(correct_stat) > 0:
				mean_rmeas = correct.mean_rmeas_calc(correct_stat)
				rmeas_dict[fname] = mean_rmeas
				status = True
			else:
				err = "%s could not be read, please check the file" %cor_name
				logger.info('Error: {}'.format(err))
				statue = False
				return status, sort_rmeas_xasci

		sort_rmeas_xasci = sorted(rmeas_dict.items(), key=lambda x:x[1])
	return status, sort_rmeas_xasci


def ref_choice(hklpaths, fom='bfac'):
	reference = None;
	if fom == 'bfac':
		status, sort_bfac_hkls = check_bfactor(hklpaths)
		if status == False and len(sort_bfac_hkls) == 0:
			return status, reference
		if status == True and len(sort_bfac_hkls) > 0:
			try:
				reference = sort_bfac_hkls[0][0]
			except IndexError, ValueError:
				err = 'bfactor selection may not work'
				logger.info('Error:{}'.format(err))
				status = False
				return status, reference
			return status, reference
	elif fom == 'rmeas':
		status, sort_rmeas_xasci = rank_rmeas(hklpaths)
		if status == False and len(sort_rmeas_xasci) == 0:
			return status, reference
		if status == True and len(sort_rmeas_xasci) > 0:
			try:
				reference = sort_rmeas_xasci[0][0]
			except IndexError, ValueError:
				err = 'Rmeas based referenceing may not have worked'
				logger.info('Error:{}'.format(err))
				status = False
				return status, reference
			return status, reference
	else:
		pass

def Bfact_sorter(hklpaths):
	bfac_sorted_hkls = [];
	status, sort_bfac_ascii = check_bfactor(hklpaths)
	if status == False and len(sort_bfac_ascii) == 0:
		return status, bfac_sorted_hkls
	elif status == True and len(sort_bfac_ascii) > 0:
		for i in range(len(sort_bfac_ascii)):
			bfac_sorted_hkls.append(sort_bfac_ascii[i][0])
		status = True
		return status, bfac_sorted_hkls
	else:
		status = False
		err = "Bfactor selection may have some problem,check"
		logging.info('Error:{}'.format(err))
		return status, bfac_sorted_hkls

def rmeas_sorter(hklpaths):
	rmeas_sorted_hkls = [];
	status, sort_rmeas_xasci = rank_rmeas(hklpaths)
	if status == False and len(sort_rmeas_xasci) == 0:
		return status, rmeas_sorted_hkls
	elif status == True and len(sort_rmeas_xasci) > 0:
		for i in range(len(sort_rmeas_xasci)):
			rmeas_sorted_hkls.append(sort_rmeas_xasci[i][0])
		status = True
		return status, rmeas_sorted_hkls
	else:
		status = False
		err = "Rmeas based sorting did not work, check"
		logger.info("Error:{}".format(err))
		return status, rmeas_sorted_hkls

def main():
	hklpaths = glob.glob(os.path.join(sys.argv[1], 'XDS_ASCII.HKL'))
	state, hkls = rmeas_sorter(hklpaths)
	print hkls

if __name__ == '__main__':
	main()
