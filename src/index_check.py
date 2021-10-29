
__author__ = ["S. Basu"]
__license__ = "M.I.T"
__date__ = "26/06/2017"
__refactordate__ = "07/05/2021"

import sys
import subprocess as sub
from src.ascii import ASCII
from src.run_command import *

logger = logging.getLogger('sxdm')

def pointless(fname, refname, dirname=None, user=None):

	success = False
	input_cmd = """pointless XDSIN %s -HKLREF %s << EOF\nTESTFIRSTFILE\nsetting symmtery-based\nEOF"""%(fname, refname)
	try:
		fh = open("index_check", 'w')
		fh.write("#!/bin/bash\n\n")
		fh.write(input_cmd)
		fh.close()
		chmod_cmd = "chmod +x index_check"; ind_cmd = "./index_check > pointless.log"
		try:
			run_command('sxdm', dirname, user, chmod_cmd, 'merge.log')
			run_command('sxdm', dirname, user, ind_cmd, 'merge.log')
		except (OSError, TypeError, Exception) as e:
			sub.call(chmod_cmd, shell=True)
			sub.call(ind_cmd, shell=True)
			success = True
	except Exception as err:
		logger.error(err)
		success = False
	return success

def read_output():
	results = {"Status": False}
	if os.path.isfile("pointless.log"):
		fh = open("pointless.log", 'r')
		_all = fh.readlines()
		fh.close()
		fkeys = ["Best Solution","Reindex operator", "Laue group probability", "Space group confidence"]

		for ii in range(len(_all)):
			for k in fkeys:
				if k in _all[ii]:
					line = _all[ii].split(":")
					val = line[-1].strip('\n')
					results[k] = val.strip(' ')

					results["Status"] = True
				elif "ERROR" in _all[ii]:
					results['Status'] = False
					return results
				else:
					pass
		return results
	else:
		err = IOError("pointless did not run")
		logger.error(err)
		results["Status"] = False
		return results


def is_correct(xdsfile, reference, dirname=None, user=None):

	success = pointless(xdsfile, reference, dirname, user)

	if not success:
		return False
	else:
		out = read_output()
		logger.info('index_info:{}'.format(out))

	if out["Status"]:
		if out["Reindex operator"] == '[h,k,l]':
			out["indexing"] = "correct"

			return True
		elif out["Reindex operator"] != '[h,k,l]':
			out["indexing"] = "incorrect"
			logger.info('operator:{}'.format(out))

			return False
		else:
			return True
	elif not out['Status']:

		return False

def similar_symmetry(hkl1, hklref):
	xdsdict = {'xds_ascii': hklref}
	refdata = ASCII(xdsdict)
	refdata.get_data(xdsdict)
	xdsdict['xds_ascii'] = hkl1
	tstdata = ASCII(xdsdict)
	tstdata.get_data(xdsdict)
	try:
		return refdata.results['symm'].is_similar_symmetry(tstdata.results['symm'], relative_length_tolerance=0.45, absolute_angle_tolerance=10.5)
	except Exception:
		return False

if __name__ == '__main__':
	print(similar_symmetry(sys.argv[1], sys.argv[2]))

