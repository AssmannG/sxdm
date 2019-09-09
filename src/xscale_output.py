'''
Created by S.Basu
On 28-June-2017
'''
import os, sys
import logging

logger = logging.getLogger('Scale&Merge')

def parse_xscale_output(filename):
	results = []; status = False;
	count = 1;
	'''
	results = {
	'Status' : False,
	'resolution_limit': [],
	'observed_reflections': [],
	'unique_reflections' : [],
	'possible_reflections' : [],
	'completeness' : [],
	'rfactor_observed' : [],
	'rfactor_expected' : [],
	'compared_reflections' : [],
	'i_sigma' : [],
	'r_meas' : [],
	'cc_half' : [],
	'cc_ano' : [],
	'sig_ano' : [],
	'anom_compared' : []
	}
	'''
	skip_line = ['RESOLUTION', 'LIMIT', 'total']

	if not filename.endswith(".LP"):
		e = TypeError("wrong file type %s, need .LP file\n" %filename)
		logger.info('TypeError: {}'.format(e))
		status = False
		return status, results
	elif not os.path.isfile(filename):
		e = IOError("%s file does not exist\n" %filename)
		logger.info('Error: {}'.format(e))
		status = False
		return status, results

	fh = open(filename, 'r')
	_all = fh.readlines(); fh.close()
	try:
		start = _all.index(' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n')
		end = _all.index(' ========== STATISTICS OF INPUT DATA SET ==========\n')
		useful = _all[start+1:end-1]; status = True

	except (TypeError, ValueError, IndexError) as err:
		logger.info('Error:there was no XDS_ASCII files, so nothing to XSCALE. Check input/output directory naming')
		logger.info(err); status = False
		return status, results

	for lines in useful:
		if any(k in lines for k in skip_line):
			pass
		else:
			logger.info(lines[:-1])
			line = lines.split()
			if len(line) == 14:
				each_row = {}
				each_row['order'] = count
				each_row['resolution_limit'] = float(line[0])
				each_row['observed_reflections'] = int(line[1])
				each_row['unique_reflections'] = int(line[2])
				each_row['possible_reflections'] = int(line[3])
				each_row['completeness'] = float(line[4].strip('%'))
				each_row['r_factor_observed'] = float(line[5].strip('%'))
				each_row['r_factor_expected'] = float(line[6].strip('%'))
				each_row['compared_reflections'] = int(line[7])
				each_row['Isigma'] = float(line[8])
				each_row['rmeas'] = float(line[9].strip('%'))
				each_row['sigAno'] = float(line[12])
				each_row['anom_compared'] = float(line[13])
				if line[10][-1] == '*':
					each_row['cc_half'] = float(line[10].strip('*'))
				else:
					each_row['cc_half'] = float(line[10])
				if line[11][-1] == '*':
					each_row['anomalous_correlation'] = float(line[11].strip('*'))
				else:
					each_row['anomalous_correlation'] = float(line[11])
				results.append(each_row)
				count += 1; status = True
			elif len(line) > 0 and len(line) < 14:
				msg = "float formatting error in xscale table"
				logger.info('xscale-output-error:{}'.format(msg))
				each_row = {}
				each_row['order'] = count
				each_row['resolution_limit'] = float(line[0])
				each_row['observed_reflections'] = int(line[1])
				each_row['unique_reflections'] = int(line[2])
				each_row['possible_reflections'] = int(line[3])
				each_row['completeness'] = float(line[4].strip('%'))
				each_row['r_factor_observed'] = 'NaN'
				each_row['r_factor_expected'] = 'NaN'
				each_row['compared_reflections'] = 'NaN'
				each_row['Isigma'] = 'NaN'
				each_row['rmeas'] = 'NaN'
				each_row['sigAno'] = 'NaN'
				each_row['anom_compared'] = 'NaN'
				each_row['cc_half'] = 'NaN'
				each_row['anomalous_correlation'] = 'NaN'
				results.append(each_row)
				count += 1; status = True
			else:
				pass

	return status, results

def display_table(results):
	ofh = open('table.txt', 'w')
	ofh.write('Resolution   observed    unique   possible  completeness    Robs      Rexpect   compared  I/sig      Rmeas     CC1/2   CCano  anom_comp  sig_ano\n')
	#ofh.write("%10s %12s %10s %10s %12s %11s %11s %9s %8s %10s %9s %7s %8s %9s\n"\
	#	%('Resolution','observed', 'unique', 'possible', \
	#	'completeness', 'Robs', 'Rexpect', 'compared', 'I/sig', 'Rmeas', 'CC1/2', 'CCano' \
	#	'anom_compared', 'sig_ano'))
	for row in results:
		ofh.write("%6.2f %12d %10d %10d %12.1f %11.1f %11.1f %9d %8.1f %10.1f %9.1f %7.1f %8d %9.3f\n"\
			%(row['resolution_limit'], row['observed_reflections'], \
			row['unique_reflections'], row['possible_reflections'], \
			row['completeness'], row['r_factor_observed'], \
			row['r_factor_expected'], row['compared_reflections'], \
			row['Isigma'], row['rmeas'], row['cc_half'], \
			row['anomalous_correlation'], row['anom_compared'], row['sigAno']))
	ofh.close()

def print_table(results):
	if results:
		print 'Resolution   observed    unique   possible  completeness      Robs      Rexpect   compared  I/sig      Rmeas     CC1/2   CCano  anom_comp  sig_ano'
		for row in results:
			print "%6.2f %12d %10d %10d %12.1f %11.1f %11.1f %9d %8.1f %10.1f %9.1f %7.1f %8d %9.3f"\
				%(row['resolution_limit'], row['observed_reflections'], \
				row['unique_reflections'], row['possible_reflections'], \
				row['completeness'], row['r_factor_observed'], \
				row['r_factor_expected'], row['compared_reflections'], \
				row['Isigma'], row['rmeas'], row['cc_half'], \
				row['anomalous_correlation'], row['anom_compared'], row['sigAno'])
	else:
		print "no stats to display, xscale did not run\n"

if __name__ == '__main__':
	logging.basicConfig(level=logging.DEBUG,
	format='%(asctime)s %(levelname)-8s %(message)s',
	datefmt='%a, %d %b %Y %H:%M:%S',
	filename='stats.log',
	filemode='w')
	status, result = parse_xscale_output(sys.argv[1])
	print_table(result)
