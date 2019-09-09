'''
Created by S. Basu
24-Jan-2018
'''
import os, sys
import logging

logger = logging.getLogger('Scale&Merge')

def parse_xds_stats(filename):
    result = []; status = False;
    count = 1;
    skip_string = ['RESOLUTION', 'LIMIT']
    search_string = ' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n'
    if not filename.endswith('.LP'):
        status = False
        err = 'Wrong file type %s, it has to be .LP file' %filename
        logger.info('TypeError: {}'.format(err))
        return status, result
    elif not os.path.isfile(filename):
        status = False
        err = 'provided %s file does not exist in the path' %filename
        logger.info('OSError: {}'.format(err))
        return status, result
    try:
        fh = open(filename, 'r')
        _all = fh.readlines()
        fh.close()
        start_idx = [idx for idx, search in enumerate(_all) if search == search_string]
        useful = _all[start_idx[-1]+1:start_idx[-1]+8]; # get only last chunk of statistics tables
        status = True
    except (ValueError, IndexError, TypeError) as err:
        logger.info('Error: {}'.format(err))
        status = False
        return status, result

    for lines in useful:
        if any(k in lines for k in skip_string):
            pass
        else:
            line = lines.split()
            if len(line) > 0:
                each_row={}
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
                result.append(each_row)
                count += 1
            else:
                pass
            status = True

    return status, result

def mean_rmeas_calc(result):
    rmeas_lst = [];
    if len(result) > 0:
        for each_row in result:
            rmeas_lst.append(each_row['rmeas']/100.0)
        mean = (sum(rmeas_lst)/len(rmeas_lst))*100.0
        return mean
    else:
        err = "CORRECT.LP file doesn't have stats"
        logger.info('ValueError: {}'.format(err))

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
	format='%(asctime)s %(levelname)-8s %(message)s',
	datefmt='%a, %d %b %Y %H:%M:%S',
	filename='stats.log',
	filemode='w')
    status, result = parse_xds_stats(sys.argv[1])
    if status == True and len(result) > 0:
        print mean_rmeas_calc(result)
