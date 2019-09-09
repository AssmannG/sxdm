'''
Created on 9-Aug-2016
by basu_s
'''
import os, sys
import subprocess as sub
import shutil
import multiprocessing as mp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--hklfile", type=str,
                        help='provide me a scaled hkl file from xscale')
parser.add_argument("--symm", type=str,
                        help='provide the symmetry or space group')
parser.add_argument("--sites", type=str, nargs='+',
                        help='provide me number of sites.. 1 2 3 so on ')
parser.add_argument("--resolution", type=str, nargs='+',
                        help='provide me the resolution parameters.. 2.2 2.4 3.0 so on')
parser.add_argument("--emins", type=str, nargs='+',
                        help='provide Emins parameters.. 1.2 1.3 1.4 so on')
parser.add_argument("--ntries", type=str, nargs='+',
                        help='provide number of trials.. 2000 3000 so on')

args = parser.parse_args()

def shelxc_infile(ref_file, outfile, cell, symm, sites):
    fh = open(outfile, 'w')
    fh.write('SAD ' + ref_file + '\n')
    fh.write('CELL ' + cell + '\n')
    fh.write('SPAG ' + symm + '\n')
    fh.write('FIND ' + str(sites) + '\n')
    fh.write('SFAC S\n')
    fh.write('MAXM 1000\n')
    fh.close()

def run_shelxc(project, inp):
    sub.call(["shelxc "+ project + "< "+inp + "| tee " + project + "-shelxc.log"], shell=True)


def create_copies(ifh, ofh):
    ifh = ifh + '_fa.hkl'
    ofh = ofh + '_fa.hkl'
    shutil.copyfile(ifh, ofh)
#    shutil.move(ifh, ofh)


def shelxd_input(outname, inname, resolution, iters, emins):
    '''function that modifies .ins file for running shelxd.
    '''

    #os.rename(filename, filename+'~')
    new = open(outname, 'w')
    select = [];   keys = ['SHEL', 'NTRY', 'MIND', 'NTPR']
    source = open(inname, 'r')
    all_lines = source.readlines()
    source.close()
   # all_lines.pop(11); all_lines.pop(13); #remove SHEL and NTRY lines..buggy. case when multiple SYM lines appear.
#    all_lines.pop(7); all_lines.pop(9); #remove SHEL and NTRY lines..

    for line in all_lines:
       # line = lines.split()
        if not any(string in line for string in keys):
           select.append(line)

    all_lines = select

    for ii in xrange(len(all_lines)):
        line = all_lines[ii]
        new.write(line)
        line = line.split()

        if 'FIND' in line:

            new.write('SHEL 999 ' + str(resolution)+'\n')
            new.write('ESEL ' + str(emins)+ '\n')
            new.write('MIND -3.5 2.8 \n')
            new.write('NTPR 600\n')
        if 'PATS' in line:

            new.write('NTRY '+ str(iters) + '\n')
    new.close()

    #sub.call(["rm "+filename+"~"], shell=True)


def create_ins(base, cell, symm, site, reso, iters):
    name = base + '.ins'
    new = open(name, 'w')
    new.write('TITL '+ name +' SAD in '+ symm + '\n')
    new.write('CELL 0.98000 ' + cell +'\n')
    new.write('LATT -1' + '\n')
    new.write('SYMM -X, 1/2+Y, -Z' + '\n')
    new.write('SFAC S' + '\n')
    new.write('UNIT 64' + '\n')
    new.write('FIND ' + str(site) + '\n')
    new.write('SHEL 999 ' + str(reso) + '\n')
    new.write('MIND -1.5 -0.1' + '\n')
    new.write('PATS' + '\n')
    new.write('NTRY ' + str(iters) + '\n')
    new.write('SEED 1' + '\n')
    new.write('HKLF 3' + '\n')
    new.write('END' + '\n')
    new.close()


def run_shelxd(filename, outname):
    sub.call(["shelxd " + filename + " | tee " + outname], shell=True)

def get_20percent(total_sites):
    val = total_sites*0.2
    return round(val)

def reso_range(start, stop, step):
    i = start
    while i <= stop:
          yield i
          i += step


def get_string(filename, string):
    file = open(filename, 'r')
    all_lines = file.readlines()
    found = False
    select = []
    for lines in all_lines:
        line = lines.split()
        if string in line:
            found = True
            select.append(line)
    return select


def par_shelxc(ref_file, cell, symm, tot_sites):
    proc = [];
    for sites in tot_sites:
        project = 'tr-'+str(sites)
        inp = str(sites)+'.inp'
        shelxc_infile(ref_file, inp, cell, symm, sites)
        run_shelxc(project,inp)


def par_shelxd(tot_sites, reso, emins, ntry):

    proc = []; #pool = mp.Pool(processes=4)
    for sites in tot_sites:
        sites = int(sites)
        for res in reso:
            res = float(res)
            for emin in emins:
                emin = float(emin)
                for trial in ntry:
                    trial = int(trial)
                    base_c = 'tr-'+str(sites)
                    project = base_c +'-'+str(res)+'-'+str(emin)+'-'+str(trial)
                    create_copies(base_c, project)
                    #create_ins((project+'_fa'), cell, symm, sites, res, trial)
                    inname = base_c + '_fa' + '.ins'
                    name = project + '_fa'
                    outname = name + '.ins'
                    logname = project +'-shelxd.log'
                    shelxd_input(outname, inname, res, trial, emin)
                    run_shelxd(name, logname)


def main():
    if args.hklfile is None:
       sys.exit('Need a hkl/ahkl file. --help\n')
    else:
       ref_file = args.hklfile
       cell = sub.check_output(["grep CELL_CONSTANT "+ref_file + " | awk '{print $2,$3,$4,$5,$6,$7}'" ], shell=True)
    '''
    if len(args.cell) != 6:
       sys.exit('Need a,b,c and all angle.. --help\n')
    else:
       cell = args.cell
       cell = " ".join(cell)
    '''
    if args.sites is None:
       sys.exit('Need at least one sites.. --help\n')
    else:
       tot_sites = args.sites

    if args.symm is None:
       sys.exit('Need a space group.. --help\n')
    else:
       symm = args.symm

    if args.resolution is None:
       sys.exit('Need atleast one resolution parameter.. --help\n')
    else:
       reso = args.resolution
    if args.emins is None:
       sys.exit('Need atleast one Emin parameter.. --help\n')
    else:
       emins = args.emins

    if args.ntries is None:
       sys.exit('Need atleast one ntries parameter.. --help\n')
    else:
       ntry = args.ntries
    #cell = "135.29    300.64    145.10  90.000 113.343  90.000"
    par_shelxc(ref_file, cell, symm, tot_sites)
    par_shelxd(tot_sites, reso, emins, ntry)

if __name__ == '__main__':
   main()
