'''
Created on 15-June-2017
by basu_s
'''
import os, sys, time
import subprocess as sub
import optparse, errno
from ascii import ASCII
import multiprocessing as mp

def options():
    parser = optparse.OptionParser()
    parser.add_option('-f', "--hkl",  type=str,
                    help="XDS_ASCII or XSCALE file absolutely needed")
    parser.add_option('-r', "--hres", type=str,
                    help="It is better if provided, default is 0.0")
    parser.add_option('-s', "--seqin", type=str,
                    help="Better if provided, else provide residue number")
    parser.add_option('-n', "--nres",  type=int,
                    help="If no sequence, then give residue number/monomer")
    parser.add_option('-H', "--hatom_pdb", type=str,
                    help="substructure can be provided optionally")
    parser.add_option('-N', "--native", type=str,
                    help="optionally, native dataset can be provided")
    parser.add_option("-a", "--hatom", type=str,
                    help="heavy atom type can be provided, default is S atom")
    parser.add_option("-S", "--sites", type=int,
                    help="heavy atom number can be provided, default is to estimate based on seq or # residue")

    (opts, args) = parser.parse_args()
    return opts, args

class Building(object):

    def __init__(self, xds_ascii, **kwargs):
        self.hklfile= xds_ascii
        self.hkldata = ASCII(self.hklfile)
        self.mtzfile=None
        self.seq_file=kwargs.get('seq_file', None)
        self.solvent_content=kwargs.get('solvent_content', '0.50')
        self.hres=kwargs.get('hres', '1.0')
        self.substr=kwargs.get('substr', "None")
        self.native = kwargs.get('native', None)
        self.fmt = ['CCP4_I+F', 'CCP4']
        self.cell = None; self.symm = None;
        self.sites = kwargs.get('sites', 0); self.atomtype = kwargs.get('atomtype', 'S')
        self.nres = kwargs.get('nres', 0);

    def guess_sites(self):
        if self.nres > 0 and self.seq_file is None and self.sites == 0:
           self.sites = int(self.nres/25.0)
        else:
          print "no. of residue should be greater than 0. defaults is 0.\n"
          print "sites: %d" %self.sites

    def get_wvl(self):
        try:
            self.wvl = self.hkldata.input_files[1][1]
        except (KeyError, IndexError, ValueError) as e:
            raise e
        '''
        wvl_cmd = ["grep ' 1 X-RAY_WAVELENGTH' "+self.hklfile+" | awk '{print $5}'"]
        self.wvl = sub.check_output(wvl_cmd, shell=True)
        self.wvl = self.wvl.strip('\n')
        '''
        print "X-ray Wavelength: %s\n" %self.wvl

    def get_friedel_cell(self, fname):
        '''
        _info = sub.check_output(["grep 'FRIEDEL' "+fname], shell=True)
        words = _info.split()
        self.friedel = words[2]
        cmd = ["grep 'CELL_CONSTANT' "+fname+" | awk '{print $2,$3,$4,$5,$6,$7}'"]
        self.cell = sub.check_output(cmd, shell=True); self.cell = self.cell.strip('\n')
        sg_cmd = ["grep 'SPACE_GROUP' "+fname+" | awk '{print $2}'"]
        self.symm = sub.check_output(sg_cmd, shell=True); self.symm = self.symm.strip('\n')
        '''
        self.friedel = self.hkldata.anom
        self.cell = self.hkldata.header['!UNIT_CELL_CONSTANTS']
        self.symm = self.hkldata.spg
        print "Unit-cell:{}".format(self.cell)
        print "Space group: %s" %self.symm

        return self.friedel

    def create_mtz(self, fname, res, format):

        if os.path.isfile(fname):
           try:
             linkname = "data.hkl"
             os.symlink(fname, linkname)
           except OSError, e:
              if e.errno == errno.EEXIST:
                  os.remove(linkname)
                  os.symlink(fname, linkname)
              else:
                 raise e
           fh = open("XDSCONV.INP", 'w')
           fh.write("OUTPUT_FILE=tmp.hkl  %s\n" %format)
           fh.write("INPUT_FILE=%s\n" %linkname)
           fh.write("INCLUDE_RESOLUTION_RANGE= 50 %s\n" %res)
           fh.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS= 0.05\n")
           fh.write("%s\n" %self.friedel)
           fh.close()
        else:
           sys.exit(1)

        print "running xdsconv\n"
        os.system("xdsconv > /dev/null")

        endfmt = ['.hkl', '.ahkl', '.HKL']
        try:
          self.rootname = os.path.basename(fname)
          for efmt in endfmt:
           if self.rootname.endswith(efmt):
             self.rootname=self.rootname.strip(efmt)
             self.mtzname=self.rootname+'-'+format+'.mtz'
        except OSError:
            print 'HKL file may not exist \n'

        if os.path.isfile("F2MTZ.INP"):
           os.system("f2mtz HKLOUT "+self.mtzname+" < F2MTZ.INP > /dev/null")
           if self.friedel == "TRUE":
             self.mtzfile = self.mtzname;
           os.remove("tmp.hkl")

    def cad_mtz(self):
        mtz1 = self.rootname+'-CCP4_I+F.mtz'
        mtz2 = self.rootname+'-CCP4.mtz'
        self.mtzfile = self.rootname+'.mtz'

        cad_cmd = "cad HKLIN1 %s HKLIN2 %s HKLOUT %s" %(mtz1,mtz2,self.mtzfile)

        fh = open("cad.inp", 'w')
        fh.write("LABIN FILE 1 ALL\n")
        fh.write("LABIN FILE 2 E1=DANO E2=SIGDANO\n")
        fh.write("END")
        fh.close()

        os.system(cad_cmd + "< cad.inp > /dev/null")
        try:
          os.remove(mtz1)
          os.remove(mtz2)
          os.remove('cad.inp')
        except OSError:
            pass

    def prep_mtz(self):

        try:
          self.get_friedel_cell(self.hklfile)
          if self.friedel == "TRUE":
             self.create_mtz(self.hklfile, self.hres, "CCP4_I+F")
          else:
             for format in self.fmt:
                self.create_mtz(self.hklfile, self.hres, format)
             self.cad_mtz()
        except OSError:
            print "hkl file %s was not given\n" %self.hklfile


        if os.path.isfile(self.mtzfile):
           #os.system("mtzdmp "+self.mtzfile)
           print "final mtz file with all required labels is created\n"
        else:
           print "something might go wrong in creating mtzfile. check\n"

    def create_natmtz(self):

        if self.native != None and os.path.isfile(self.native):
            self.get_friedel_cell(self.native)
            self.create_mtz(self.native, "1.0", "CCP4_I+F")
            cad_cmd = "cad HKLIN1 %s HKLOUT %s" %(self.mtzname,"nat.mtz")

            fh = open("nat.inp", 'w')
            fh.write("LABIN FILE 1 E1=FP E2=SIGFP\n")
            fh.write("LABOUT E1=FNAT E2=SIGFNAT\n")
            fh.write("SYMM %s\n" %self.symm)
            fh.write("END")
            fh.close()

            os.system(cad_cmd + "< nat.inp > /dev/null")

        else:
            pass

    def residue_counter(self):
        if self.seq_file != None and os.path.isfile(self.seq_file):
            fh = open(self.seq_file, 'r')
            _all = fh.readlines()
            fh.close()
            key = ['>','|']; res = [];
            n_cys = 0; n_met = 0; self.nres = 0;
        else:
            print "sequence file was not provided or did not exist"
            return
        for lines in _all:
            if not any(k in lines for k in key):
                res.append(lines)
        for vals in res:
            val = vals.split()
            self.nres += len(val[0])
            n_cys += val[0].count('C')
            n_met += val[0].count('M')

        self.sites = n_cys + n_met
        print "No. of residues: %d" %self.nres
        print "No. of sites: %d (CYS: %d, MET: %d)" %(self.sites, n_cys, n_met)



    def matt_calc(self):

        self.get_friedel_cell(self.hklfile)
        try:
           self.residue_counter()
        except OSError, TypeError:
            try:
                self.guess_sites()
            except ValueError:
                raise ValueError("Either Nresidue or sequence needed\n")

        fh = open('matt.inp', 'w')
        fh.write('CELL  %s\n' %self.cell)
        fh.write('SYMM  %s\n' %self.symm)
        fh.write('NRES  %d\n' %self.nres)
        fh.close()
        os.system('matthews_coef < matt.inp > matt.log')

        try:
            fh = open('matt.log', 'r')
            _all = fh.readlines()
            fh.close()
        except OSError:
            raise OSError("matthews_coefficient did not run properly\n")
        keys = ["The Matthews Coefficient is", "the solvent % is"]

        for lines in _all:
            if any(k in lines for k in keys):
                print lines,
                line = lines.split(':')
                self.solvent_content = float(line[-1].strip('\n'))/100.0

        print "Solvent content: %.3f" %self.solvent_content
        os.remove('matt.inp')
        os.remove('matt.log')

    def calc_fps(self):

        fh =open("crossec.inp", 'w')
        fh.write("ATOM %s\n" %self.atomtype)
        fh.write("NWAV 1 %s\n" %self.wvl)
        fh.close()

        os.system("crossec < crossec.inp > out.log")
        fh = open('out.log', 'r')
        _all = fh.readlines(); fh.close()
        val = [];
        for lines in _all:
            if lines.startswith(self.atomtype):
                line = lines.split()
                self.fp = line[2]
                self.fpp = line[3]

        os.remove('crossec.inp')
        os.remove('out.log')

    def crank_(self, **kwargs):

        self.calc_fps()

        fh = open("run-crank", 'w')
        fh.write('#!/bin/bash\n\n')
        fh.write('crank=$CCP4/share/ccp4i/crank2/crank2.py\n')
        fh.write('ccp4-python $crank dirout crank_result << EOF\n')
        fh.write('target::SAD\n')
        if not os.path.isfile(self.substr):

            fh.write('faest\n')
            fh.write('substrdet num_trials::2000 shelxd\n')
            fh.write('model typ=substr atomtype=%s fp=%s fpp=%s exp_num_atoms=%s\n' %(self.atomtype, self.fp, self.fpp, self.sites))
            fh.write('phas\n')
            fh.write('handdet\n')
            fh.write('dmfull\n')
            fh.write('comb_phdmmb\n')
            fh.write('ref target::MLHL\n')
        else:
            fh.write('phas\n')
            fh.write('handdet\n')
            fh.write('dmfull\n')
            fh.write('comb_phdmmb\n')
            fh.write('ref target::MLHL\n')
            fh.write('model typ=substr atomtype=%s fp=%s fpp=%s file=%s\n' %(self.atomtype, self.fp, self.fpp, self.substr))

        fh.write('fsigf typ=plus  f=F(+) sigf=SIGF(+)  file=%s\n' %self.mtzfile)
        fh.write('fsigf typ=minus f=F(-) sigf=SIGF(-)\n')
        fh.write('exclude typ=freeR free=FreeRflag\n')


        if not self.native is None and os.path.isfile(self.native):
            fh.write('fsigf typ=average f=FNAT sigf=SIGFNAT  file=nat.mtz\n')

        if not self.seq_file is None and os.path.isfile(self.seq_file):
            fh.write('sequence file=%s\n' %self.seq_file)
        else:
            fh.write('sequence residues_mon=%s solvent_content=%s\n' %(self.nres,self.solvent_content))

        fh.write('EOF\n')

        fh.close()
        os.system('chmod +x run-crank')
        os.system('./run-crank &')

    def sharp_(self, **kwargs):

        if not self.seq_file is None:
            sharp_cmd = "run_autoSHARP.sh -seq %s -ha %s -nsit %d -wvl %s -mtz %s" %(self.seq_file, self.atomtype, self.sites, self.wvl, self.mtzfile)
        elif self.nres != 0 and self.seq_file is None:
            sharp_cmd = "run_autoSHARP.sh -nres %s -ha %s -nsit %d -wvl %s -mtz %s" %(self.nres, self.atomtype, self.sites, self.wvl, self.mtzfile)
        else:
            raise ValueError("Either sequence or residue numbers needed")
        if not self.native is None and os.path.isfile(self.native):
            sharp_cmd = sharp_cmd + " -nat -mtz %s" %(self.mtzname)

        os.system(sharp_cmd +"&")

    def run_firefox(self):
        cwd = os.getcwd()
        path = os.path.join(cwd,"autoSHARP","LISTautoSHARP.html")
        time.sleep(60)
        if os.path.isfile(path):
            os.system('firefox %s' %path)
        else:
            pass


def main():
    (opts, args) = options()

    if opts.hkl is None:
        sys.exit(1)
    keywords = {}
    if opts.hres != None:
        keywords['hres'] = opts.hres
    if opts.seqin != None:
        keywords['seq_file'] = opts.seqin
    elif opts.nres != None and opts.seqin is None:
        keywords['nres'] = opts.nres
    else:
        print "Need at least seq file or # residues to run the program"
        sys.exit(1)

    if opts.hatom_pdb != None:
        keywords['substr'] = opts.hatom_pdb
    if opts.hatom != None:
        keywords['atomtype'] = opts.hatom
    if opts.native != None:
        keywords['native'] = opts.native
    if opts.sites != None:
        keywords['sites'] = opts.sites



    print keywords
    model = Building(opts.hkl, **keywords)
    model.prep_mtz()
    model.matt_calc()
    model.get_wvl()
    if 'native' in keywords.keys():
        model.create_natmtz()

    id1 = mp.Process(target=model.crank_)
    id2 = mp.Process(target=model.sharp_)
    id3 = mp.Process(target=model.run_firefox)
    jobs = [id1, id2, id3]
    for p in jobs:
        p.start()
    for p in jobs:
        p.join()

#model = Building(sys.argv[1], native=sys.argv[2])
#model.create_natmtz()

if __name__ == '__main__':
    main()
