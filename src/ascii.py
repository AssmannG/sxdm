__author__ = ["S. Basu"]
__date__ = "Sept-2017"
__refactordate__ =  "11/05/2021"

import time
import errno
import re
from src.run_command import *
import subprocess as sub
from src.abstract import Abstract


logger = logging.getLogger("sxdm")
try:
    from cctbx.array_family import flex
    from cctbx import crystal
    from cctbx import miller
except (ImportError, RuntimeError) as err:
    logger.error(err)
class ASCII(Abstract):

    @property
    def getInDataScheme(self):
        return {
            "type": "object",
            "required": ['xds_ascii'],
            "properties": {
                'xds_ascii': {"type": "string"},
                "resolution": {"type": "number"},
                "mtz_format": {"type": "string"},
                "user": {"type": "string"}
            }
        }

    @property
    def getOutDataScheme(self):
        return {
            "type": "object",
            "properties": {
                "anom": {"type": "string"},
                "spg": {"type": "string"},
                "wave": {"type": "string"},
                "unit_cell": {"type": "string"},
                "symm": {"type": "string"},
                "indices": {
                    "type": "array",
                    "items": {"type": "integer"}
                },
                "iobs": {
                    "type": "array",
                    "items": {"type": "number"}
                },
                "sigI": {
                    "type": "array",
                    "items": {"type": "number"}
                },
                "input_files": {
                    "type": "array",
                    "items": {
                        "type": "array",
                        "items": [
                            {"type": "string"},
                            {"type": "number"}
                        ]
                    }
                },
                "header": {
                    "type": "object",
                    "patternProperties": {
                        "^[a-zA-Z]+$": {"type": "string"}
                    }
                },
                "mulitiplicity": {
                    "type": "array",
                    "items": {"types": "number"}
                }
            }
        }

    def get_data(self, inData):
        if os.path.isfile(inData['xds_ascii']):
            self.readfile(inData)
            self._extract()
        else:
            err = "XDS_ASCII file cannot be read\n"
            logger.error(err)
        return


    def readfile(self, inData):
        self.results['header'] = dict()
        self.results['input_files'] = dict()
        self.results['unit_cell'] = dict()
        self.results['indices'] = []
        self.results['iobs'] = []
        self.results['sigI'] = []

        fname = inData['xds_ascii']
        fh = open(fname, 'r')
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


                if "INPUT_FILE" in params:
                    self.results['input_files'][iset] = params["INPUT_FILE"]
                # elif "X-RAY_WAVELENGTH" in params:
                #     self.results['input_files'][iset][1] = params["X-RAY_WAVELENGTH"]
                else:
                    pass


            elif "!" in lines:
                lst = kwd.findall(lines)
                try:
                    for item in lst:
                        self.results['header'][item[0]] = item[1]
                except IndexError:
                    pass
            else:
                pass

        for field in chunk:
            line = field.split()
            index = [int(line[0]), int(line[1]), int(line[2])]
            self.results['indices'].append(index)
            self.results['iobs'].append(float(line[3]))
            self.results['sigI'].append(float(line[4]))

        self.results['indices'] = flex.miller_index(self.results['indices'])
        self.results['iobs'], self.results['sigI'] = map(lambda x:flex.double(x), (self.results['iobs'], self.results['sigI']))
        return

    def _extract(self):
        try:
            self.results['anom'] = self.results['header']["FRIEDEL'S_LAW"].strip(' ')
            self.results['spg'] = self.results['header']["!SPACE_GROUP_NUMBER"]
            cell = self.results['header']["!UNIT_CELL_CONSTANTS"]
            cell = cell.split()
            self.results['unit_cell']['a'] = cell[0]
            self.results['unit_cell']['b'] = cell[1]
            self.results['unit_cell']['c'] = cell[2]
            self.results['unit_cell']['al'] = cell[3]
            self.results['unit_cell']['be'] = cell[4]
            self.results['unit_cell']['ga'] = cell[5]
            #cctbx crystal symmetry into flex array type
            a, b, c, al, be, ga = map(lambda x:float(x), cell)
            self.results['symm'] = crystal.symmetry(unit_cell=(a, b, c, al, be, ga), space_group=int(self.results['spg']))
            self.results['wave'] = self.results['header']['!X-RAY_WAVELENGTH']
        except (KeyError, IndexError) as e:
            # err = "ASCII-class: could be XSCALE file or check file for keyErrors\n"
            pass
            # logger.info(e)
        return

    def create_miller_set(self):
        if self.results['anom'] == "TRUE":
            anom_flag = True
        elif self.results['anom'] == "FALSE":
            anom_flag = False
        else:
            anom_flag = None

        return miller.set(crystal_symmetry=self.results['symm'], indices=self.results['indices'], anomalous_flag=anom_flag)

    def i_obs(self):
        array_info = miller.array_info(source_type="xds_ascii")
        iobs = miller.array(self.create_miller_set(),
                            data=self.results['iobs'], sigmas=self.results['sigI']).set_info(array_info).set_observation_type_xray_intensity()

        return iobs

    def multiplicity_check(self):
        data = self.i_obs()
        merge = data.merge_equivalents(use_internal_variance=False)
        mult_list = []
        redundancy = merge.redundancies().data()
        try:
            for x in sorted(set(redundancy)):
                each_data_point = [x, redundancy.count(x)]
                mult_list.append(each_data_point)
            #import matplotlib.pyplot as plt
            #plt.hist(self.mult_dict.values())
            #plt.savefig('multiplicity-distribution.png')

        except Exception as err:
            logger.info('multiplicity_check_error:{}'.format(err))
        self.results['multiplicity'] = mult_list
        return

    def xdsconv(self, inData):
        fname = inData['xds_ascii']
        linkname = "data.hkl"
        res = inData.get('resolution', 0.0)
        form = inData['mtz_format']
        if os.path.isfile(fname):
            try:
                os.symlink(fname, linkname)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    os.remove(linkname)
                    os.symlink(fname, linkname)
                else:
                    logger.info('xdsconv-error:{}'.format(e))
        else:
            logger.info('xdsconv_error:{}'.format('file does not exist'))

        namelist = os.path.basename(fname).split('.')
        rootname = namelist[0]
        mtzname = str(rootname + "_" + form + ".mtz")

        xdsconv_cmd = 'xdsconv > /dev/null'
        f2mtz_cmd = 'f2mtz HKLOUT %s < F2MTZ.INP > /dev/null' %mtzname

        fh = open('XDSCONV.INP', 'w')
        fh.write("OUTPUT_FILE=temp.hkl  %s\n" %form)
        fh.write("INPUT_FILE=%s\n" %linkname)
        fh.write("FRIEDEL'S_LAW=%s\n" %self.results['anom'])
        fh.write("INCLUDE_RESOLUTION_RANGE=50 %f\n" %res)
        fh.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS= 0.05\n")
        fh.close()
        logger.info("running xdsconv")

        try:
            run_command("sxdm", self.getOutputDirectory(), inData['user'], xdsconv_cmd, 'xdsconv.log')
        except KeyError:
            sub.call(xdsconv_cmd, shell=True)
        time.sleep(2)
        if os.path.isfile('F2MTZ.INP'):
            try:
                run_command("sxdm", self.getOutputDirectory(), inData['user'], f2mtz_cmd, 'xdsconv.log')
            except KeyError:
                sub.call(f2mtz_cmd, shell=True)
            time.sleep(2)
            os.remove('temp.hkl')
        else:
            logger.info('f2mtz_error:{}'.format('xdsconv did not run\n'))
        return

if __name__ == '__main__':
    indict = {"xds_ascii": "/nfs/ssx/shbasu/MEmmery/processed/adm_serial-xtal/adm_45/ISa_Select.HKL"}
    xscale = ASCII(indict)
    xscale.get_data(indict)
    print(xscale.results)