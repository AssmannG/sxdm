import os, sys
import pathlib
import logging
from src.Merging import Merging
logger = logging.getLogger('sxdm')

def finder(root, expt):
    #print(len(root))
    path_list =[]
    if expt == 'serial-xtal':
        for ii in range(len(root)):
            dirs = pathlib.Path(root[ii])
            posix = list(dirs.glob('*/*'))
            path_list += list(map(lambda x: str(x), posix))

    elif expt == 'native-sad':
        for ii in range(len(root)):
            dirs = pathlib.Path(root[ii])
            paths = dirs.glob('set*')
            path_list.append(paths)
    return path_list

def get_paths_xscale():
    xscale_file = sys.argv[1]
    paths = []
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

def optargs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=str, nargs='+')
    parser.add_argument("--expt", type=str)
    parser.add_argument("--isa_cut", type=str, default='3.0')
    parser.add_argument("--reference", type=str)
    parser.add_argument("--res_cut", type=str, default='2.0')
    parser.add_argument("--friedel", type=str, default="FALSE")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    datefmt='%y-%m-%d %H:%M',
    filename='merge.log',
    filemode='a')

    op = optargs()

    if op.root is not None and op.expt is not None:
        hklpath_list = finder(op.root, op.expt)

    else:
        logger.error("command line arguments are missing, quit!")
        sys.exit()

    inData=dict()
    # orginally it was pathlist, now xtallist or dirlist (see email shibom
    inData['xtallist'] = hklpath_list
    inData['experiment'] = op.expt
    inData['reference'] = op.reference
    inData['resolution'] = op.res_cut
    inData['friedels_law'] = op.friedel
    inData['running_folder'] = None

    
    mm = Merging(inData)

    mm.run_()
    if mm.is_success():
        mm.writeOutputData(mm.results)
