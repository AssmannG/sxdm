import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import glob
from matplotlib.backends.backend_pdf import PdfPages


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

def list_file(dirname):

    file_list = [];
    os.chdir(dirname)
    for file in glob.glob("tr-*shelxd.log"):
        file_list.append(file)
    file_list.sort()

    return file_list

def get_CCs(filename):

    CCs = get_string(filename, 'All/Weak')
    CCall = []; CCweak = [];

    # a nasty way of sanity checks on CCs list if it got some funny string pattern

    for line in CCs:
        sub = 'CFOM-'
        for s in line:
            if sub in s:
               CCs.remove(line)
    # hope now it should be clean..
    for line in CCs:
        if len(line) == 15:
           try:
             CCall.append(float(line[6].strip(',')))
             CCweak.append(float(line[8].strip(',')))
           except ValueError:
               pass
        else:
           try:
             CCall.append(float(line[5].strip(',')))
             CCweak.append(float(line[7].strip(',')))
           except ValueError:
               pass

    return CCall, CCweak

class shelxd_plotter(object):
    def __init__(self, folder):
        self.directory = folder
        self.frame_list = [];
       # self.create_list_frame()
       # self.plot_cc()
    def create_list_frame(self):
        lists = list_file(self.directory)
        nums = len(lists)
        frame = int(nums/25) + 1;
        for jj in range(frame):
            start = 25*jj; stop = 25*(jj+1)
            try:
               file_frame = lists[start:stop]
               self.frame_list.append(file_frame)
            except IndexError:
                file_frame = lists[start:(nums-start)]
                self.frame_list.append(file_frame)


    def plot_cc(self):
        with PdfPages("CCplot.pdf") as pdf:
            for each_frame in self.frame_list:
                fig, axs = plt.subplots(5,5, sharex=True, sharey=True, squeeze=True, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.5, wspace=0.1)

                fig.text(0.5, 0.04, 'CCweak', ha='center')
                fig.text(0.04, 0.5, 'CCall', va='center', rotation='vertical')
                fig.set_size_inches(8,11)

                axs = axs.flatten()

                for i in range(len(each_frame)):
                   ccall, ccweak = get_CCs(each_frame[i])
                   if len(ccweak) > len(ccall):
                     diff = len(ccweak) - len(ccall)
                     axs[i].plot(ccweak[0:(len(ccweak)-diff)], ccall, 'o', linewidth=1, rasterized=True)
                   elif len(ccweak) < len(ccall):
                      diff = len(ccall) - len(ccweak)
                      axs[i].plot(ccweak, ccall[0:(len(ccall)-diff)], 'o', linewidth=1, rasterized=True)
                   else:
                      axs[i].plot(ccweak, ccall, 'o', linewidth=1, rasterized=True)

                   titles = each_frame[i].strip('-shelxd.log')
                   axs[i].set_title(titles, fontsize=8)
                pdf.savefig(fig)
                plt.close()

def main():
    w = shelxd_plotter(sys.argv[1])
    w.create_list_frame()
    w.plot_cc()

if __name__ == '__main__':
    main()
            
