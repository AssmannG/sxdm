'''
Created on 13-Oct-2016
by basu_s
'''
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import PyQt4.QtGui as Qt
import glob
import subprocess as sub
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

def plot_shelxc_multi(file_list):
    in_list = open(file_list)
    for line in in_list:
        line = line.strip('\n')
        xc = utils(line, None)
        stats = xc.get_shelxc_stats()

        plt.plot(stats[0,:], stats[6,:], '-o', linewidth=2)
        plt.xlabel('Resolution',fontsize=12, fontweight='bold')
        plt.ylabel('d"/sig',fontsize=12, fontweight='bold')
    plt.show()

class utils():
    def __init__(self, inputs=None, folder=None):
        self.filename = inputs;
        self.string = None
        self.dirname = folder

    def extract_xscale(self):
        try:
          fh = open(self.filename, 'r')
          line_num = sub.check_output(["grep -n 'SUBSET' XSCALE.LP | awk -F: '{print $1}'"], shell=True)
          line_end = int(line_num) + 22;
          #print 'awk "{NR>= %s && NR<= %d {print}}" XSCALE.LP' %(line_num, line_end)
          #sub.call(['awk "NR>= %s && NR<= %d {print}" XSCALE.LP' %(line_num, line_end)], shell=True)

          all_lines = fh.readlines()
          fh.close()
        except OSError:
          print "file not found"
        #look for indexes of the specific lines

        start = all_lines.index(' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n')
        end = all_lines.index(' ========== STATISTICS OF INPUT DATA SET ==========\n')
        # read only the selected chunk ..
        select = all_lines[start:end]
        # prune a bit of junks and then store everything else in a list
        useful = [];
        for lines in select:
           if 'SUBSET' in lines:
              pass
           elif 'RESOLUTION' in lines or 'LIMIT' in lines:
               pass
           elif 'total' in lines:
               pass
           else:
              line = lines.split()
              if len(line) > 0:
                 useful.append(line)
        return useful


    def get_stats(self, data_list, min_range, max_range):
        reso = []; Isig = []; Ano = []; Rfactor = [];
        norm = min_range + 2
        for ii in xrange(min_range, max_range):
          if len(data_list) > 0:

            reso.append(float(data_list[ii][0]))
            Isig.append(float(data_list[ii][8]))
            Ano.append(float(data_list[ii][12]))

            Rfactor.append(float(data_list[ii][5].strip('%')))

        return reso, Isig, Ano, Rfactor

    def get_string(self):
        file = open(self.filename, 'r')
        all_lines = file.readlines()
        found = False
        select = []
        for lines in all_lines:
            line = lines.split()
            if self.string in line:
                found = True
                select.append(line)
        return select

    def list_file(self):
        file_list = [];
        os.chdir(self.dirname)
        for file in glob.glob("tr-*shelxd.log"):
            file_list.append(file)
        file_list.sort()

        return file_list

    def get_contrast(self):
        self.string = "Contrast"
        contrast = self.get_string()
        cont_lst = []
        for letter in contrast:
            tmp = letter[5].strip(',')
            cont_lst.append(tmp)
        return cont_lst

    def get_CCs(self):
        self.string = 'All/Weak'
        CCs = self.get_string()
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

    def get_sites(self):
#        self.string = '1'
#        sites = self.get_string()
        self.string = 'CFOM'
        cfom = self.get_string()
        for vals in cfom:
            CCall = vals[5].strip(',') + '%'
            CCweak = vals[7].strip(',')+ '%'
            fom = vals[9].strip(',')+ '%'

        occu = []; peaks = [];
        count = 0;
        self.string = '1'; sites = self.get_string()
        for line in sites:

            tmp = line[5].strip(',')
            occu.append(tmp)
            count += 1
            peaks.append(count)

        return occu, peaks, CCall, CCweak, fom

    def get_shelxc_stats(self):
        metrics = ['Resl.','Chi-sq', '<I/sig>', '%Complete', 'Multipl.', 'R(pim)%', '<d"/sig>', 'CC(1/2)']
        stats = []; tmp = []; reso = [];
        for ele in metrics:
            self.string = ele
            line = self.get_string()
            stats.append(line)
        for vals in stats:
            for val in vals:
                tmp.append(val)

        tmp[0].pop(0); tmp[0].pop(0);
        tmp[1].pop(0); tmp[2].pop(0);
        tmp[3].pop(0); tmp[4].pop(0);
        tmp[5].pop(0); tmp[6].pop(0);
        tmp[8].pop(0); tmp.pop(7)

        stat_array = np.array(tmp, dtype=float)

        return stat_array

    def create_list_frame(self):
        lists = self.list_file()
        nums = len(lists)
        frame = int(nums/25) + 1;
        frame_list = [];
        for jj in range(frame):
            start = 25*jj; stop = 25*(jj+1)
            try:
               file_frame = lists[start:stop]
               frame_list.append(file_frame)
            except IndexError:
                file_frame = lists[start:(nums-start)]
                frame_list.append(file_frame)
        return frame_list

class Canvas(Qt.QWidget):
    def __init__(self, stringX, stringY):
        Qt.QWidget.__init__(self)

        self.stringX = stringX
        self.stringY = stringY

        self.layout = Qt.QVBoxLayout()
    def other_xaxis(self, axis):

       axs2 = axis.twiny()
       axs_ticks = axis.get_xticks(); ds = [];
       for val in axs_ticks:
          ds.append("%3.1f" %(1/val))
       axs2.set_xticks(axs_ticks)
       axs2.set_xbound(axis.get_xbound())
       axs2.set_xticklabels(ds)
       axs2.set_xlabel("Resolution ($\mathbf{\AA}$)")

    def plot_xscale(self, reso, Isig, Anom, Rfactor):
        fig, axs = plt.subplots(3,1)
        axs = axs.flatten()
        dstar = [];
        for val in reso:
            dstar.append((1.0/val))
        axs[0].plot(dstar, Isig, '-o', linewidth=2)
        axs[0].set_xlim(max(reso), min(reso))
        axs[0].set_xlabel('1/d ($\mathbf{\AA^{-1}}$)')
        axs[0].set_ylabel("I/sig(I)")

        axs[1].plot(dstar, Anom, '-o', linewidth=2)
        axs[1].set_xlim(max(reso), min(reso))
        axs[1].set_xlabel('1/d ($\mathbf{\AA^{-1}}$)')
        axs[1].set_ylabel("Anomalous_signal")

        axs[2].plot(reso, Rfactor, '-o', linewidth=2)
        axs[2].set_xlim(max(reso), min(reso))
        axs[2].set_xlabel('1/d ($\mathbf{\AA^{-1}}$)')
        axs[2].set_ylabel("Rfactor(%)")

        for ii in range(len(axs)):
            other_xaxis(axs[ii])

        self.layout.addWidget(FigureCanvas(fig))
        self.image = FigureCanvas(fig)


    def plot_shelxd(self, data1, data2, opts):
        fig = plt.Figure()
        self.axes = fig.add_subplot(1,1,1)
        if opts == 'site_occupancy':
            self.axes.plot(data1, data2, '-o', linewidth=2)
        elif opts == 'CCall_CCweak':
            if len(data1) > len(data2):
               diff = len(data1) - len(data2)
               self.axes.plot(data1[0:(len(data1)-diff)], data2, 'o', linewidth=2)
            elif len(data1) < len(data2):
               diff = len(data2) - len(data1)
               self.axes.plot(data1, data2[0:(len(data2)-diff)], 'o', linewidth=2)
            else:
               self.axes.plot(data1, data2, 'o', linewidth=2)
        self.axes.set_xlabel(self.stringX, fontsize=12, fontweight='bold')
        self.axes.set_ylabel(self.stringY, fontsize=12, fontweight='bold')
        self.layout.addWidget(FigureCanvas(fig))
        self.image = FigureCanvas(fig)

    def plot_shelxc(self, data1, data2):
        fig = plt.Figure()
        self.axes = fig.add_subplot(1,1,1)
        self.axes.plot(data1, data2, '-o', linewidth=2)
        plt.gca().invert_xaxis()
        self.axes.set_xlabel(self.stringX, fontsize=12, fontweight='bold')
        self.axes.set_ylabel(self.stringY, fontsize=12, fontweight='bold')
        self.layout.addWidget(FigureCanvas(fig))
        self.image = FigureCanvas(fig)

    def plot_cc(self, folder, rows, cols):
        cluster = utils(None, folder)
        lists = cluster.list_file()
        nums = len(lists)
        fig, axs = plt.subplots(rows,cols, sharex=True, sharey=True, squeeze=True, facecolor='w', edgecolor='k')
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.5, wspace=0.1)

        fig.text(0.5, 0.04, self.stringX, ha='center', fontsize=8)
        fig.text(0.04, 0.5, self.stringY, va='center', fontsize=8, rotation='vertical')

        axs = axs.flatten()
        print axs.shape

        for i in range(nums):
            cluster.filename = lists[i]
            ccall, ccweak = cluster.get_CCs()
            if len(ccweak) > len(ccall):
               diff = len(ccweak) - len(ccall)
               axs[i].plot(ccweak[0:(len(ccweak)-diff)], ccall, 'o', linewidth=2)
            elif len(ccweak) < len(ccall):
                diff = len(ccall) - len(ccweak)
                axs[i].plot(ccweak, ccall[0:(len(ccall)-diff)], 'o', linewidth=2)
            else:
                axs[i].plot(ccweak, ccall, 'o', linewidth=2)
            titles = lists[i].strip('-shelxd.log')
            axs[i].set_title(titles, fontsize=8)
        self.layout.addWidget(FigureCanvas(fig))
        self.image = FigureCanvas(fig)

    def plot_pdf(self, folder):
        cc_pdf = utils(None, folder)
        self.frame_list = cc_pdf.create_list_frame()

        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages("CCplot.pdf") as pdf:
            for each_frame in self.frame_list:
                fig, axs = plt.subplots(5,5, sharex=True, sharey=True, squeeze=True, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.5, wspace=0.1)

                fig.text(0.5, 0.04, self.stringX, ha='center')
                fig.text(0.04, 0.5, self.stringY, va='center', rotation='vertical')
                fig.set_size_inches(8,11)

                axs = axs.flatten()

                for i in range(len(each_frame)):
                    cc_pdf.filename = each_frame[i]
                    ccall, ccweak = cc_pdf.get_CCs()
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



class Main(Qt.QMainWindow):
    def __init__(self, inputs, feature, *args):
        Qt.QMainWindow.__init__(self)
        self.filename = inputs; self.opts = feature;
        try:
           self.rows = int(args[0]); self.cols = int(args[1]);
        except IndexError:
            self.rows = None; self.cols = None;

        self.main_widget = Qt.QWidget(self)
        saveFile = Qt.QAction("&Save File", self)
        saveFile.setShortcut("Ctrl+S")
        saveFile.triggered.connect(self.file_save)
        self.statusBar()
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("&File")
        fileMenu.addAction(saveFile)

        self.layout = Qt.QVBoxLayout(self.main_widget)
        self.plotter()
        self.main_widget.setFocus()
        self.main_widget.setLayout(self.layout)
        self.setCentralWidget(self.main_widget)

    def file_save(self):
        try:
           name = unicode(Qt.QFileDialog.getSaveFileName(self, 'Save File'))
           self.canvas.image.print_figure(name, dpi=200)
        except AttributeError:
             name1 = unicode(Qt.QFileDialog.getSaveFileName(self, 'Save File'))
             self.canvas1.image.print_figure(name1, dpi=200)
             name2 = unicode(Qt.QFileDialog.getSaveFileName(self, 'Save File'))
             self.canvas2.image.print_figure(name2, dpi=200)

    def plotter(self):
        if self.opts == 'cluster':
           self.canvas = Canvas("CCweak", "CCall")
           self.canvas.plot_pdf(self.filename)
           self.canvas.setLayout(self.canvas.layout)
           self.layout.addWidget(self.canvas)
           print "saving plot as pdf..\n"
           self.canvas.plot_pdf(self.filename)
           print "done.."

        elif self.opts == 'xscale':
          xscale_stat = utils("XSCALE.LP")
          useful = xscale_stat.extract_xscale()
          if len(useful) > 12:
            res, Isig, Ano, Rf = xscale_stat.get_stats(useful,0,19)
          else:
             res, Isig, Ano, Rf = xscale_stat.get_stats(useful,0,10)

          self.canvas = Canvas('', '')
          self.canvas.plot_xscale(res, Isig, Ano, Rf)

          self.canvas.setLayout(self.canvas.layout)
          self.layout.addWidget(self.canvas)

        else:
           cc = utils(self.filename)
           if self.opts == 'CCall_CCweak':
              CCall, CCweak = cc.get_CCs()
              self.canvas = Canvas("CCweak", "CCall")
              self.canvas.plot_shelxd(CCweak, CCall, self.opts)
              self.canvas.setLayout(self.canvas.layout)
              self.layout.addWidget(self.canvas)
           elif self.opts == 'site_occupancy':
               occu, peaks, ccall,ccweak, fom = cc.get_sites()
               self.canvas = Canvas("Peak Number", "Occupancy")
               self.canvas.plot_shelxd(peaks, occu, self.opts)
               self.canvas.setLayout(self.canvas.layout)
               self.layout.addWidget(self.canvas)
           elif self.opts == 'shelxc':
               shelxc_array = cc.get_shelxc_stats()
               self.canvas1 = Canvas('Resolution','d"/sig')
               self.canvas1.plot_shelxc(shelxc_array[0,:], shelxc_array[6,:])
               self.canvas1.setLayout(self.canvas1.layout)

               self.canvas2 = Canvas('Resolution','CC1/2')
               self.canvas2.plot_shelxc(shelxc_array[0,:], shelxc_array[7,:])
               self.canvas2.setLayout(self.canvas2.layout)

               self.layout.addWidget(self.canvas1)
               self.layout.addWidget(self.canvas2)



def main():
    app = Qt.QApplication(sys.argv)
    w = Main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    w.show()
    return app.exec_()
    plot_shelxc_multi(sys.argv[1])

if __name__ == '__main__':
   main()
