import sys, os

# prevent creation of compiled bytecode files
sys.dont_write_bytecode = True

import HTInspect.hitime.HTS_resutls_file_parser as HTS_FP
import HTInspect.hitime.HT_search_postprocessing as HTS_PP

from HTInspect.res.utils import *
from HTInspect.gui import HT_search

from PyQt5 import QtCore, QtGui
import pyqtgraph as pg
from pyqtgraph.Point import Point

import numpy as np

'''
HT_search
'''
class HT_search (QtGui.QDialog, HT_search.Ui_Dialog):
    '''
    HT_search inherits from HT_search.Ui_Dialog
    '''
    def __init__(self, darkMode = True, parent = None, runner = None):
        '''
        Initialise class and call __init__ from super classes
        '''
        # TODO: Upgrade HT_search UI

        # Set PG plot background to white before GUI class init
        if not darkMode:
            pg.setConfigOption('background','w')
            #pg.setConfigOption('foreground', 'k')

        super(HT_search, self).__init__(parent)

        self.setupUi(self)
        self.setWindowTitle('HTInspect')

        if runner:
            self.runner, self.runner_thread, self.q = runner

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.check_response)

        self.RP_min_intensity.setText(str(0))
        self.RP_mzDelta.setText(str(3.010051111))
        self.RP_EIC_width.setText(str(0.03))

        self.HT_search_list = []
        self.hits = [] # stores hitime hits for RV tab
        self.raw_data = None # raw hitime output data for RV tab
        self.raw_HT_RP_data = None # HT_RP raw data holding list
        self.refined_HT_RP_data = None # RP refined data

        self.RP_HT_file = None
        self.RP_mzML_file = None
        self.HT_RP_outputFileName = None

        ''' set plot config options '''
        # HM
        self.HT_RV_HM = pg.ScatterPlotItem()
        # nb labels/gridlines are applied to the widget canvas
        # rather than the scatterplotitem
        self.HT_RV_HM_widget.showGrid(x = True, y = True)
        self.HT_RV_HM_widget.showLabel('bottom', show = True)
        self.HT_RV_HM_widget.showLabel('top', show = True)
        self.HT_RV_HM_widget.showLabel('left', show = True)
        self.HT_RV_HM_widget.showLabel('right', show = True)
        self.HT_RV_HM_widget.setLabel(axis = 'bottom', text = 'm/z')
        self.HT_RV_HM_widget.setLabel(axis = 'left', text = 'Retention Time (s)')
        self.HT_RV_HM_widget.setLabel(axis = 'top', text = '')
        self.HT_RV_HM_widget.setLabel(axis = 'right', text = '')
        self.HT_RV_HM_widget.addItem(self.HT_RV_HM)

        # EIC
        self.HT_RV_EIC.showGrid(x = True, y = True)
        self.HT_RV_EIC.showLabel('bottom', show = True)
        self.HT_RV_EIC.showLabel('left', show = True)
        self.HT_RV_EIC.showLabel('right', show = True)
        self.HT_RV_EIC.setLabel(axis = 'bottom', text = 'Retention Time (s)')
        self.HT_RV_EIC.setLabel(axis = 'left', text = 'Intensity')
        self.HT_RV_EIC.setLabel(axis = 'top', text = '')
        self.HT_RV_EIC.setLabel(axis = 'right', text = '')

        # MS
        self.HT_RV_MS.showGrid(x = True, y = True)
        self.HT_RV_MS.showLabel('bottom', show = True)
        self.HT_RV_MS.showLabel('left', show = True)
        self.HT_RV_MS.showLabel('right', show = True)
        self.HT_RV_MS.setLabel(axis = 'bottom', text = 'm/z')
        self.HT_RV_MS.setLabel(axis = 'left', text = 'Intensity')
        self.HT_RV_MS.setLabel(axis = 'top', text = '')
        self.HT_RV_MS.setLabel(axis = 'right', text = '')

        ''' set tablewidget geometries and selection behaviours '''
        self.HT_RV_hitlist.resizeRowsToContents()
        self.HT_RV_hitlist.resizeColumnsToContents()
        self.HT_RV_hitlist.horizontalHeader().setStretchLastSection(True)

        self.HT_RV_accepted_list.resizeRowsToContents()
        self.HT_RV_accepted_list.resizeColumnsToContents()
        self.HT_RV_accepted_list.horizontalHeader().setStretchLastSection(True)

        hitHeaders = self.HT_RV_hitlist.horizontalHeader()
        hitHeaders.setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        hitHeaders.setResizeMode(1, QtGui.QHeaderView.Stretch)
        hitHeaders.setResizeMode(2, QtGui.QHeaderView.Stretch)
        hitHeaders.setResizeMode(3, QtGui.QHeaderView.ResizeToContents)

        acceptedHeaders = self.HT_RV_accepted_list.horizontalHeader()
        acceptedHeaders.setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
        acceptedHeaders.setResizeMode(1, QtGui.QHeaderView.Stretch)
        acceptedHeaders.setResizeMode(2, QtGui.QHeaderView.Stretch)
        acceptedHeaders.setResizeMode(3, QtGui.QHeaderView.ResizeToContents)

        ''' set selection model for TableWidget items - want to select the entire row when a cell is clicked '''
        self.HT_RV_hitlist.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.HT_RV_accepted_list.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

        ''' set tablewidget view - remove row counter column '''
        self.HT_RV_hitlist.verticalHeader().setVisible(False)
        self.HT_RV_accepted_list.verticalHeader().setVisible(False)

        ''' set infinite lines for MS L/H markers '''
        self.loLine = pg.InfiniteLine(angle = 90, movable = False, pen = 'r')
        self.hiLine = pg.InfiniteLine(angle = 90, movable = False, pen = 'r')

        '''
                Plotting functios and data for HT postprocessing
                '''
        # HM
        self.RP_heat_map.showGrid(x = True, y = True)
        self.RP_heat_map.showLabel('bottom', show = True)
        self.RP_heat_map.showLabel('top', show = True)
        self.RP_heat_map.showLabel('left', show = True)
        self.RP_heat_map.showLabel('right', show = True)
        self.RP_heat_map.setLabel(axis = 'bottom', text = 'm/z')
        self.RP_heat_map.setLabel(axis = 'left', text = 'Retention Time (s)')
        self.RP_heat_map.setLabel(axis = 'top', text = '')
        self.RP_heat_map.setLabel(axis = 'right', text = '')


        self.RP_min_intensity.setText(str(2))
        ''' form GUI connections '''
        self.HT_RP_connections()
        self.HT_RV_connections()

        ''' field activation defaults '''
        self.fileOpenDialogPath = os.path.expanduser('~')

        self.RP_HT_input_field.setText(str('/home/mleeming/Code/HTInspect/HTInspect/data/BSA_FULL.mzML.max'))
        self.RP_output_file_field.setText(str('/home/mleeming/Code/HTInspect/HTInspect/test'))
        self.RP_mzML_input_field.setText(str('/home/mleeming/Code/HTInspect/HTInspect/data/BSA_FULL.mzML'))

        screenShape = QtGui.QDesktopWidget().screenGeometry()
        self.resize(screenShape.width()*0.8, screenShape.height()*0.8)
        self.splitter.setSizes([100,300])

    def check_response(self):
        '''
        On timer tick - check response from NTPS
                - if process returns done, load results into RV GUI tab
        '''
        while not self.q.empty():
            update = self.q.get()
            if update == 'done':
                time.sleep(0.1)
                self.runner.p.terminate()
                self.timer.stop()
                self.runner_thread.exit()
                print('\nResponse returned to guiprocs')
            else:
                self.textEditAppend(update)
                print (update)
        return


    '''
        DATA POSTPROCESSING FUNCTIONS
        '''
    def HT_RP_connections(self):
        self.RP_HT_input_button.clicked.connect(self.get_HT_input_file)
        self.RP_mascot_input_button.clicked.connect(self.select_mzml_file)
        self.RP_min_intensity.textChanged.connect(self.plot_map)
        self.RP_output_file_button.clicked.connect(self.set_output_file)
        self.RP_run_button.clicked.connect(self.run_data_postprocessing)
        self.submit.clicked.connect(self.parse_HT_data)
        return

    def set_output_file(self):
        self.HT_RP_outputFileName, _ = QtGui.QFileDialog.getSaveFileName(self, 'Select Output File Name', self.fileOpenDialogPath)
        self.fileOpenDialogPath = os.path.dirname(str(self.HT_RP_outputFileName))
        self.RP_output_file_field.setText(self.HT_RP_outputFileName)
        return

    def select_mzml_file(self):
        self.RP_mzML_file, _ = QtGui.QFileDialog.getOpenFileName(self, 'Select mzML Data File', self.fileOpenDialogPath)
        self.fileOpenDialogPath = os.path.dirname(str(self.RP_mzML_file))
        self.RP_mzML_input_field.setText(str(self.RP_mzML_file))
        return

    def get_HT_input_file(self):
        # get input file
        self.RP_HT_file, _ = QtGui.QFileDialog.getOpenFileName(
            self, 'Select HiTIME Results File', '/home/mleeming/Code/HTInspect/HTInspect'
        )
        self.RP_HT_input_field.setText(str(self.RP_HT_file))
        if self.RP_HT_file != '':
            self.parse_HT_data()
        return

    def parse_HT_data(self):
        self.RP_HT_file = str(self.RP_HT_input_field.text())

        self.fileOpenDialogPath = os.path.dirname(str(self.RP_HT_file))

        # run file parser
        self.raw_HT_RP_data = HTS_FP.reader(self.RP_HT_file)

        if not self.raw_HT_RP_data: return

        # parse data and plot histogram
        self.plot_map()
        return


    def run_data_postprocessing(self):
        '''
        Pull postprocessing input parameters and run script
        '''

        minIntensity = str(self.RP_min_intensity.text())
        mzDelta = str(self.RP_mzDelta.text())
        eicWidth = str(self.RP_EIC_width.text())

        minIntensity = float(minIntensity) if minIntensity != '' else 0
        mzDelta = float(mzDelta) if mzDelta != '' else 0
        eicWidth = float(eicWidth) if eicWidth != '' else 0.03

        args = {
                'htIn' :str(self.RP_HT_input_field.text()),
                'mzmlFile' :self.RP_mzML_input_field.text(),
                'mzDelta' : mzDelta,
                'eicWidth' : eicWidth,
                'minIntensity' : minIntensity,
                'outFile' : str(self.RP_output_file_field.text()),
                'peakList' : True,
                'plotEICs' : True,
        }

        run_job(self.runner, self.runner_thread, self.q, HTS_PP.guiRun, args)
        self.timer.start(1000)
        return


    def plot_map(self):

        '''
        Score histogram
        --------------------------------------
        '''

        # TODO >> this is a mess > fixme please

        # clear plots
        self.RP_heat_map.clear()

        minIntensity = float(str(self.RP_min_intensity.text()))
        '''
                Hit Heatmap
                Note - increase plotting speed > only plot highest scoring 10,000 points
                '''
        rts = [float(x.rt) for x in self.raw_HT_RP_data if x.amp > minIntensity]
        mzs = [float(x.mz) for x in self.raw_HT_RP_data if x.amp > minIntensity]

        self.RP_heat_map.plot(x = mzs, y = rts, pen = None, symbol = 'o')
        return

    '''
        CONNECTIONS AND FUNCTIONS FOR HITIME RESULTS TAB
        '''
    def HT_RV_connections(self):
        ''' form GUI connections '''
        self.HT_RV_load_results.clicked.connect(self.parse_HT_results_file)
        self.HT_RV_hitlist.itemSelectionChanged.connect(lambda: self.HT_RV_update_plots(None))
        self.HT_RV_accepted_list.itemSelectionChanged.connect(lambda: self.HT_RV_update_plots(None, acceptedHit = True))
        self.HT_RV_up.clicked.connect(lambda: self.HT_RV_update_plots('up'))
        self.HT_RV_down.clicked.connect(lambda: self.HT_RV_update_plots('down'))
        self.HT_RV_reject.clicked.connect(lambda: self.HT_RV_update_plots('down'))
        self.HT_RV_accept.clicked.connect(self.add_hit_to_accepted_list)
        self.HT_RV_accept.clicked.connect(lambda: self.HT_RV_update_plots('down'))
        self.HT_RV_remove.clicked.connect(self.remove_hit_from_accepted_list)
        self.HT_RV_write.clicked.connect(self.write_accepted_to_file)
        self.RV_reset.clicked.connect(self.reset)
        self.HT_RV_HM.sigClicked.connect(self.pointClicked)
        return

    def pointClicked(self, plot, points):

        # find hit list index of clicked point
        for p in points:
            x = str(points[0]._data[0])
            y = str(points[0]._data[1])
            for i, h in enumerate(self.hits):
                if x == str(h.mz) and y == str(h.rt):
                    # relevant result found - highlight row
                    self.HT_RV_hitlist.selectRow(i)
                    # update plots
                    #self.HT_RV_update_plots(None)
                    break
        return

    def keyPressEvent(self, event):
        # update plot ranges when the 'r' key is pressed

        keystroke = str(event.text())
        if keystroke == 'r' or keystroke == 'R':
            self.HT_RV_update_plots(None)
            return
        elif keystroke == '':
            self.close()
        else:
            return

    def parse_HT_results_file(self):
        ''' Called when a file is selected with load results button '''


        inFile, _ = QtGui.QFileDialog.getOpenFileName(self,
                    'Specify Postprocessed HiTIME Reults File',
                    '/home/mleeming/Code/HTInspect/HTInspect'
                    )

        if inFile == '': return

        # safe browser dirname to defaults
        self.fileOpenDialogPath = os.path.dirname(str(inFile))

        print ('Reading results file......')
        # run file parser
        raw_HT_data, HT_hits, self.headers = HTS_FP.processed_results_reader(inFile)

        # self.headers is a dictionary with keys:
        # 'htFile'
        # 'mzmlFile'
        # 'mzWidth'
        # 'rtWidth'
        # 'scoreCutoff';
        # 'mzDelta'
        # 'eicWidth'
        # 'mzML_time_unit'

        # store results in class level variables
        self.raw_data = raw_HT_data
        self.hits = [x for x in HT_hits]
        del raw_HT_data, HT_hits

        self.mzDelta = self.headers['mzDelta']

        # populate resutls table
        for i, hit in enumerate(self.hits):
            self.HT_RV_hitlist.insertRow(i)
            self.HT_RV_hitlist.setItem(i, 0, QtGui.QTableWidgetItem(str(hit.index)))
            self.HT_RV_hitlist.setItem(i, 1, QtGui.QTableWidgetItem(str(hit.rt)))
            self.HT_RV_hitlist.setItem(i, 2, QtGui.QTableWidgetItem('%.4f' %float(hit.mz)))
            self.HT_RV_hitlist.setItem(i, 3, QtGui.QTableWidgetItem('%.1f' %float(hit.score)))

            # centre text in tablewidget columns
            self.HT_RV_hitlist.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
            self.HT_RV_hitlist.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)
            self.HT_RV_hitlist.item(i,2).setTextAlignment(QtCore.Qt.AlignCenter)
            self.HT_RV_hitlist.item(i,3).setTextAlignment(QtCore.Qt.AlignCenter)

        self.HT_HMX = [float(_.mz) for _ in self.hits]
        self.HT_HMY = [float(_.rt) for _ in self.hits]

        # set first row as highlighted
        self.HT_RV_hitlist.selectRow(0) #
        # NOTE: calling selectRow on the tableWidget item also calls the linked updatePlot method
        # ---> separately calling update_plots is unnecessary
        # # run update_plots to plot EIC and MS
        # self.HT_RV_update_plots(None)

        # update table measurements to fit data
        self.HT_RV_hitlist.resizeRowsToContents()
        self.HT_RV_hitlist.resizeColumnsToContents()
        self.HT_RV_hitlist.horizontalHeader().setStretchLastSection(True)
        return

    def HT_RV_update_plots(self, direction, acceptedHit = False):
        '''
        Update data analysis plots and stats when
                1) up or down button is clicked, or
                2) a new peptide entry is clicked in the results list
                3) MS spectrum type is changed using radio buttons

        Parameters:
                direction:
                        'up' if up button clicked
                        'down' if down dubbon clicked
                        'None' if new entry clicked in results list
        '''

        if not acceptedHit:

            # remove any highlighting in accepted list
            self.HT_RV_accepted_list.clearSelection()

            # get selected row from RV_peptides table
            if direction is not None:
                index = self.HT_RV_hitlist.selectionModel().selectedRows()[0].row()
                if direction == 'up':
                    self.HT_RV_hitlist.selectRow(index-1)
                elif direction == 'down':
                    self.HT_RV_hitlist.selectRow(index+1)

            highlighted_row = self.HT_RV_hitlist.selectionModel().selectedRows()[0].row()
            hitNumber = int(str(self.HT_RV_hitlist.item(highlighted_row, 0).text()))

        else:
            try:
                # item selected in accepted list
                highlighted_row = self.HT_RV_accepted_list.selectionModel().selectedRows()[0].row()
                hitNumber = int(str(self.HT_RV_accepted_list.item(highlighted_row, 0).text()))
            except IndexError:
                # land here when a row is removed from accepted list
                hitNumber = 1

        # get matching hit
        hit = None
        for x in self.hits:
            if x.index == hitNumber:
                hit = x
                break

        # clear existing plots
        self.HT_RV_MS.clear()
        self.HT_RV_EIC.clear()

        MS_mz = np.asarray(hit.MS_mz)
        MS_int = np.asarray(hit.MS_int)

        EIC_RT = np.asarray(hit.EIC_RT)
        EIC_int_light = np.asarray(hit.EIC_int_light)
        EIC_int_heavy = np.asarray(hit.EIC_int_heavy)

        # set ranges
        xRangeMS = self.getXRange(hit.mz, 10, rangeType = 'fixed')
        xRangeRT = self.getXRange(hit.rt, 10, rangeType = 'pc')

        try:
            yRangeMS = self.getYRange(MS_mz, MS_int, xRangeMS)
            yRangeRTL = self.getYRange(EIC_RT, EIC_int_light, xRangeRT)
            yRangeRTH = self.getYRange(EIC_RT, EIC_int_heavy, xRangeRT)
            yRangeRT = max([yRangeRTL, yRangeRTH])
        except:
            return

        yRangeMS = ( 0 - yRangeMS * 0.05, yRangeMS)
        yRangeRT = ( 0 - yRangeRT * 0.05, yRangeRT)

        self.HT_RV_EIC.setLimits(xMin = float(min(EIC_RT)), xMax = float(max(EIC_RT)), yMin = yRangeMS[0], yMax = yRangeRT[1]*1.5)
        self.HT_RV_MS.setLimits(xMin = float(min(MS_mz)), xMax = float(max(MS_mz)), yMin = yRangeRT[0], yMax = yRangeMS[1]*1.5)

        self.HT_RV_MS.setRange(xRange = xRangeMS, yRange = yRangeMS)
        self.HT_RV_EIC.setRange(xRange = xRangeRT, yRange = yRangeRT)

        # set iLine positions and add to MS
        self.loLine.setValue(hit.mz)
        self.hiLine.setValue(hit.mz + self.mzDelta)

        self.HT_RV_MS.addItem(self.loLine, ignoreBounds = True)
        self.HT_RV_MS.addItem(self.hiLine, ignoreBounds = True)

        # plot MS
        self.HT_RV_MS.plot(
                            x = np.repeat(MS_mz, 3),
                            y = np.dstack((np.zeros(MS_int.shape[0]), MS_int, np.zeros(MS_int.shape[0]))).flatten(),
                            pen = 'b',
                            connect = 'all'
                            )

        # plot EIC
        self.HT_RV_EIC.plot(x = EIC_RT, y = EIC_int_light, pen = 'b')
        self.HT_RV_EIC.plot(x = EIC_RT, y = EIC_int_heavy, pen = 'r')

        # plot heatmap
        self.HT_RV_HM.clear()
        self.HT_RV_HM.setData(x = self.HT_HMX, y = self.HT_HMY, pen = None, symbol = 'o', brush = 'b')
        self.HT_RV_HM.addPoints(x = [float(hit.mz)], y = [float(hit.rt)], symbol = 'o', brush = 'r')

        return

    def add_hit_to_accepted_list(self):
        ''' take highlighted row from hit list and add to accepted hits '''

        # find highlighted row
        highlighted_row = self.HT_RV_hitlist.selectionModel().selectedRows()[0].row()
        hitNumber = int(self.HT_RV_hitlist.item(highlighted_row, 0).text())

        for hit in self.hits:
            if hit.index == hitNumber:
                rowCount = self.HT_RV_accepted_list.rowCount()
                self.HT_RV_accepted_list.insertRow(rowCount)
                self.HT_RV_accepted_list.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(hit.index)))
                self.HT_RV_accepted_list.setItem(rowCount, 1, QtGui.QTableWidgetItem(str(hit.rt)))
                self.HT_RV_accepted_list.setItem(rowCount, 2, QtGui.QTableWidgetItem('%.4f' %float(hit.mz)))
                self.HT_RV_accepted_list.setItem(rowCount, 3, QtGui.QTableWidgetItem('%.1f' %float(hit.score)))

                # centre text in tablewidget columns
                self.HT_RV_accepted_list.item(rowCount, 0).setTextAlignment(QtCore.Qt.AlignCenter)
                self.HT_RV_accepted_list.item(rowCount, 1).setTextAlignment(QtCore.Qt.AlignCenter)
                self.HT_RV_accepted_list.item(rowCount, 2).setTextAlignment(QtCore.Qt.AlignCenter)
                self.HT_RV_accepted_list.item(rowCount, 3).setTextAlignment(QtCore.Qt.AlignCenter)
                break

        # update table measurements to fit data
        self.HT_RV_accepted_list.resizeRowsToContents()
        return

    def remove_hit_from_accepted_list(self):
        # find highlighted row
        rows = self.HT_RV_accepted_list.selectionModel().selectedRows()

        if len(rows) == 0: return

        highlighted_row = rows[0].row()
        self.HT_RV_accepted_list.removeRow(highlighted_row)

        return

    def write_accepted_to_file(self):
        # get output file name and location
        ofname, _ = QtGui.QFileDialog.getSaveFileName(self, 'Set Output File', self.fileOpenDialogPath)

        if ofname == '': return

        # write data
        of1 = open(str(ofname), 'wt')

        of1.write('#- Validated HiTIME hit list:\n')

        # write postprocessing headers to results file
        of1.write('#- Postprocessing parameters\n')
        for k, v in self.headers.items():
            of1.write('# %s::: %s\n'%(k,v))
        of1.write('# validated::: 1\n') # indicates that these results have been user-checked

        of1.write('#-\n')

        of1.write('## mz, rt, score\n')
        # get data from rows of accepted hit table
        for r in range(self.HT_RV_accepted_list.rowCount()):
            hitNumber = int(str(self.HT_RV_accepted_list.item(r, 0).text()))
            rt = str(self.HT_RV_accepted_list.item(r, 1).text())
            mz = str(self.HT_RV_accepted_list.item(r, 2).text())
            score = str(self.HT_RV_accepted_list.item(r, 3).text())
            of1.write('%s, %s, %s\n' %(mz, rt, score))
        of1.close()
        return

    def getXRange(self, hit, window, rangeType = 'pc'):
        '''
        Get the xrange values at +/- Xpc of a target value
        '''
        hit = float(hit)
        window = float(window)

        if rangeType == 'pc':
            xmin = hit - hit * window/100
            xmax = hit + hit * window/100
            return (xmin, xmax)

        elif rangeType == 'fixed':
            center = hit + self.mzDelta/2
            xmin = center - window
            xmax = center + window
            return (xmin, xmax)

    def getYRange(self, xData, yData, xRange):

        xmin, xmax = xRange

        # find indices of xData entries within xrange
        ind = np.where((xData > xmin) & (xData < xmax))

        # find y values corresponding to these points
        yDataInXR = yData[ind]

        # return max Y value
        return float(np.max(yDataInXR))

    def reset(self):
        # clear contents of tablewidgets

        # clear hit list
        for i in reversed(range(self.HT_RV_hitlist.rowCount())):
            self.HT_RV_hitlist.removeRow(i)

        # clear accepted list
        for i in reversed(range(self.HT_RV_accepted_list.rowCount())):
            self.HT_RV_accepted_list.removeRow(i)

        # clear plots
        self.HT_RV_HM.clear()
        self.HT_RV_MS.clear()
        self.HT_RV_EIC.clear()

        return
