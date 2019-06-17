import sys
import time
import os, re, math
import numpy as np
import pyqtgraph as pg

from PyQt5 import QtCore, QtGui
from multiprocessing import Process, Queue
from subprocess import Popen, PIPE, STDOUT, call

class Runner(QtCore.QObject):
    '''
    Runs a job in a separate process and forwards messages from the job to the
    main thread through a pyqtSignal.
    '''
    msg_from_job = QtCore.pyqtSignal(object)

    def __init__(self, start_signal):
        '''
        :param start_signal: the pyqtSignal that starts the job
        '''
        super(Runner, self).__init__()
        self.job_input = None
        self.job_function = None
        self.comm_queue = None
        start_signal.connect(self._run)

    def _run(self):
        self.p = Process(target=self.job_function, args=(self.comm_queue, self.job_input,))
        self.p.start()

class hitime_hit (object):
    def __init__(self, mz, rt, score):
        self.mz = mz
        self.rt = rt
        self.score = score
    def __repr__(self):
        return ( 3 * '%s, ' %(self.mz, self.rt, self.score))

def run_job(runner, runner_thread, queue, function, input):
    """ Call this to start a new job """
    runner.job_function = function
    runner.job_input = input
    runner.comm_queue = queue
    runner_thread.start()

def prepare_HM_color_grad(data):

    pos = np.array([0.0, 1.0])

    color = np.array([[255,255,255,1], [255,0,0,1]], dtype=np.ubyte)
    a = pg.ColorMap(pos, color)
    lut = a.getLookupTable(0.0, 1.0, 256)

    points = np.asarray([x.score for x in data])
    #	dmax = np.max(points)

    points = pg.applyLookupTable(points, lut)

    dmax0 = np.max(points[:,0])
    dmax1 = np.max(points[:,1])
    dmax2 = np.max(points[:,2])
    pens = []
    for i in points:
        r, g, b = float(i[0])/dmax0*255, float(i[1])/dmax1*255, float(i[2])/dmax2*255
        colour = pg.mkColor(r,g,b)

        #print i, r, g, b
        pens.append(colour)

    return pens

def plot_EIC(results, options):

    mzml_file = str(options.mzmlFile)
    neutral_mod_delta = options.mzDelta
    EIC_width = options.eicWidth

    # build list of light EIC targets including ranges for heavy/light isotopes
    for result in results:
        result.light_ll = result.mz - EIC_width
        result.light_hl = result.mz + EIC_width
        result.heavy_ll = result.mz - EIC_width + neutral_mod_delta
        result.heavy_hl = result.mz + EIC_width + neutral_mod_delta
        result.EIC_rt = []
        result.EIC_int_light = []
        result.EIC_int_heavy = []

    # get spectra - returns tuple of time, mzs, ints
    spectra = readSpectra(mzml_file, 1)

    #print results[0].ht_rt_index
    print('Extracting EIC...')
    for n, spectrum in enumerate(spectra):
        time, mzs, ints, lvl = spectrum

        # test if spectrum corresponds to HT hit point
        for res in results:
            res.EIC_rt.append(float(time))
            res.EIC_int_light.append(0)
            res.EIC_int_heavy.append(0)

            if res.ht_rt_index == n:
                # print 'found spectrum: spec rt: %s, result rt: %s' %(time, res.rt)
                # if so, record data
                res.HT_MS_mz = mzs
                res.HT_MS_int = ints

        for res in results:
            # get mzs indices within window for light and heavy EICs
            light = np.where((mzs > res.light_ll) & (mzs < res.light_hl))
            heavy = np.where((mzs > res.heavy_ll) & (mzs < res.heavy_hl))

            # check that the array is non-empty
            if light[0].shape[0] > 0: res.EIC_int_light[-1] += np.sum(ints[light])
            if heavy[0].shape[0] > 0: res.EIC_int_heavy[-1] += np.sum(ints[heavy])

    return results

def readSpectra(mzml_file, msLevel = None):
    '''
    read MS_spectra from an mzML file of the given msLevel
    '''
    import pymzml

    msrun = pymzml.run.Reader(str(mzml_file))
    for n, spectrum in enumerate(msrun):

        # only consider MS1 level
        if msLevel:
            if spectrum['ms level'] != msLevel: continue

        lvl = spectrum['ms level']

        try:
            time = spectrum['scan time']
        except:
            try:
                time = spectrum['scan start time']
            except Exception as e:
                print ('Warning, skipping spectrum %s' %n)
                print ('Stack trace:')
                print (str(e))
                continue

        try:
            mzs = np.array(spectrum.mz, dtype = "float32")
            ints = np.array(spectrum.i, dtype = 'float32')
            assert mzs.shape == ints.shape
            yield time, mzs, ints, lvl

        except Exception as e:
            print ('Warning, skipping spectrum %s' %n)
            print ('Stack trace:')
            print (str(e))
            continue


def getXRange(hit, pc):
    '''
    Get the xrange values at +/- Xpc of a target value
    '''

    hit = float(hit)
    pc = float(pc)

    xmin = hit - hit * pc/100
    xmax = hit + hit * pc/100
    return (xmin, xmax)

def retrieve_file_matching_string(string, files):
    '''
    Input: list of filename objects
    '''

    for f in files:
        if f.ifname == str(string):
            return f.ifpath

def zero_fill(xData, yData):

    x = np.repeat(xData, 3)
    y = np.dstack((np.zeros(yData.shape[0]), yData, np.zeros(yData.shape[0]))).flatten()

    return x, y
