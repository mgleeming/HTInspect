import os, re, sys, operator
from datetime import datetime
import numpy as np
import argparse

from HTInspect.res.utils import *
import HTInspect.hitime.HTS_resutls_file_parser as HTS_FP

parser = argparse.ArgumentParser(description = 'find hits in a HiTIME reuslts file')

parser.add_argument('-htIn',
                    help = 'HiTIME results file or peak list for processing',
                    required = True,
                    type = str)
parser.add_argument('-outFile',
                    help = 'Output file name',
                    type = str)
parser.add_argument('--mzmlFile',
                    help = 'mzML file used in HiTIME searching',
                    type = str)
parser.add_argument('--peakList',
                    help = 'Input file provided is a list of HT peaks',
                    action = 'store_true')
parser.add_argument('--mzDelta',
                    help = 'mass difference between heavy and light isotopes. Default = 6.0201',
                    default = 6.0201,
                    type = float)
parser.add_argument('--rtWidth',
                    help = 'retention time width for local maxima detection',
                    type = float)
parser.add_argument('--mzWidth',
                    help = 'm/z width for local maxima detection',
                    type = float)
parser.add_argument('--eicWidth',
                    help = 'm/z width used for plotting Extracted Ion Chromatograms. Default = 0.03',
                    default = 0.03,
                    type = float)
parser.add_argument('--scoreCutoff',
                    help = 'HiTIME score threshold. Default = 0',
                    default = 0,
                    type = float)
parser.add_argument('--rtExclusion',
                    help = 'Exclude points for specified time following accepted hit',
                    default = 0,
                    type = float)


def write_headers(options, unit, of1):
    of1.write('#- BEGIN_POSTPROCESSING_PARAMETERS\n')
    of1.write('# htFile::: %s\n' %options.htIn)
    of1.write('# mzmlFile::: %s\n' %options.mzmlFile)
#    of1.write('# mzWidth::: %s\n' %options.mzWidth)
#    of1.write('# rtWidth::: %s\n' %options.rtWidth)
    of1.write('# minIntensity::: %s\n' %options.minIntensity)
    of1.write('# mzDelta::: %s\n' %options.mzDelta)
    of1.write('# eicWidth::: %s\n' %options.eicWidth)
    of1.write('# mzML_time_unit::: %s\n' %unit)
    of1.write('#- END_POSTPROCESSING_PARAMETERS\n')
    return

def write_ht_hits_to_file(results, of1):
    of1.write("#- BEGIN_RAW_HITIME_RESULTS\n")
    of1.write('## rt, mz, amp\n')
    for i in results:
        of1.write('%s, %s, %s\n' %(i.rt, i.mz, i.amp))
    of1.write("#- END_RAW_HITIME_RESULTS\n\n")
    return

def write_data_to_file(results, of1):
    of1.write('BEGIN_HITIME_HIT_RESULTS\n')
    for index, i in enumerate(results):
        of1.write('\nBEGIN_Results_for_hit %s\n' %index)
        of1.write('<RT> %.2f\n' %i.rt)
        of1.write('<MZ> %s\n' %i.mz)
        of1.write('<SCORE> %.2f\n' %i.amp)
        of1.write('<EIC_RT> %s\n' %['%.3f' %x for x in i.EIC_rt])
        of1.write('<EIC_int_light> %s\n' %i.EIC_int_light)
        of1.write('<EIC_int_heavy> %s\n' %i.EIC_int_heavy)
        of1.write('<MS_mz> %s\n' %[('%.4f' %x) for x in i.HT_MS_mz])
        of1.write('<MS_int> %s\n' %[int(x) for x in i.HT_MS_int])
        of1.write('<END_Results_for_hit> %s\n' %index)
    return

def build_rt_index(mzml_file):
    import pymzml

    print ('building rt index for mzML file:')

    msrun = pymzml.run.Reader(str(mzml_file))
    times = []
    unit = None
    for counter, spectrum in enumerate(msrun):
        level = spectrum['ms level']

        try:
            time = spectrum['scan time']
        except:
            try:
                time = spectrum['scan start time']
            except:
                if 'total ion current chromatogram' in spectrum.keys(): continue
                else:
                    raise Exception('mzML spectrum read error')
                continue
        times.append([level, time])

    times = np.asarray(times)

    return times, unit

def resolve_rt_targets_to_spectra(results, rtIndexArray = None, msLevel = None, mstype = None):
    '''
    Find RT of spectrum in mzML file used for HT scanning that matches a given HT
    '''

    #print rtIndexArray.shape
    for count, res in enumerate(results):

        if hasattr(res, 'HT_rt') or hasattr(res, 'pep_rt'):
            # for targeted search
            if msLevel == 1:
                ht_rt = res.HT_rt
            elif msLevel == 2:
                ht_rt = float(res.pep_rt)

        else:
            # for simple HT results postprocessing
            ht_rt = res.rt

        # get indices of all rows containing msLevel spectra
        msindices = np.where((rtIndexArray[:,0] > (msLevel - 0.1)) & (rtIndexArray[:,0] < (msLevel + 0.1)) )

        # create array of rows containing appropriate msLevel scans
        msrts = rtIndexArray[msindices]

        # get absolute value of difference between these RTs and ht_rt
        abs_min = np.absolute(msrts - ht_rt)

        # index of minimum value
        min_rt =  np.argmin(abs_min, axis = 0)[1]

        # get minimum value
        min_rt_val = msrts[min_rt]

        # get index of this value in original rtIndexArray
        full_array_index = np.where(rtIndexArray[:,1] == msrts[min_rt][1])[0][0]
        full_array_value = rtIndexArray[full_array_index]

        if mstype == 'MS1':
            # assign to res.rt
            #res.rt = min_rt # might not actually want to change this - not necessary anyway
            res.ht_rt_index = full_array_index

        elif mstype == 'MS2':
            # assign the index of matching ms2 spectra to pep_r
            # print 'MS2 correlation stats'
            # print 'rtIndexArray size is:'
            # print rtIndexArray.shape
            # print full_array_index, full_array_value, res.HT_rt, min_rt_val, min_rt
            res.pep_rt_index = full_array_index

            search = 200
            for i in xrange(1, search):

                prev_index = full_array_index - i
                next_index = full_array_index + i

                # test if either of these are ms1 spectra
                if rtIndexArray[next_index][0] == 1:
                    res.nearest_ms1 = next_index
                    break
                else:
                    pass

    return results

#@profile
def main(options):

    # parse raw HT datafile
    #results = sort_HT_results(options, of1)
    print ('Reading HiITME file')
    results = HTS_FP.reader(options.htIn)
    results = [x for x in results if x.amp > options.minIntensity]

    print ('%s HiTIME hits detected' %len(results))

    ms1RTIndex, unit = build_rt_index(options.mzmlFile)

    results = resolve_rt_targets_to_spectra(
                                    results,
                                    rtIndexArray = ms1RTIndex,
                                    msLevel = 1,
                                    mstype = 'MS1'
                                    )

    # extract EIC and MS spectrumfor each hit
    results = plot_EIC(results, options)

    print ('Writing data to file')
    # write data to file

    # Create output file
    of1 = open(options.outFile, 'wt')

    write_headers(options, unit, of1)

    # write HT htis to file
    write_ht_hits_to_file(results, of1)

    # write hit MS and EIC data to file
    write_data_to_file(results, of1)

    of1.close()

    print ( 'Finished postprocessing')
    return

class GuiOptions(object):
    def __init__(self, args):
        for k,v in args.items():
            setattr(self, k, v)

def guiRun(q, args):
    options = GuiOptions(args)
    main(options)
    q.put('done')
    return

if __name__ == '__main__':
    options = parser.parse_args()
    sys.exit(main(options))
