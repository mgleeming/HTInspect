import sys
import time
import os, re, math
import numpy as np
from rtree import index
import pyqtgraph as pg

from PyQt4 import QtCore, QtGui
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

def launch_HT_search(q, args):
	try:
		for arg in args:

			q.put ('Processing file %s using HT-CPP with %s threads' %(arg['inputFile'], arg['threads']))

			# assemble HT-CPP exec line
			execLine = 'hitime-score'

			execLine += ' -i %s' %arg['inputFile']
			execLine += ' -o %s' %(os.path.splitext(arg['outputFile'])[0] + '.mzML')

			execLine += ' -j %s' %arg['threads']
			execLine += ' -a %s' %arg['intensityRatio']
			execLine += ' -r %s' %arg['rtWidth']
			execLine += ' -t %s' %arg['rtWidth']
			execLine += ' -m %s' %arg['mzWidth']
			execLine += ' -z %s' %arg['mzSigma']
			execLine += ' -d %s' %arg['mzDelta']

			proc = Popen(
				execLine,
				stdout = PIPE,
				stderr = PIPE,
				shell = True
				)

			output, error = proc.communicate()

			if len(error.strip()) == 0:
				q.put('Finished processing file: %s' %arg['inputFile'])
				q.put('done')
				return
			else: raise Exception(error)

	except Exception, e:
		q.put('HT-CPP exception: %s' %e)
		print str(e)
		try:
			# run regular HT search
			import viewHTResults.hitime.hitime_methods as HTM
			HTM.gui_init(q, args)

		# it's all over
		except Exception, e: print str(e)
	return


def read_hitime_files(inFile, returnNp = False, **kwargs):
	# get hitime file type
	# options are:
	#	- raw data python
	# 	- raw data cpp
	# 	- subtracted data
	# 	- validated hist from gui
	# 	- invalid

	'''
	raw data python
		4 columns	- rt, mz, amp, score

	raw data cpp
		3 columns	- mz, rt, score

	validated data
		headers
			- begin with '#'
			-
		data formatted as mz, rt, score

	subtracted data
		headers -


	general approach
	----------------------------------

	use numpy np.getfromtxt function
		has arguments: delimiter, names (column names), skip_header .... + more


	if headers encountered in first x lines
		- use this infor to parameterise numpy matrix
		- count number of header lines

	read data into numpy array

	sort by score

	apply score cutoff

	create hitime_hit objects if needed

	'''

	htHeaders = dict()

	if '.mzML' in inFile:
		print 'reading mzML file'

		spectra = readSpectra(inFile)

		Times, Mzs, Ints = None, None, None
		for spectrum in spectra:
			times, mzs, ints, lvl = spectrum
			times = np.repeat(times,mzs.shape[0])

			# data = np.vstack((time, mzs, ints))
			if Times is not None:
				Times = np.hstack((Times, times))
				Mzs = np.hstack((Mzs, mzs))
				Ints = np.hstack((Ints, ints))
			else:
				Times = times
				Mzs = mzs
				Ints = ints

		data = [Times, Mzs, Ints]
		npData = np.core.records.fromarrays(
											data,
											names = 'rt, mz, score',
											formats = '<f8, <f8, <f8'
											)
	else:
		print 'reading text file'
		f = open(str(inFile),'r')
		# count header rows
		headerRows = 0
		htHeaders = dict()
		while True:
			line = f.readline().strip()
			if line.startswith('#') and not line.startswith('##') and not line.startswith('#-'):
				headerRows += 1 # parameter - value pair
				k,v = line.replace('#','').replace(' ','') .split(':::')
				htHeaders[k] = v
				print k,v

			elif line.startswith('#-'): headerRows += 1 # blank/info only line

			elif line.startswith('##'): # column header line
				names = line.replace('##','').replace(' ','').split(',')
				usecols = (0,1,2)
			else:
				# read until the first non-header line
				if ',' in line:	delimiter = ','
				f.seek(0)
				break

		if headerRows == 0:
			# raw data file
			# - get raw data type

			while True:
				line = f.readline().strip()
				if line != '': break

			# get delimiter
			if ',' in line: delimiter = ','
			else: delimiter = ' '

			data = [x for x in line.replace(',', ' ') if x != '']

			if len(data) == 3:
				htRawDataType = 'cpp'
				names = ['mz', 'rt', 'score']
				usecols = (0,1,2)
			else:
				htRawDataType = 'py'
				names = ['rt', 'mz', 'amp', 'score']
				usecols = (0,1,3)
		# get data and place in numpy array
		f.seek(0)

		# TODO >> assign data types here
		npData = np.genfromtxt(
								f,
								delimiter = delimiter,
								names = names,
								skip_header = headerRows,
								usecols = usecols,
								comments = '#', # skip rows starting with #
							)
		f.close()

	scoreCutoff = kwargs.get('scoreCutoff', None)

	# apply score cutoff
	if scoreCutoff:
		npData  = npData[np.where(npData['score'] > scoreCutoff)]

	# sort entries in order of score - highest first
	npData = npData[npData['score'].argsort()[::-1]]

	# note:
	# to sort lowest first
	# npData[npData['score'].argsort()]

	retrieveTop = kwargs.get('retrieveTop', None)
	if retrieveTop: npData = npData[:retrieveTop]

	# return numpy array if requested
	if returnNp: return htHeaders, npData

	# dictionary specifying hitime_hit attributes other than rt,mz,score
	# these will be applied uniformly to all instances
	addHtAttributes = kwargs.get('addHtAttributes', None)

	# make list of hitime hit objects
	data = []
	for x in npData:
		hit = hitime_hit(x['mz'], x['rt'], x['score'])

		if addHtAttributes:
			for k,v in addHtAttributes.iteritems(): setattr(hit, k, v)

		data.append(hit)

	# clean up
	del npData


	return htHeaders, data

def get_HT_regions(HT_data, Drt, Dmz, rtExclusion = 0, minScore = None):

	count = 0
	HT_regions = []

	idx = index.Index()

	if minScore:
		HT_data = [x for x in HT_data if x.score > minScore]

	if Drt > 0 and Dmz > 0:
		for x in HT_data:
			mz = x.mz
			rt = x.rt
			score = x.score
			coord = (rt-Drt, mz-Dmz, rt+Drt+rtExclusion, mz+Dmz)
			if idx.count(coord) == 0:
				idx.insert(count, coord)
				HT_regions.append(x)
				count += 1
	else:
		return HT_data

	return HT_regions

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

def geztTimeUnitsFromMZML(ifile):
	import pymzml
	msrun = pymzml.run.Reader(str(ifile))
	for spectrum in msrun:
		for element in spectrum.xmlTree:
			elements = element.items()
			for i, e in enumerate(elements):
				# is there a nicer way to do this?
				if e[0] == 'name' and e[1] == 'scan time' or e[1] == 'scan start time':
					unit = elements[i+1][1]
					return unit

def build_rt_index(mzml_file):
	import pymzml

	print 'building rt index for mzML file:'

	msrun = pymzml.run.Reader(str(mzml_file))
	times = []
	unit = None
	for counter, spectrum in enumerate(msrun):

		if counter == 0:
			for element in spectrum.xmlTree:
				elements = element.items()
				for i, e in enumerate(elements):
					if e[0] == 'name' and e[1] == 'scan time' or e[1] == 'scan start time':
						unit = elements[i+1][1]

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

def plot_EIC(results, options):

	mzml_file = str(options.mzmlFile)
	neutral_mod_delta = options.mzDelta
	EIC_width = options.eicWidth

	# build list of light EIC targets including ranges for heavy/light isotopes
	for result in results:
		if hasattr(result, 'res_type'):
			if result.res_type != 'direct': continue

		else:
			setattr(result, 'charge', 1)
			setattr(result, 'HT_mz', result.mz)

		result.light_ll = result.HT_mz - EIC_width
		result.light_hl = result.HT_mz + EIC_width
		result.heavy_ll = result.HT_mz - EIC_width + neutral_mod_delta / result.charge
		result.heavy_hl = result.HT_mz + EIC_width + neutral_mod_delta / result.charge
		result.EIC_rt = []
		result.EIC_int_light = []
		result.EIC_int_heavy = []

	# get spectra - returns tuple of time, mzs, ints
	spectra = readSpectra(mzml_file, 1)

	last = None

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

		if last is not None and n-1 != last:
			print 'error in spectral indexing'
			print n, last
			sys.exit()
		last = n
	print('Done!')

	return results

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

def readSpectra(mzml_file, msLevel = None):
	'''
	read MS_spectra from an mzML file of the given msLevel
	'''
	import pymzml

	unit = geztTimeUnitsFromMZML(mzml_file)
	print unit
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
			except Exception, e:
				print 'Warning, skipping spectrum %s' %n
				print 'Stack trace:'
				print str(e)
				continue

		try:
			mzs = np.array(spectrum.mz, dtype = "float32")
			ints = np.array(spectrum.i, dtype = 'float32')
			if 'second' in unit:
				time /= 60 # HTCPP output files are in seconds
			assert mzs.shape == ints.shape
			yield time, mzs, ints, lvl

		except Exception, e:
			print 'Warning, skipping spectrum %s' %n
			print 'Stack trace:'
			print str(e)
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
