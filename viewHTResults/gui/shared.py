import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

'''
HiTIME tab classes
'''
# class filenames(object):
# 	''' class handling HT I/O filenames '''
# 	def __init__(self, ifpath, ftype):
# 		ofname, ofpath, ifname = self.get_fnames(ifpath)
# 		self.ifpath = str(ifpath)
# 		self.ifname = str(ifname)
# 		self.ofname = str(ofname)
# 		self.ofpath = str(ofpath)
# 		self.ftype = str(ftype)
#
# 	def __repr__(self):
# 		return( 5* '%s, ' %(self.ftype, self.ifname, self.ifpath, self.ofname, self.ofpath))

class process_parameters(object):
	''' prepare input data for running hitime'''
	def __init__(self, format, intensityRatio, rtWidth, rtSigma, ppm, mzWidth, mzSigma, inputFile, outputFile, logFile, mzDelta, removeLow, outDir, noScore, minSample):

		self.format = format
		self.intensityRatio = intensityRatio
		self.rtWidth = rtWidth
		self.rtSigma = rtSigma
		self.ppm = ppm
		self.mzWidth = mzWidth
		self.mzSigma = mzSigma
		self.inputFile = inputFile
		self.outputFile = outputFile
		self.logFile = logFile
		self.mzDelta = mzDelta
		self.removeLow = removeLow
		self.outDir = outDir
		self.noScore = noScore
		self.minSample = minSample

	def __repr__(self):
		return(15 * '%s, ' %(self.format, self.intensityRatio, self.rtWidth, self.rtSigma, self.ppm, self.mzWidth, self.mzSigma, self.inputFile, self.outputFile, self.logFile, self.mzDelta, self.removeLow, self.outDir, self.noScore, self.minSample))

'''
NTPS tab classes
'''
class cStructure(object):
	def __init__(self, structure_type, smile, img):
		self.structure_type = structure_type
		self.smile = smile
		self.img = img # rel path to saved png image file
		self.get_structure_data()

	def get_structure_data(self):
		self.mol = Chem.MolFromSmiles(self.smile)
		self.formula_dict = self.get_formula(Chem.rdMolDescriptors.CalcMolFormula(self.mol))
		self.formula_plain = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
		self.MW = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)

	def get_formula(self, f_string):
		print('f_string is: %s' %f_string)
		self.f_string = f_string.replace('*','')
		i = re.findall(r'([A-Z][a-z]*)(\d*)', f_string)

		formula = {}
		for x in i:

			if x[1] != '':
				formula[x[0]] = x[1]
			else:
				formula[x[0]] = 1

		self.formula = formula
		return formula

	def __repr__(self):
		return (5*'%s, ' %(self.structure_type, self.smile, self.MW, self.formula_plain, self.formula_dict))

class element_ranges(object):
	def __init__(self, mols):
		self.min_mz, self.max_mz, self.atom_dict = self.get_atom_dict(mols)

	def get_atom_dict(self, mols):
		atom_dict = {}
		min_mass = None
		max_mass = None
		for mol in mols:
			#print(mol)
			if min_mass is not None:
				if float(mol.MW) < min_mass:
					min_mass = float(mol.MW)
			else: min_mass = float(mol.MW)

			if mol.structure_type == 'precursor':
				max_mass = mol.MW

			for key, value in mol.formula_dict.iteritems():
				#print('key, value are:')
				#print(key, value)
				if key in atom_dict:
				#	print('Key found')
				#	print(key, value)
				#	print('atom dict entry is')
					ad_min =  atom_dict[key][0]
					ad_max =  atom_dict[key][1]
					if int(value) < int(atom_dict[key][0]):
						atom_dict[key][0] = int(value)
					if int(value) > int(atom_dict[key][1]):
						atom_dict[key][1] = int(value)
				else:
					atom_dict[key] = [int(value),int(value)] # min / max
		return min_mass, max_mass, atom_dict

	def __repr__(self):
		return (3 * '%s, ' %(self.min_mz, self.max_mz, self.atom_dict))
