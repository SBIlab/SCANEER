import os, sys
from __init__ import lib_path

def check_conf_file():
	f = open("./SCANEER/coe/Energetics.properties")
	line_list = f.readlines()
	HOME_DIR = line_list[4].strip().split("=")[1]
	f.close()
	if not os.path.isdir(HOME_DIR):
		new_HOME_DIR = os.path.abspath("./SCANEER/coe")
		fo = open("./SCANEER/coe/Energetics.properties", 'w')
		for i, line in enumerate(line_list):
			if i == 4: print >> fo, "HOME_DIRECTORY=%s" %new_HOME_DIR
			else: print >> fo, line.strip()
		fo.close()

check_conf_file()

class ProcCoe:
	# input aln, out coe file path are required
	def __init__(self, aln_file, coe_file, algorithm="McBASC"):
		self.algorithm = algorithm
		self.input_aln_file = os.path.abspath(aln_file)
		self.output_coe_file = os.path.abspath(coe_file)
		self.result = []

	# Available Algorithm: McBASC, ELSC, SCA, MI, OMES
	def run(self):
		if sys.platform == 'win32': # Windows Environment
			jexePath = 'java.exe'
			cmd_apd = ' & '
		else: # Linux Environment
			jexePath = os.path.join(lib_path, "jre/linux/bin/java")
			jexePath = "/usr/bin/java"
			cmd_apd = ';'
		jClassPath = os.path.join(lib_path, "coe")
		jClassAlgorithm = "covariance.algorithms.%s" % self.algorithm
		exejClass = "%s %s %s %s" % (jexePath, jClassAlgorithm, self.input_aln_file, self.output_coe_file) 
		os.system(cmd_apd.join(['cd %s' % jClassPath, exejClass]))
		
	# Parse coe-calculation result
	def parse(self):
		f = open(self.output_coe_file, 'r')
		f.next()
		for line in f.xreadlines():
			fields = line.split()
			self.result.append((int(fields[0]), int(fields[1]), float(fields[2])))
		f.close()
		return self.result

	def convertResPos(self, res_pos_dic):
		f = open(self.output_coe_file, 'w')
		print >> f, "residue_1\tresidue_2\tscore"
		for cell in self.result:
			print >> f, '%d\t%d\t%s' %\
			 (res_pos_dic[cell[0]]+1, res_pos_dic[cell[1]]+1, str(cell[2]))
		f.close()
			
