from SCANEER.SCANEER import *
from SCANEER import msa, coe, iscalc

input_path = "./input/"
output_path = "./output/"
if not os.path.isdir(output_path):
	os.makedirs(output_path)

file_list = glob.glob(os.path.join(input_path, "*.aln"))
for msa_path in file_list:
	msa_file = os.path.basename(msa_path)
	prefix = msa_file.split('.')[0]

	base_pth = os.path.join(output_path, prefix)
	if not os.path.isdir(base_pth):
		os.makedirs(base_pth)
	elif len(os.listdir(base_pth)) == 5:
		print "already calculated! skip!"
		continue

	# Getting all single AA substitution list
	subst_list = get_subst_list(msa_path)

	# Loading MSA
	pm = msa.ProcMsa("tmp", msa_path, "tmp", "tmp")
	pm.parse()
	msa_dic = build_msa(msa_path)
	
	# Analyzing co-evolutionary network
	CN_dic, w_CN_dic, coenet_dic = get_CN(base_pth, pm, prefix)

	# Calculating SCI
	calc_SCI(base_pth, subst_list, prefix, msa_dic, CN_dic, w_CN_dic, coenet_dic)

