# IS Calculation Script
import os, glob, sys, math
import msa, coe, iscalc
import numpy as np

def get_subst_list(msa_path):
	f = open(msa_path)
	data = f.readlines()
	f.close()
	query = data[3].strip().split()[0]
	seq = "".join(map(lambda x: x.strip().split(' ')[-1] if query in x else '', data))
	seq = str(filter(lambda x: x!='-', seq))

	subst_list = []
	AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for i, AA1 in enumerate(seq):
		for AA2 in AA_list:
			if AA1==AA2: continue
			subst_list.append("%s%s%s" %(AA1, i+1, AA2))
	return subst_list

def build_coevolution_score_dic(prefix, cn_kind):
    output_dic = {}
    f = open("./output/%s/%s.coe_out_%s" %(prefix, prefix, cn_kind))
    f.readline()
    for line in f.xreadlines():
        line = line.strip().split("\t")
        output_dic[(line[0], line[1])] = float(line[2])
    f.close()
    return output_dic

def get_percentile_rank(aList, bReverse=False):
    step = 1.0 / (len(aList) - 1)
    bList = []
    for v in aList:
        bList.append(v)
    bList.sort(reverse=bReverse)
    rDic, cur = ({}, 0)
    for v in bList:
        if v not in rDic:
            rDic[v] = cur 
        cur += step
    rList = []
    for v in aList:
        rList.append(rDic[v])
    return (rList, rDic)

def calc_CN(sorted_key, len_seq, len_threshold):
    CN_dic = {}
    for res in range(1, len_seq+1):
        CN_dic[res] = 0 
    for i in range(int(len_seq*len_threshold)):
        res1, res2 = sorted_key[i]
        CN_dic[int(res1)] += 1
        CN_dic[int(res2)] += 1
    values = CN_dic.values()
    tmp_list, per_dic = get_percentile_rank(values)
    per_CN_dic = {}
    for res in CN_dic.keys():
        per_CN_dic[res] = per_dic[CN_dic[res]]
    return CN_dic, per_CN_dic

def build_coupling_dic(sorted_key, len_seq, len_threshold):
    coupling_dic = {}
    for res in range(1, len_seq+1):
        coupling_dic[res] = []
    for i in range(int(len_seq*len_threshold)):
        res1, res2 = sorted_key[i]
        coupling_dic[int(res1)].append(int(res2))
        coupling_dic[int(res2)].append(int(res1))
    return coupling_dic

def build_msa(msa_path):
	prefix = msa_path.split("/")[-1].split(".")[0]
	tmp_dic = {}
	f = open(msa_path)
	for line in f.xreadlines():
		line = line.strip().split()
		if len(line) == 2:
			if prefix in line[0]: line[0] = prefix
			if not tmp_dic.has_key(line[0]):
				tmp_dic[line[0]] = ""
			tmp_dic[line[0]] += line[1]
	f.close()
	msa_dic = {}
	gene = prefix.split("_")[0]
	for i in range(len(tmp_dic[gene])):
		if tmp_dic[gene][i] != "-":
			for key in tmp_dic.keys():
				if not msa_dic.has_key(key):
					msa_dic[key] = ""
				msa_dic[key] += tmp_dic[key][i]
	return msa_dic

def calc_ent_diff_coupled(msa_dic, gene, coupling_dic, residue, AA1, AA2, cn_kind="McBasc"):
	try:
		if msa_dic[gene][residue-1] == AA1:
			ent_diff_list = []
			for residue2 in coupling_dic[residue]:
				residue_list = map(lambda x: x[residue-1]+x[residue2-1], msa_dic.values())
				residue_list = filter(lambda x: not("-" in x), residue_list)
				n_AA1 = residue_list.count(AA1 + msa_dic[gene][residue2-1])
				n_AA2 = residue_list.count(AA2 + msa_dic[gene][residue2-1])
				entropy_difference = -math.log((n_AA2+1)/float(n_AA1))
				ent_diff_list.append(entropy_difference)
			return np.mean(ent_diff_list)
		else:
			print "Not accurate information; The AA in %s (%d) is %s not %s" %(gene, residue, msa_dic[gene][residue-1], AA1)
			return np.nan
	except TypeError:
		return np.nan

def get_CN(base_pth, pm, prefix):
	# Get CN (default: McBASC)
	coe_output_mcbasc = os.path.join(base_pth, "%s.coe_out_mcbasc" % prefix)
	pcn = iscalc.ProcCN(pm, coe_output_mcbasc, pm.result[0].id, cn_cutoff=2.0, coe_algorithm="McBASCCovariance")
	cn_result = pcn.calc()
	tmp_coe = build_coevolution_score_dic(prefix, 'mcbasc')
	sorted_coe = sorted(tmp_coe, key=lambda x: tmp_coe[x], reverse=True)
	CN_dic, w_CN_dic = calc_CN(sorted_coe, pcn.query_len, 2.0)
	
	# Print CN results
	CNScores, query_len = (pcn.CNScores, pcn.query_len)
	f = open(os.path.join(base_pth, "%s.cn") %(prefix), 'w')
	print >> f, "\t".join(['res', 'CN_McBasc', 'w_CN'])
	for res in range(1, query_len+1):
		output = [res, CN_dic.get(res,0), w_CN_dic.get(res, 0)]
		print >> f, "\t".join(map(str, output))
	f.close()

	# Build and extract coevolutionary network
	coenet_dic = build_coupling_dic(sorted_coe, pcn.query_len, 2.0)
	f = open(os.path.join(base_pth, "%s.coenet") %(prefix), 'w')
	print >> f, "res1\tres2"
	for res1 in coenet_dic.keys():
		for res2 in coenet_dic[res1]:
			if res1 > res2:
				print >> f, "\t".join(map(str, [res1, res2]))
	f.close()
	return CN_dic, w_CN_dic, coenet_dic

def calc_SCI(base_pth, subst_list, prefix, msa_dic, CN_dic, w_CN_dic, coenet_dic):
	SCI_list = []
	for subst in subst_list:
		AA1, AA2, res = subst[0], subst[-1], int(subst[1:-1])
		avg_ent_diff = calc_ent_diff_coupled(msa_dic, prefix.split("_")[0], coenet_dic, res, AA1, AA2)
		cn, w_cn = CN_dic.get(res,0), w_CN_dic.get(res,0)
		sci = w_cn*avg_ent_diff*(-1)
		SCI_list.append([res, AA1, AA2, sci, cn])
	sorted_SCI_list = sorted(filter(lambda x: ~np.isnan(x[3]), SCI_list), key=lambda x: x[3], reverse=True)
	fo = open(os.path.join(base_pth, "%s.txt" %prefix), 'w')
	print >> fo, "\t".join(["res", "AA1", "AA2", "SCI", "Number of co-evolutionary relationships"])
	for line in sorted_SCI_list: print >> fo, "\t".join(map(str, line))
	fo.close()

