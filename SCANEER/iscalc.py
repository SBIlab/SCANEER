import sys, os
from Bio import AlignIO
import coe

# Global Function: get percentile rank
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
        
# Process of CN-Calculation        
class ProcCN:
    def __init__(self, proc_msa, coe_out_path, query_id, cn_cutoff=2.0, coe_algorithm="McBASC"):
        # Options (derivated from ProcMsa)
        self.proc_msa = proc_msa 
        self.aln_file_path = proc_msa.output_aln
        self.alignment = proc_msa.parse()
        self.coe_out_path = coe_out_path
        if not os.path.isabs(coe_out_path):
            self.coe_out_path = os.path.join(os.path.abspath('.'), coe_out_path)
        
        # Options (coe and cn options)
        self.query_id = query_id
        self.cn_cutoff = cn_cutoff
        self.coe_algorithm = coe_algorithm
        
        # Output
        self.CNScores = None
        self.CNPercentiles = None
        
    def calc(self):
        # Get Percentile Value of CNScores    
        def percentile_cn(self, raw_data, res_pos_dic, query_len):
            cn_dic = {}
            num_cutoff = int(self.cn_cutoff * query_len)
            
            raw_data.sort(key=lambda item:item[2], reverse=True)
            for (_ri, _rj, _) in raw_data[:num_cutoff]:
                ri, rj = (res_pos_dic[_ri], res_pos_dic[_rj])
                cn_dic[ri] = cn_dic.get(ri, 0) + 1
                cn_dic[rj] = cn_dic.get(rj, 0) + 1
                
            rList, rDic = get_percentile_rank(cn_dic.values())
            result = {}
            for a in cn_dic:
                result[a] = rDic[cn_dic[a]]
            return (cn_dic, result)                
        
        # calculation CN
        coe_in_path = os.path.join(os.path.dirname(self.coe_out_path), os.path.basename(self.aln_file_path) + '_cn')
        res_pos_dic, query_len = self.proc_msa.makeCNinput(coe_in_path, self.query_id, cutoff=0.2)
        Coe = coe.ProcCoe(coe_in_path, self.coe_out_path, self.coe_algorithm)
        Coe.run()
        cn_raw_data = Coe.parse()
        Coe.convertResPos(res_pos_dic)
        self.CNScores, self.CNPercentiles = percentile_cn(self, cn_raw_data, res_pos_dic, query_len)
        self.query_len = query_len
        return (self.CNScores, self.CNPercentiles)

