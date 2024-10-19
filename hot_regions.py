import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import pearsonr
import math
from scipy.stats import fisher_exact
import os
import json
import datetime
import time


class JsonEncoder(json.JSONEncoder):
    """Convert numpy classes to JSON serializable objects."""

    def default(self, obj):
        if isinstance(obj, (np.integer, np.floating, np.bool_)):
            return obj.item()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(JsonEncoder, self).default(obj)



def calc_distance(start1, end1, start2, end2):
    if end1 >= start2:
        return end1 - start2
    if end2 >= start1:
        return end2 - start1
    if end1 < start2 < start1:
        return start2 - end1
    if end2 < start1 < start2:
        return start1 - end2



class Gene:
    def __init__(self):
        self.type = ''
        self.chr = ''
        self.gene_id= ''
        self.start = 0
        self.end = 0
        self.strand = ''
        self.part = {'intron':[]}
        self.first_rna_data = {'id': '', 'start': '', 'end': ''}
        self.back_gene = {}
        self.next_gene = {}
        self.nearest_gene = {}

    def find_nearest_gene(self):
        if self.back_gene == {} and self.next_gene != {}:
            self.nearest_gene = self.next_gene
            distance = abs(self.nearest_gene['start'] - self.end)
            self.nearest_gene['distance'] = distance
        elif self.back_gene != {} and self.next_gene == {}:
            self.nearest_gene = self.back_gene
            distance = abs(self.nearest_gene['end'] - self.start)
            self.nearest_gene['distance'] = distance
        else:
            distance_back = abs(self.back_gene['end'] - self.start)
            distance_next = abs(self.next_gene['start'] - self.end)
            if distance_back < distance_next:
                self.nearest_gene = self.back_gene
                self.nearest_gene['distance'] = distance_back
            else:
                self.nearest_gene = self.nearest_gene
                self.nearest_gene['distance'] = distance_next






class Classify:

    def __init__(self, gfilename, pfilename, pefilename, zfile, vfile, model='Live'):
        """
        :param gfilename: gff3数据文件名
        :param pfilename: peptide数据文件名
        """
        self.gfilename = gfilename
        self.pfilename = pfilename
        self.pefilename = pefilename
        self.zfile=zfile
        self.vfile=vfile
        self.z_d = {}
        self.v_d = {}
        self.gff_d = {}
        self.pep_d = {}
        self.pep_express_d = {}
        self.model = model

    def load_petide_expression_data(self):
        if self.model != 'Test':
            df = pd.read_excel(self.pefilename)
            columes = list(df)
            for i in tqdm(range(len(df)), desc=f'解析{self.pefilename}文件中'):
                li = []
                for colume in columes:
                    li.append(df.iloc[i][colume])
                d = {
                    'Tas': (li[1] + li[2] + li[3]) / 3,
                    'GP': (li[4] + li[5] + li[6]) / 3,
                    'SiU': (li[7] + li[8] + li[9]) / 3,
                    'EP': (li[10] + li[11] + li[12]) / 3,
                    'FS': (li[13] + li[14] + li[15]) / 3,
                    'JLB3': (li[16] + li[17] + li[18]) / 3,
                    'ML8': (li[19] + li[20] + li[21]) / 3,
                    'VM19D': (li[22] + li[23] + li[24]) / 3,
                    'En8DAP': (li[25] + li[26] + li[27]) / 3,
                    'Em20DAP': (li[28] + li[29] + li[30]) / 3,
                    'IND': (li[31] + li[32] + li[33]) / 3,
                    'PR5D': (li[34] + li[35] + li[36]) / 3,
                    'SR7D': (li[37] + li[38] + li[39]) / 3,

                }
                self.pep_express_d[li[0]] = d
            with open('pep_express_d.json', 'w') as f:
                j_d = json.dumps(self.pep_express_d, cls=JsonEncoder)
                f.write(j_d)
        else:
            if os.path.exists('pep_express_d.json'):
                with open('pep_express_d.json', 'r') as f:
                    self.pep_express_d = json.load(f)
            else:
                df = pd.read_excel(self.pefilename)
                columes = list(df)
                for i in tqdm(range(len(df)), desc=f'解析{self.pefilename}文件中'):
                    li = []
                    for colume in columes:
                        li.append(df.iloc[i][colume])
                    d = {
                        'Tas': (li[1]+li[2]+li[3])/3,
                        'GP':(li[4]+li[5]+li[6])/3,
                        'SiU': (li[7] + li[8] + li[9]) / 3,
                        'EP': (li[10] + li[11] + li[12]) / 3,
                        'FS': (li[13] + li[14] + li[15]) / 3,
                        'JLB3': (li[16] + li[17] + li[18]) / 3,
                        'ML8': (li[19] + li[20] + li[21]) / 3,
                        'VM19D': (li[22] + li[23] + li[24]) / 3,
                        'En8DAP': (li[25] + li[26] + li[27]) / 3,
                        'Em20DAP': (li[28] + li[29] + li[30]) / 3,
                        'IND': (li[31] + li[32] + li[33]) / 3,
                        'PR5D': (li[34] + li[35] + li[36]) / 3,
                        'SR7D': (li[37] + li[38] + li[39]) / 3,

                    }
                    self.pep_express_d[li[0]] = d
                with open('pep_express_d.json', 'w') as f:
                    j_d = json.dumps(self.pep_express_d, cls=JsonEncoder)
                    f.write(j_d)

    def load_petide_data(self):
        if self.model != 'Test':
            df = pd.read_excel(self.pfilename)
            columes = list(df)
            for i in tqdm(range(len(df)), desc=f'解析{self.pfilename}文件中'):
                li = []
                for colume in columes:
                    li.append(df.iloc[i][colume])
                if 'Chr' not in str(li[3]):
                    li[3] = 'Chr' + str(li[3])
                li[3] = li[3].strip()
                if li[3] not in self.pep_d:
                    self.pep_d[li[3]] = [li]
                else:
                    self.pep_d[li[3]].append(li)
            with open('pep_d.json', 'w') as f:
                j_d = json.dumps(self.pep_d, cls=JsonEncoder)
                f.write(j_d)
        else:
            if os.path.exists('pep_d.json'):
                with open('pep_d.json', 'r') as f:
                    self.pep_d = json.load(f)
            else:
                df = pd.read_excel(self.pfilename)
                columes = list(df)
                for i in tqdm(range(len(df)), desc=f'解析{self.pfilename}文件中'):
                    li = []
                    for colume in columes:
                        li.append(df.iloc[i][colume])
                    if 'Chr' not in str(li[3]):
                        li[3] = 'Chr'+str(li[3])
                    li[3] = li[3].strip()
                    # if li[3] == 'Chr1':
                    if li[3] not in self.pep_d:
                        self.pep_d[li[3]] = [li]
                    else:
                        self.pep_d[li[3]].append(li)
                with open('pep_d.json', 'w') as f:
                    j_d = json.dumps(self.pep_d, cls=JsonEncoder)
                    f.write(j_d)

    def load_vfile(self):
        with open(self.vfile, 'r') as f:
            lines = f.readlines()
            for i in lines:
                if 'Transcritp' not in i:
                    li = i.split('\t')
                    self.v_d[li[0]] = li[1]


    def load_zfile(self):
        with open(self.zfile, 'r') as f:
            lines = f.readlines()
            for i in lines:
                if 'TF_ID' not in i:
                    li = i.split('\t')
                    if not self.z_d.get(li[2].strip()):
                        self.z_d[li[2].strip()] = [li[1]]
                    else:
                        self.z_d[li[2].strip()].append(li[1])



    def classify_chr(self):
        with open(self.gfilename) as f:
            lines = f.readlines()
            for i in tqdm(lines,desc='解析gff3文件中'):
                if 'wareLab' in i and 'Chr' in i:
                    i = i.strip()
                    li = i.split('\t')
                    self.gff_d[li[0]] = {'start': int(li[3]), 'end': int(li[4]), 'ranges': [[1,self.start,0]]}


    def calc_hot_region(self):
        with open('Chr_length.csv', 'w') as f:
            f.write('Chr, physical_location, mid\n')
            for k,v in self.gff_d.items():
                f.write(f'{k}, {v["end"]},{(v["end"]+v["start"])/2}\n')
                while True:
                    start = v['ranges'][-1][0]+self.step
                    end = v['ranges'][-1][1]+self.step
                    if end > v['end']:
                        v['ranges'].append([start,v['end'],0])
                        break
                    else:
                        v['ranges'].append([start,end,0])

                for pep in self.pep_d[k]:
                    pep_start = pep[4]
                    pep_end = pep[6]
                    pep_middle = (pep_start+pep_end)/2
                    for r in v['ranges']:
                        if r[0]<pep_middle<r[1]:
                            r[2]+=1
        with open('hot_regions.csv', 'w') as f:
            f.write('Chr, position,count\n')
            for k,v in self.gff_d.items():
                for r in v['ranges']:
                    if r[2] >= self.num:
                        p = (r[0]+r[1])/2
                        f.write(f'{k}, {p},{r[2]}\n')


    def hot_region(self):
        self.load_petide_data()
        self.classify_chr()
        self.calc_hot_region()

    def heatmapA(self):
        self.load_vfile()
        self.load_zfile()
        self.load_petide_data()
        self.load_petide_expression_data()
        with open('cp_tf.csv', 'w') as f1:
            s1 = ','.join(LIST)
            f1.write('Peptide, pep_gene, pep_tf, pep_f,'+s1+'\n')
            for chr, chr_li in self.pep_d.items():
                for pep in chr_li:
                    pep_name = pep[0]
                    pep_gene = pep[2]
                    tf_id = self.v_d.get(pep_gene, 'Na')
                    for f,tfs in self.z_d.items():
                        for tf in tfs:
                            if tf == tf_id:
                                pep_f = f
                                pep_press_d = self.pep_express_d[pep_name]
                                s = f'{pep_name},{pep_gene},{tf_id}, {pep_f}'
                                for i in LIST:
                                    v = pep_press_d[i]
                                    s += f',{v}'
                                f1.write(s+'\n')

    def heatmapBC(self):
        d = {}
        d1 = {}
        line_num = 0
        with open('cp_tf.csv', 'r') as f:
            lines = f.readlines()
            for i in lines:
                if 'Tas' not in i:
                    line_num+=1
                    li = i.strip().split(',')
                    if not d.get(li[3]):
                        d[li[3]]={}
                    li1 = li[4:]
                    for num in range(len(LIST)):
                        if not d[li[3]].get(LIST[num]):
                            d[li[3]][LIST[num]] = [1,0,float(li1[num])]
                        else:
                            d[li[3]][LIST[num]][0] += 1
                            d[li[3]][LIST[num]][2] += float(li1[num])
                        if float(li1[num]) != 0:
                            d[li[3]][LIST[num]][1] += 1
            d2 = {}
            for i in LIST:
                d2[i] = line_num
            for i in lines:
                if 'Tas' not in i:
                    li = i.strip().split(',')
                    li1 = li[4:]
                    for num in range(len(LIST)):
                        if float(li1[num]) == 0:
                            d2[LIST[num]] -= 1
            for f,v in d.items():
                for f1,v1 in v.items():
                    # 平均数
                    avg = v1[2]/v1[0]
                    # tas去掉0的rav数量
                    num0 = v1[1]
                    # 玉米中rav的总数量
                    num1 = v1[0]
                    # 所有Tas去掉0的数量
                    num2 = d2[f1]
                    # 所有CP数量
                    num3 = line_num
                    data = [[num0, num1], [num2, num3]]
                    p = fisher_exact(data)[1]
                    if p >= 0.05:
                        p = 1
                    p = -math.log(p,10)
                    if not d1.get(f1):
                        d1[f1]={f:{'avg':avg, 'p': p}}
                    else:
                        d1[f1][f] = {'avg': avg, 'p':p}
            s1 = ''
            s = ''
            s2 = ''
            with open('cp_tf_avg.csv', 'w') as f:
                d2 = {}
                for i in d.keys():
                    d2[i]=''
                for i in LIST:
                    s1+=i+','
                f.write('label,'+s1+'\n')
                for num in range(len(LIST)):
                    for key in d1[LIST[num]]:
                        if d2[key] == '':
                            d2[key] += f'{key}, {d1[LIST[num]][key]["avg"]}'
                        else:
                            d2[key] += f',{d1[LIST[num]][key]["avg"]}'
                for v in d2.values():
                    f.write(f'{v.strip()}\n')

            s1 = ''
            with open('cp_tf_p.csv', 'w') as f:
                d2 = {}
                for i in d.keys():
                    d2[i] = ''
                for i in LIST:
                    s1 += i + ','
                f.write('label,' + s1 + '\n')
                for num in range(len(LIST)):
                    for key in d1[LIST[num]]:
                        if d2[key] == '':
                            d2[key] += f'{key}, {d1[LIST[num]][key]["p"]}'
                        else:
                            d2[key] += f',{d1[LIST[num]][key]["p"]}'
                for v in d2.values():
                    f.write(f'{v.strip()}\n')












if __name__ == '__main__':
    LIST = ['Tas', 'GP', 'SiU', 'EP', 'FS', 'JLB3', 'ML8', 'VM19D', 'En8DAP', 'Em20DAP', 'IND', 'PR5D', 'SR7D']
    c = Classify('Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3', 'pair-source-protein-peptide.xlsx',
                 'peptide-expression.xlsx',zfile='Zma_TF_list.txt', vfile='V4V3ANNO0928.txt', model='Live')
    c.heatmapA()
    c.heatmapBC()

