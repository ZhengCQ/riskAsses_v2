#!/usr/bin/env  python
# -*- coding:UTF-8 -*-
# @Author: Zheng ChenQing
# @Date: 2021.10.18
# @E-mail: zhengchenqing@qq.com

import sys
import re
import random
import json
import multiprocessing as mp
#sys.path.append("../reohit2/bin/lib/")
import bin.Read as Read
import datetime
rank_score = {
	'stop_gained': 10,
	'stop_retained_variant': 10, 
	'stop_lost': 10,
	'start_lost': 10,
	'splice_acceptor_variant': 9,
	'splice_donor_variant': 9,
	'frameshift_variant': 9,
	'missense_variant': 8,
	'missense_deleterious': 9,
	'missense_benign': 7,
	'disruptive_inframe_deletion': 4,
	'splice_region_variant': 4,
	'intron_variant': 3,
	'3_prime_UTR_variant': 3,      
	'inframe_deletion': 3,
	'inframe_insertion': 3, 
	'5_prime_UTR_premature_start_codon_gain_variant': 2.5,
	'5_prime_UTR_variant': 2.5,
	'upstream_gene_variant': 2.5,
	'downstream_gene_variant': 2.5,
	'TF_binding_site_variant': 2,
	'non_coding_exon_variant': 2,
	'synonymous_variant': 1.5,
	'intergenic_region': 1,
	'intragenic_variant': 1.5,
	'sequence_feature': 1,
}

Gscores = json.load(open('/BIGDATA1/sysuls_yliu_1/linhongzhou/biosoft/riskAsses_v2-master/db/Grantham_Scores.json'))

class HighVar(object):
	"""docstring for HighVariant"""
	def __init__(self, anno):
		self.anno = anno
		self.high_variant(self.anno.functions)

	def high_variant(self, funcLst):
		#值初始化
		self.func = ''
		self.hgvs_c = ''
		self.hgvs_p = ''
		self.gene = ''
		self.transcript = ''
		index = ''
		max_func = self.func
		#print(funcLst)
		for idx, val in enumerate(funcLst):
			if re.search(r'&',val): # splice_region_variant&intron_variant 存在 &，且不为空值
				tmpLst = val.split('&') #存在一个& 及多个&的情况
				for i in range(len(tmpLst) - 1):
					func = self.high_func(tmpLst[i], tmpLst[i + 1])
					max_func = self.high_func(self.func, func)
			else: # func不为空值的情况
				max_func = self.high_func(self.func, val)
			#获得的最大危害功能与原来不等时，重新赋值
			if self.func != max_func:
				self.func = max_func
				index = idx	
				#print(index)
				#print(self.func)
		# 根据危害度最大的功能的值，获取相应的其他值
		self.hgvs_c = self.anno.hgvs_cs[index]
		self.hgvs_p = self.anno.hgvs_ps[index]
		self.gene = self.anno.genes[index]
		self.transcript = self.anno.transcripts[index]
		return

	def high_func(self, a, b): #判断两个功能类型，返回危害度最大的
		self.rank_score = rank_score
		if a and b: # a,b均不为空值
			try: 
				self.rank_score[a] > self.rank_score[b] #存在字典并大于
				return a
			except:
				return b 
		elif a: # b为空值
			return a
		elif b: # a为空值
			return b
		else: # a,b均为空值
			return ''

class HandleGroup(object):
	def __init__(self, invcf, grp_dict, grp_name, grp_lst):
		self.invcf = invcf
		self.grp_dict = grp_dict
		self.grp_name = grp_name
		self.grp_lst = grp_lst
		self.handle_grp()

	def handle_grp(self):
	#将存在于vcf中的样本名称和idx存于字典
		samples_info = Read.Readvcf(self.invcf).samples
		for idx, val in enumerate(samples_info.split('|')):
			if val in self.grp_lst:
				self.grp_dict.setdefault(self.grp_name, {}).setdefault('name', []).append(val)
				self.grp_dict.setdefault(self.grp_name, {}).setdefault('index', []).append(idx)

		#检查样本在vcf是否存在
		for each in self.grp_lst:
			try:
				if each not in self.grp_dict[self.grp_name]['name']:
					print ('Warnings: %s not in vcf'%(each))
			except:
				print ('Warnings: %s haven\'t sample in vcf'%(self.grp_name))


class FuncRisk(object):
	"""count functional，获得每一个功能的数目"""
	def __init__(self, i, invcf, A, B,
		         A_samples, B_samples, C_samples,
		         fix_sites, all_sites, work_dir):
		"""
		invcf: snpEff注释后的vcf
		A: A population name
		B: B population name
		A_samples: A samples 数组
		B_samples: B samples 数组
		"""
		self.i = i
		self.invcf = invcf
		self.A = A
		self.B = B
		self.C = 'C_OutGroup'
		self.A_samples = A_samples
		self.B_samples = B_samples
		self.C_samples = C_samples
		self.fix_sites = fix_sites
		self.all_sites = all_sites
		self.work_dir  = work_dir
		self.grp_dict = {}
		""" #根据组别信息，得到组别相关的字典, 
		grp_dict[grp]['index'] = 【2，3，4]
		grp_dict[grp]['name'] = 【'SYSb6745', 'SYSb6746', 'SYSb6747']
		"""
		self.get_grp_info()
		#查看是否正确获取样本名和idx
		self.main()

	def get_grp_info(self):
		"""
		@return grp_dict: {'GCT': ['GCT_AB0002974', 'GCT_AB0002977', 'GCT_AB0002978', 'GCT_4libs']}
		"""
		HandleGroup(self.invcf, self.grp_dict, self.A, self.A_samples)
		HandleGroup(self.invcf, self.grp_dict, self.B, self.B_samples)
		HandleGroup(self.invcf, self.grp_dict, self.C, self.C_samples)

	def _gt_count(self, gtinfo, grp):
		"""
		  count each mut genotype for grp and each amples
		"""
		def add_num(grp_name, sample_name, mut_type):
			#初始化每个样本及每个突变基因型的值，
			for i in ['ref_hom', 'mut_hom', 'mut_mut', 'mut_het', 'miss']:
				gt_mut_dict.setdefault(sample_name, {}).setdefault(i, 0)
				try:
					gt_mut_dict[grp_name][i]
				except:
					gt_mut_dict.setdefault(grp_name, {}).setdefault(i, 0)
			#指定类型的值增加
			gt_mut_dict[grp_name][mut_type] += 1
			gt_mut_dict[sample_name][mut_type] += 1

		grp_gt_list = [gtinfo.split('|')[i] for i in self.grp_dict[grp]['index']]
		gt_mut_dict = {}

		for idx, val in enumerate(grp_gt_list):
			sample_name = self.grp_dict[grp]['name'][idx]
			gt = val.split('/')
			if gt[0] == '.' and gt[1] == '.':
				add_num(grp, sample_name, 'miss') #./.				
			elif gt[0] == '0' and gt[1] == '0':#字符为0,非数字
				add_num(grp, sample_name, 'ref_hom') #0/0
			elif gt[0] != '0' and gt[1] != '0':
				if gt[0] == gt[1]:
					add_num(grp, sample_name, 'mut_hom') #1/1
				else:
					add_num(grp, sample_name, 'mut_mut') #类似突变为1/2
			elif gt[0] != '0' or gt[1] != '0':
				add_num(grp, sample_name, 'mut_het') #0/1，1/0, 2/0等
		return gt_mut_dict

	def _get_allele_count(self, gt_mut_dict, name):
		mut = gt_mut_dict[name]['mut_het']*1 + gt_mut_dict[name]['mut_mut']*2 +  gt_mut_dict[name]['mut_hom']*2
		ref = gt_mut_dict[name]['mut_het']*1 + gt_mut_dict[name]['ref_hom']*2
		return mut, ref

	def functional_count(self, gt_mut_dict, A_mut, func, A):
		#计算每组以及每个样本的功能累计
		if A_mut > 0: #每个组的突变信息在外层获取，是因为前面有3个组的mut要先获取，用来比较。避免重复获取
			try:
				self.func_count_dict[func][A] += 1 			
			except:
				#for k in rank_score:
				self.func_count_dict.setdefault(func, {}).setdefault(A, 0) #每个功能值初始化为0
		try:
			#计算每个样本的功能累计
			for sample in self.grp_dict[A]['name']:
				mut, ref = self._get_allele_count(gt_mut_dict, sample) #gt_mut_dict用来获取样本的突变信息
				self.functional_count(gt_mut_dict, mut, func, sample)
		except:
			pass
			#print 'functional_count error for sample count'

	def functional_count_fi(self, A_fi_mut_dict, B_fi_mut_dict, func, A, B):
		"""
		A_fi_mut_dict: {'A':val1, 'sample1':val2}
		B_fi_mut_dict: {'B':val1, 'sample1':val2}
		A: 组
		B: 组
		两个组都要计算，derived相对值需要A和B一起计算。
		"""
		def add_dict(func, A, fA):
			try:
				self.func_di_count_dict[func][A] += fA
			except:
				self.func_di_count_dict.setdefault(func, {}).setdefault(A, fA)

		fA = A_fi_mut_dict[A]['di_alf']
		fB = B_fi_mut_dict[B]['di_alf']
		add_dict(func, 'AvsB', fA*(1 - fB))
		add_dict(func, 'BvsA', fB*(1 - fA))
		#self.func_di_count_dict.setdefault(func, {}).setdefault('AvsB', []).append(fA*(1 - fB))
		#self.func_di_count_dict.setdefault(func, {}).setdefault('BvsA', []).append(fB*(1 - fA))
		add_dict(func, A, fA)
		add_dict(func, B, fB)
		for sample in self.grp_dict[A]['name']:
			add_dict(func, sample, A_fi_mut_dict[sample]['di_alf'])
		for sample in self.grp_dict[B]['name']:
			add_dict(func, sample, B_fi_mut_dict[sample]['di_alf'])

	def functional_count_homhet(self, A_fi_mut_dict, B_fi_mut_dict, func, A, B):
		"""
		{'A': {'di_alf': 0.75, 'mut_hom': 2, 'mut_het': 3}}
		"""
		def add_dict(homhet_dict, A_fi_mut_dict, func, A, mut_type):
			try:
				homhet_dict[mut_type][func][A] +=  A_fi_mut_dict[A][mut_type]
			except:
				homhet_dict.setdefault(mut_type,{})\
					       .setdefault(func, {}).setdefault(A, A_fi_mut_dict[A][mut_type])

		for mut_dict, grp in zip([A_fi_mut_dict, B_fi_mut_dict],[A, B]):
			### add all mut hom info
			for mut_type in ['mut_hom', 'mut_het', 'ref_hom']:
				add_dict(self.func_homhet_count_dict, mut_dict, func, grp, mut_type)
				## 每个样本的
				for sample in self.grp_dict[grp]['name']:
					add_dict(self.func_homhet_count_dict, mut_dict, func, sample, mut_type)
			### add derived mut and hom info
			if mut_dict[grp]['di_type'] != 'NONE':
				### 添加di het 
				add_dict(self.func_di_homhet_count_dict, mut_dict, func, grp, 'mut_het')
				for sample in self.grp_dict[grp]['name']:
					add_dict(self.func_di_homhet_count_dict, mut_dict, func, sample, 'mut_het')
				### 添加di hom
				if mut_dict[grp]['di_type'] == 'ALT':
					add_dict(self.func_di_homhet_count_dict, mut_dict, func, grp, 'mut_hom')
					for sample in self.grp_dict[grp]['name']:
						add_dict(self.func_di_homhet_count_dict, mut_dict, func, sample, 'mut_hom')
				elif mut_dict[grp]['di_type'] == 'REF':
					add_dict(self.func_di_homhet_count_dict, mut_dict, func, grp, 'ref_hom')
					for sample in self.grp_dict[grp]['name']:
						add_dict(self.func_di_homhet_count_dict, mut_dict, func, sample, 'ref_hom')

	def _mut_fi_count(self, gtinfo, func, A, B, C): 
		"""
		gtinfo: 读取到的基因型数据 0/0|0/1|1/1
		grp: 组别信息
		grp_dict： 全局变量的组别字典

		each site i we write the observed derived allele frequency in population A as fAi = dAi / nAi , 
		where nAi is the total number of alleles called there in population A 
		and dAi is the number of derived alleles called
		"""
		A_gt_mut_dict = self._gt_count(gtinfo, A)
		A_mut, A_ref = self._get_allele_count(A_gt_mut_dict, A)
		self.functional_count(A_gt_mut_dict, A_mut, func, A)

		B_gt_mut_dict = self._gt_count(gtinfo, B)
		B_mut, B_ref = self._get_allele_count(B_gt_mut_dict, B)
		self.functional_count(B_gt_mut_dict, B_mut, func, B)

		C_gt_mut_dict = self._gt_count(gtinfo, C)
		C_mut, C_ref = self._get_allele_count(C_gt_mut_dict, C)
		#print C_mut, C_ref
		self.functional_count(C_gt_mut_dict, C_mut, func, C)

		A_fi_mut_dict = {}
		B_fi_mut_dict = {}
		self._fi_count(A_fi_mut_dict, A_gt_mut_dict, A_mut, A_ref, B_mut, B_ref, C_mut, C_ref, A)
		self._fi_count(B_fi_mut_dict, B_gt_mut_dict, B_mut, B_ref, A_mut, A_ref, C_mut, C_ref, B)
        
		return (A_fi_mut_dict, B_fi_mut_dict)
	
	def _fi_count(self, fi_mut_dict, gt_mut_dict, A_mut, A_ref, B_mut, B_ref, C_mut, C_ref, A):
		#计算 每个组及每个样本的 derived_allele_freq，数目总和，并存在gt_mut_dict中，
		fA, di_type = self._derived_allele_freq(A_mut, A_ref, B_mut, B_ref, C_mut, C_ref)
		fi_mut_dict.setdefault(A, {}).setdefault('di_alf', fA)
		fi_mut_dict.setdefault(A, {}).setdefault('di_type', di_type)
		##突变基因型的数目
		fi_mut_dict.setdefault(A, {}).setdefault('mut_hom', gt_mut_dict[A]['mut_hom'])
		fi_mut_dict.setdefault(A, {}).setdefault('mut_het', gt_mut_dict[A]['mut_het'])
		fi_mut_dict.setdefault(A, {}).setdefault('ref_hom', gt_mut_dict[A]['ref_hom'])
		
		if A in self.grp_dict:
			## 增加每一样本的di_alf信息
			for sample in self.grp_dict[A]['name']:
				mut, ref = self._get_allele_count(gt_mut_dict, sample) #gt_mut_dict用来获取样本的突变信息
				self._fi_count(fi_mut_dict, gt_mut_dict, mut, ref, B_mut, B_ref, C_mut, C_ref, sample)

	def _derived_allele_freq(self, A_mut, A_ref, B_mut, B_ref, C_mut, C_ref, freq_cut = 1): #A群体，B群体及外群C
		#derived allele 的计算公式
		di_type = 'NONE'
		ni = A_mut + A_ref
		di = 0
		if A_mut > 0 and B_mut == 0 and C_mut == 0: # 这里ref为祖先碱基。B 群体均为ref，C群体均为ref, A群体携带有alt
			di = A_mut
			di_type = 'ALT'
		elif A_ref > 0 and B_ref == 0  and C_ref == 0: #这里alt为祖先碱基。B群体均为alt，C群体均为alt，A群体携带有ref
			di = A_ref #
			di_type = 'REF'
		if ni == 0:
			fi = 0
		else:
			fi = float(di)/float(ni)
		#print di,ni,
		if fi <= freq_cut and fi > 0: #
			#print A_grp_gt_list,B_grp_gt_list,C_grp_gt_list,fi,di,ni
			return fi, di_type
		else:
			return 0, di_type
		

	def _risk_count(self, freq_dict, func1, func2):
		"""
		RA/B(C) = LA,B(C) / LB,A(C)
		"""
		#获取AvsB，及BvsA对应的idx位置
		AvsB_idx = ''
		BvsA_idx = ''
		for idx, val in enumerate(freq_dict['functional']):
			if val == 'AvsB':
				AvsB_idx = idx
			if val == 'BvsA':
				BvsA_idx = idx
		try:
			L_AB_C = freq_dict[func1][AvsB_idx]/freq_dict[func2][AvsB_idx]
			L_BA_C = freq_dict[func1][BvsA_idx]/freq_dict[func2][BvsA_idx]
			#print category,sum(freq_dict['AvsB'][category]), sum(freq_dict['AvsB']['intergenic']),sum(freq_dict['BvsA'][category]), sum(freq_dict['BvsA']['intergenic'])
		except:
			L_AB_C = 0
			L_BA_C = 0
		try:
			return float(L_AB_C)/float(L_BA_C)
		except:
			return 0		

	def _jackknifes(self, fix_sites, all_sites):
		"""
		刀切固定长度
		@param fix_sites: 固定位点数 
		@param all_sites: 总的数目
		"""
		start = random.randint(0, all_sites) # get start from 0 and all_sites
		#end = random.randint(start, all_sites) # get end from start and all_sites
		end = start + fix_sites
		if end > all_sites:
			start, end = self._jackknifes(fix_sites, all_sites) #重新得到的start, end
		return start, end

	def main(self):
		self.rank_score = rank_score
		self.func_count_dict = {}
		self.func_di_count_dict = {}
		self.func_homhet_count_dict = {}
		self.func_di_homhet_count_dict = {}

		start = 0
		end = self.all_sites
		if self.fix_sites < self.all_sites:
			start, end = self._jackknifes(self.fix_sites, self.all_sites) # block jackknifes on the set of sites
		vcfinfo = Read.Readvcf(self.invcf).extract #读取到注释vcf的信息
		samples = Read.Readvcf(self.invcf).samples #读取样本信息
		num = 0
		derived_freq_fi = open('%s/derived_freq_siteinfo_%s.txt'%(self.work_dir, self.i),'w')
		header = 'Chr\tPos\tFunc\tHgv_p\t{0}_di_freq\t{1}_di_freq\t{0}_di_type\t{1}_di_type\tScore\t{0}_gt_refhom\t{1}_gt_refhom\t{0}_gt_muthom\t{1}_gt_muthom\t{0}_gt_muthet\t{1}_gt_muthet'.format(self.A,self.B)
		
		header_lst = header.split('\t') + samples.split('|')
		derived_freq_fi.write("\t".join(header_lst) + '\n')

		for each in vcfinfo:
			num +=1
			if int(num) < int(start) or int(num) > int(end): continue

			high_var = HighVar(each.ann) #遇到多个功能时，取危害度最高的
			#计算不同catergory的di/dn频率值,并统计每一个功能的在群体及个体间数目
			## print gtinfo, high_var.func, high_var.hgvs_p, 'CCT', A_fi_mut_dict['CCT']['di_alf'],'GCT',B_fi_mut_dict['GCT']['di_alf']
			A_fi_mut_dict, B_fi_mut_dict = self._mut_fi_count(each.gt, high_var.func, self.A, self.B, self.C)
			
			## func_di_count_dict 统计di
			self.functional_count_fi(A_fi_mut_dict, B_fi_mut_dict, high_var.func, self.A, self.B)
			## 统计hom 和het数目
			self.functional_count_homhet(A_fi_mut_dict, B_fi_mut_dict, high_var.func, self.A, self.B)

			#### 这里记录错义突变的分值
			missense_score = '.'
			if high_var.func == 'missense_variant':
				misssense_class = ''
				aa1,aa2 = re.search(r'p\.(\w{3})\d+(\w{3})',high_var.hgvs_p).groups()[:] ## p.Arg122Met
				try:
					missense_score = Gscores[aa1][aa2]
				except:
					missense_score = Gscores[aa2][aa1]
				
				if missense_score >= 150:
					misssense_class = 'missense_deleterious'
				else:
					misssense_class = 'missense_benign'
				
				A_fi_mut_dict, B_fi_mut_dict = self._mut_fi_count(each.gt, misssense_class, self.A, self.B, self.C)

				## func_di_count_dict 统计di
				self.functional_count_fi(A_fi_mut_dict, B_fi_mut_dict, misssense_class, self.A, self.B)
				## 统计hom 和het数目
				self.functional_count_homhet(A_fi_mut_dict, B_fi_mut_dict, misssense_class, self.A, self.B)
            
			### 输出每一条记录数据
			results = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".\
				format(each.chrom, each.pos, high_var.func, \
				high_var.hgvs_p, A_fi_mut_dict[self.A]['di_alf'],\
				B_fi_mut_dict[self.B]['di_alf'], A_fi_mut_dict[self.A]['di_type'],\
				B_fi_mut_dict[self.B]['di_type'],missense_score,\
				A_fi_mut_dict[self.A]['ref_hom'],B_fi_mut_dict[self.B]['ref_hom'],\
				A_fi_mut_dict[self.A]['mut_hom'],B_fi_mut_dict[self.B]['mut_hom'],\
				A_fi_mut_dict[self.A]['mut_het'],B_fi_mut_dict[self.B]['mut_het'])
			
			results_lst = results.split('\t') + each.gt.split('|')
			derived_freq_fi.write("\t".join(results_lst) + '\n')

		derived_freq_fi.close()
		return
