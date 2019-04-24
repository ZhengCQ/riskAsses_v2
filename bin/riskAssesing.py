#!/usr/bin/env  python
# -*- coding:UTF-8 -*-

import sys
import re
import random
#sys.path.append("../reohit2/bin/lib/")
import Read as Read

rank_score = {
	'stop_gained': 10,
	'stop_retained_variant': 10, 
	'stop_lost': 10,
	'start_lost': 10,
	'splice_acceptor_variant': 9,
	'splice_donor_variant': 9,
	'frameshift_variant': 9,
	'missense_variant': 8,
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


class FuncRisk(object):
	"""count functionnal，获得每一个功能的数目"""
	def __init__(self, invcf, A, B,
		         A_samples, B_samples, C_samples,
		         fix_sites, all_sites):
		"""
		invcf: snpEff注释后的vcf
		A: A population name
		B: B population name
		A_samples: A samples 数组
		B_samples: B samples 数组
		"""
		self.invcf = invcf
		self.A = A
		self.B = B
		self.C = 'C'
		self.A_samples = A_samples
		self.B_samples = B_samples
		self.C_samples = C_samples
		self.fix_sites = fix_sites
		self.all_sites = all_sites	
		self.grp_dict = {}
		""" #根据组别信息，得到组别相关的字典, 
		grp_dict[grp]['index'] = 【2，3，4]
		grp_dict[grp]['name'] = 【'SYSb6745', 'SYSb6746', 'SYSb6747']
		"""
		self.get_grp_info(self.grp_dict)
		#A = 'CCT'
		#B = 'GCT'
		self.count()

	def handle_grp(self, grp_dict, grp_name, grp_lst):
		grp_dict.setdefault(grp_name, {}).setdefault('name', grp_lst)
		samples_info = Read.Readvcf(self.invcf).samples
		for idx, val in enumerate(samples_info.split('|')):
			if val in grp_dict[grp_name]['name']:
				grp_dict[grp_name].setdefault('index', []).append(idx)

	def get_grp_info(self, grp_dict):
		self.handle_grp(grp_dict, self.A, self.A_samples)
		self.handle_grp(grp_dict, self.B, self.B_samples)
		self.handle_grp(grp_dict, self.C, self.C_samples)
		#self.handle_grp(grp_dict, 'GCT', ['GCT_AB0002974', 'GCT_AB0002977', 'GCT_AB0002978', 'GCT_4libs'])
		#self.handle_grp(grp_dict, 'CCT', ['SYSb6745', 'SYSb6746', 'SYSb6747', 'CCT_4libs'])


	def _allele_count(self, gtinfo, grp):
		grp_gt_list = [gtinfo.split('|')[i] for i in self.grp_dict[grp]['index']]
		mut = 0
		for each in grp_gt_list:
			for allele in each.split('/'):
				if int(allele) != 0: # 不等于0则计算，1，2都算
					mut += 1
		#print grp_gt_list, mut, len(grp_gt_list)*2,
		return mut, len(grp_gt_list)*2

	def _derived_allele_freq(self, gtinfo, A, B, C): 
		"""
		gtinfo: 读取到的基因型数据 0/0|0/1|1/1
		grp: 组别信息
		grp_dict： 全局变量的组别字典

		each site i we write the observed derived allele frequency in population A as fAi = dAi / nAi , 
		where nAi is the total number of alleles called there in population A 
		and dAi is the number of derived alleles called
		"""
		A_mut, A_allele = self._allele_count(gtinfo, A)
		B_mut, B_allele = self._allele_count(gtinfo, B)
		C_mut, C_allele = self._allele_count(gtinfo, C)
		ni = A_allele
		di = 0
		if A_mut > 0 and B_mut == 0 and C_mut == 0: # 这里ref为祖先碱基。B 群体均为ref，C群体均为ref, A群体携带有alt
			di = A_mut
		elif A_mut < A_allele and B_mut == B_allele and C_mut == C_allele: #这里alt为祖先碱基。B群体均为alt，C群体均为alt，A群体携带有ref
			di = A_allele - A_mut #
		fi = float(di)/float(ni)
		#print di,ni,
		return fi

	def _risk_count(self, freq_dict, category):
		"""
		RA/B(C) = LA,B(C) / LB,A(C)
		"""
		try:
			L_AB_C = sum(freq_dict['AvsB'][category])/sum(freq_dict['AvsB']['intergenic'])
			L_BA_C = sum(freq_dict['BvsA'][category])/sum(freq_dict['BvsA']['intergenic'])
			print category,sum(freq_dict['AvsB'][category]), sum(freq_dict['AvsB']['intergenic']),sum(freq_dict['BvsA'][category]), sum(freq_dict['BvsA']['intergenic'])
		except:
			L_AB_C = 0
			L_BA_C = 0
		try:
			return L_AB_C/L_BA_C
		except:
			return 0		

	def _jackknifes(self, fix_sites, all_sites):
		start = random.randint(0, all_sites) # get start from 0 and all_sites
		#end = random.randint(start, all_sites) # get end from start and all_sites
		end = start + fix_sites
		if end > all_sites:
			start, end = self._jackknifes(fix_sites, all_sites) #重新得到的start, end
		return start, end

	def dn_ds_count(self, dnds_dict, gtinfo, func, A):
		A_mut, A_allele = self._allele_count(gtinfo, A)
		if A_mut > 0:
			if func == 'missense_variant':
				try:
					dnds_dict[A]['dn'] = dnds_dict[A]['dn'] + 1
				except:
					dnds_dict.setdefault(A, {}).setdefault('dn', 1)
			if func == 'synonymous_variant':
				try:
					dnds_dict[A]['ds'] = dnds_dict[A]['ds'] + 1
				except:
					dnds_dict.setdefault(A, {}).setdefault('ds', 1)

	def _freq_count_catalogue(self, freq_dict, gtinfo, func, A, B, C):	
		def freq_count(category):
			"""
			分别计算A population 和 B population derived alleles frequence
			分别添加AvsB,和BvsA的，不同categrory的数组
			freq_dict: 字典
			gtinfo：基因型数据
			A: 群体A
			B: 群体B
			category：不同类型，如missse, lof, intergenic等
			"""
			fA = self._derived_allele_freq(gtinfo, A, B, C)
			#print 'fA',category,fA
			fB = self._derived_allele_freq(gtinfo, B, A, C)
			#print 'fB',category,fB
			freq_dict.setdefault('AvsB', {}).setdefault(category, []).append(fA*(1 - fB))
			freq_dict.setdefault('BvsA', {}).setdefault(category, []).append(fB*(1 - fA))

        #计算不同catergory的di/dn频率值
		if func == 'missense_variant':
			#计算A vs B 及B vs A的，每一个错义突变位点的频率计算总和
			freq_count('missense')
		if func == 'synonymous_variant':
			#计算A vs B 及B vs A的，每一个同义突变位点的频率计算总和
			freq_count('synonymous')
		if func in  ['stop_gained', 'stop_lost', 'start_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant']:
			"""计算A vs B 及B vs A的，每一个lof位点的频率计算总和
			'stop_gained': 10,  
			'stop_lost': 10,
			'start_lost': 10,
			'splice_acceptor_variant': 9,
			'splice_donor_variant': 9,
			'frameshift_variant': 9,
			"""
			freq_count('lof')
		if func == 'intergenic_region':
     		#计算A vs B 及B vs A的，每一个基因减去突变位点的频率计算总和
			freq_count('intergenic')
		
	def count(self):
		self.rank_score = rank_score
		self.func_count_hash = {}
		freq_dict = {}
		dnds_dict = {}
		start, end = self._jackknifes(self.fix_sites, self.all_sites) # block jackknifes on the set of sites
		#print start, end
		for k in rank_score:
			self.func_count_hash.setdefault(k, 0) #每个功能值初始化为0
		vcfinfo = Read.Readvcf(self.invcf).extract #读取到注释vcf的信息
		num = 0
		for each in vcfinfo: #遍历每一行数据
			num +=1
			if int(num) < int(start) or int(num) > int(end): continue
			
			high_var = HighVar(each.ann) #遇到多个功能时，取危害度最高的
			
			#计算不同catergory的di/dn频率值
			self._freq_count_catalogue(freq_dict, each.gt, high_var.func, self.A, self.B, self.C)
			#分别计算群体A和群体B的dnds值
			self.dn_ds_count(dnds_dict, each.gt, high_var.func, self.A)
			self.dn_ds_count(dnds_dict, each.gt, high_var.func, self.B)
			#self.func_count_hash[high_var.func] += 1 #计算每一个功能
        
        #计算missense, synoymouns及lof的值
		self.sites = end - start

        #群体A和B的dn/ds值
		self.dn_ds_A = float(dnds_dict[self.A]['dn'])/float(dnds_dict[self.A]['ds'])
		self.dn_ds_B = float(dnds_dict[self.B]['dn'])/float(dnds_dict[self.B]['ds'])	
		#计算A群体和B群体相对风险值
		self.risk_missense = self._risk_count(freq_dict, 'missense')
		self.risk_synonymous = self._risk_count(freq_dict, 'synonymous')
		self.risk_lof = self._risk_count(freq_dict, 'lof')
		return
