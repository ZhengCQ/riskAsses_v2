#!/usr/bin/env  python
# -*- coding:UTF-8 -*-

import sys
import re
import random
import multiprocessing as mp
#sys.path.append("../reohit2/bin/lib/")
import Read as Read
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
		#将存在于vcf中的样本名称和idx存于字典
		samples_info = Read.Readvcf(self.invcf).samples
		for idx, val in enumerate(samples_info.split('|')):
			if val in grp_lst:
				grp_dict.setdefault(grp_name, {}).setdefault('name', []).append(val)
				grp_dict.setdefault(grp_name, {}).setdefault('index', []).append(idx)

		#检查样本在vcf是否存在
		for each in grp_lst:
			try:
				if each not in grp_dict[grp_name]['name']:
					print 'Warnings: %s not in vcf'%(each)
			except:
				print 'Warnings: %s haven\'t sample in vcf'%(grp_name)

	def get_grp_info(self, grp_dict):
		self.handle_grp(grp_dict, self.A, self.A_samples)
		self.handle_grp(grp_dict, self.B, self.B_samples)
		self.handle_grp(grp_dict, self.C, self.C_samples)
		#self.handle_grp(grp_dict, 'GCT', ['GCT_AB0002974', 'GCT_AB0002977', 'GCT_AB0002978', 'GCT_4libs'])
		#self.handle_grp(grp_dict, 'CCT', ['SYSb6745', 'SYSb6746', 'SYSb6747', 'CCT_4libs'])


	def _gt_count(self, gtinfo, grp):
		"""
		  count each mut genotype for grp and each amples
		"""
		def add_num(grp_name, sample_name, mut_type):
			#初始化每个样本及每个突变基因型的值，
			for i in ['ref_hom', 'mut_hom', 'mut_mut', 'mut_het']:
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
			if gt[0] == '0' and gt[1] == '0':#字符为0,非数字
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
		total = mut + ref
		return mut, ref, total

	def functionnal_count(self, gt_mut_dict, A_mut, func, A):
		#计算每组以及每个样本的功能累计
		if A_mut > 0:
			try:
				self.func_count_dict[func][A] += 1 #每个功能值初始化为0				
			except:
				#for k in rank_score:
				self.func_count_dict.setdefault(func, {}).setdefault(A, 0) #每个功能值初始化为0
		try:
			for sample in self.grp_dict[A]['name']:
				mut, ref, allele = self._get_allele_count(gt_mut_dict, sample)
				self.functionnal_count(gt_mut_dict, mut, func, sample)
		except:
			pass

	def _derived_allele_freq(self, gtinfo, func, A, B, C): 
		"""
		gtinfo: 读取到的基因型数据 0/0|0/1|1/1
		grp: 组别信息
		grp_dict： 全局变量的组别字典

		each site i we write the observed derived allele frequency in population A as fAi = dAi / nAi , 
		where nAi is the total number of alleles called there in population A 
		and dAi is the number of derived alleles called
		"""
		A_gt_mut_dict = self._gt_count(gtinfo, A)
		A_mut, A_ref, A_allele = self._get_allele_count(A_gt_mut_dict, A)
		self.functionnal_count(A_gt_mut_dict, A_mut, func, A)

		B_gt_mut_dict = self._gt_count(gtinfo, B)
		B_mut, B_ref, B_allele = self._get_allele_count(B_gt_mut_dict, B)
		self.functionnal_count(B_gt_mut_dict, B_mut, func, B)

		C_gt_mut_dict = self._gt_count(gtinfo, C)
		C_mut, C_ref, C_allele = self._get_allele_count(C_gt_mut_dict, C)
		self.functionnal_count(C_gt_mut_dict, C_mut, func, C)
		
		ni = A_allele
		di = 0
		if A_mut > 0 and B_mut == 0 and C_mut == 0: # 这里ref为祖先碱基。B 群体均为ref，C群体均为ref, A群体携带有alt
			di = A_mut
		elif A_ref > 0 and B_ref == 0  and C_ref == 0: #这里alt为祖先碱基。B群体均为alt，C群体均为alt，A群体携带有ref
			di = A_ref #
		fi = float(di)/float(ni)
		#print di,ni,
		if fi < 1 and fi > 0:
			#print A_grp_gt_list,B_grp_gt_list,C_grp_gt_list,fi,di,ni
			return fi
		else:
			return 0

	def _risk_count(self, freq_dict, category):
		"""
		RA/B(C) = LA,B(C) / LB,A(C)
		"""
		try:
			L_AB_C = sum(freq_dict['AvsB'][category])/sum(freq_dict['AvsB']['intergenic'])
			L_BA_C = sum(freq_dict['BvsA'][category])/sum(freq_dict['BvsA']['intergenic'])
			#print category,sum(freq_dict['AvsB'][category]), sum(freq_dict['AvsB']['intergenic']),sum(freq_dict['BvsA'][category]), sum(freq_dict['BvsA']['intergenic'])
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
			fA = self._derived_allele_freq(gtinfo, func, A, B, C)
			#print 'fA',category,fA,
			fB = self._derived_allele_freq(gtinfo, func, B, A, C)
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

	def format_sum_stat(self):
		func_dict_new = {}
		func_dict_new.setdefault('header', [])
		flag = True
		#对功能根据危害度分值排序
		sort_func_dict = sorted(self.rank_score.items(),key = lambda x:x[1],reverse = True)

		#整理成行为组别及样本对应的数据，列为所有的功能
		for func in [ i[0] for i in sort_func_dict ]:
			for grp in self.grp_dict:
				if flag:
					func_dict_new['header'].append(grp)
				try:
					func_dict_new.setdefault(func, []).append(self.func_count_dict[func][grp])
				except:
					func_dict_new.setdefault(func, []).append(0)
				for sample in self.grp_dict[grp]['name']:
					if flag:
						func_dict_new['header'].append(sample)
					try:
						func_dict_new.setdefault(func, []).append(self.func_count_dict[func][sample])
					except:
						func_dict_new.setdefault(func, []).append(0)
			flag = False		

		funclist = list(func_dict_new.keys()) #后续字典会持续增加sum，dn_ds等key,因此需要先把功能相关的key存为list
		
		for idx, val in enumerate(func_dict_new['header']):
			sumlist = []
			loflist = []
			for k in funclist:
				if k == 'header': continue
				sumlist.append(func_dict_new[k][idx])
				if k in  ['stop_gained', 'stop_lost', 'start_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant']:
					loflist.append(func_dict_new[k][idx])
			total = sum(sumlist)
			lof = sum(loflist)
			dn_ds = float(func_dict_new['missense_variant'][idx])/float(func_dict_new['synonymous_variant'][idx])
			lof_ds = float(lof)/float(func_dict_new['synonymous_variant'][idx])			
			func_dict_new.setdefault('total',[]).append(total)
			func_dict_new.setdefault('dn_ds',[]).append(dn_ds)
			func_dict_new.setdefault('lof_ds',[]).append(lof_ds)

		print 'functionnal' + '\t' + '\t'.join(func_dict_new['header'])
		for k in func_dict_new:
			if k not in self.func_count_dict: continue
			print k,'\t'.join([ str(i) for i in func_dict_new[k]])
		for k in [ 'dn_ds', 'lof_ds', 'total' ]:
			print k,'\t'.join([ str(i) for i in func_dict_new[k]])


	def count(self):
		self.rank_score = rank_score
		self.func_count_dict = {}
		freq_dict = {}
		dnds_dict = {}
		start, end = self._jackknifes(self.fix_sites, self.all_sites) # block jackknifes on the set of sites

		vcfinfo = Read.Readvcf(self.invcf).extract #读取到注释vcf的信息
		num = 0
		for each in vcfinfo:
			num +=1
			if int(num) < int(start) or int(num) > int(end): continue
			high_var = HighVar(each.ann) #遇到多个功能时，取危害度最高的
			#计算不同catergory的di/dn频率值,并统计每一个功能的在群体及个体间数目
			self._freq_count_catalogue(freq_dict, each.gt, high_var.func, self.A, self.B, self.C)
        
        #计算missense, synoymouns及lof的值
		self.sites = end - start
		self.format_sum_stat()
        #群体A和B的dn/ds值
		self.dn_ds_A = float(self.func_count_dict['missense_variant'][self.A])/float(self.func_count_dict['synonymous_variant'][self.A])
		self.dn_ds_B = float(self.func_count_dict['missense_variant'][self.B])/float(self.func_count_dict['synonymous_variant'][self.B])	
		#计算A群体和B群体相对风险值
		self.risk_missense = self._risk_count(freq_dict, 'missense')
		self.risk_synonymous = self._risk_count(freq_dict, 'synonymous')
		self.risk_lof = self._risk_count(freq_dict, 'lof')
		return
