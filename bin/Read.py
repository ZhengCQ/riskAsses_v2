#!/usr/bin/env	python
# -*- coding:UTF-8 -*-

#    /\_/\
#  =( °w° )=
#    )   (  //
#   (__ __)//
#  前人挖坑
# Author: Yeqi Zhang
# Date: 2019.03.01

"""
Read vcf file
"""

import sys, re
import pysam
import subprocess
import logging
import gzip
import datetime

class Snpeff(object):
	"""read INFO and return ANN, LOF and NMD."""
	def __init__(self, info):
		"""
		@param info:             INFO column in vcf, include interpretation from SnpEff
		"""
		self.info = info
		try:
			self.ann = re.search('ANN=([^;]+)', self.info).group()
		except:
			self.ann = None
		self.fun()
		self.lof_nmd()
		#print(self.ann)

	def fun(self):
		"""from snpEff 'ANN' string, get function, gene, transcript, rank, HGVS.c, HGVS.p for each variant
		@return:             lists of alt, functions, genes, transcripts, ranks, hgvs_cs, hgvs_ps"""
		self.ann_alts = []  # alternative allele for each transcript
		self.functions = []  # function for each transcript
		self.genes = []  # gene name for each transcript
		self.transcripts = []
		self.ranks = []
		self.hgvs_cs = []
		self.hgvs_ps = []
		if not self.ann:
			return # alts, functions, genes, transcripts, ranks, hgvs_cs, hgvs_ps
		for ann_str in self.ann.split(','):
			anns = ann_str.split('|')
			# 1-5: Allele, Annotation, Annotation_Impact, Gene_Name, Gene_ID,
			# 6-10: Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c,
			# 11-16: HGVS.p, cDNA.pos/cDNA.length, CDS.pos/CDS.length, AA.pos/AA.length, Distance, ERRORS/WARNINGS/INFO
			alt = anns[0] if anns[0] else '.'  # e.g. T
			self.ann_alts.append(alt)
			function = anns[1] if anns[1] else '.'  # e.g. missense_variant
			self.functions.append(function)
			gene = anns[3] if anns[3] else '.'  # e.g. ATAD3C
			self.genes.append(gene)
			transcript = anns[6] if anns[6] else '.'  # e.g. ENST00000378785
			self.transcripts.append(transcript)
			rank = anns[8] if anns[8] else '.'  # e.g. 3/12
			self.ranks.append(rank)
			hgvs_c = anns[9] if anns[9] else '.'  # e.g. c.184C>T
			self.hgvs_cs.append(hgvs_c)
			hgvs_p = anns[10] if anns[10] else '.'  # e.g. p.Arg62Cys
			self.hgvs_ps.append(hgvs_p)
		return # alts, functions, genes, transcripts, ranks, hgvs_cs, hgvs_ps

	def lof_nmd(self):
		""" get Predicted loss of function / nonsense mediated decay
		@return pred_lof:       Predicted loss of function: [Number_of_transcripts_in_gene, Percent_of_transcripts_affected]
									e.g: ['2', '0.5']
		@return pred_nmd:       Predicted nonsense mediated decay: [Number_of_transcripts_in_gene,
									Percent_of_transcripts_affected]. e.g: ['1', '1']
		"""
		self.pred_lof = '.'  # Predicted loss of function: [Number_of_transcripts_in_gene, Percent_of_transcripts_affected]
		self.pred_nmd = '.'  # Predicted nonsense mediated decay: [Number_of_transcripts_in_gene, Percent_of_transcripts_affected]
		search_lof = re.search('LOF=\(([^;]+)\)', self.ann)
		search_nmd = re.search('NMD=\(([^;]+)\)', self.ann)
		if search_lof:  # Gene_Name | Gene_ID | Number_of_transcripts_in_gene |Percent_of_transcripts_affected
			self.pred_lof = search_lof.group(1).split('|')[2:]
		if search_nmd:
			self.pred_nmd = search_nmd.group(1).split('|')[2:]
		return # pred_lof, pred_nmd


class Vcfline(object):
	"analyse vcf line info"
	def __init__(self, line):
		self.line = line
		self.infos = self.line.strip().split('\t')
		self.analyse()

	def analyse(self):
		"""add splice infos to self.
		splice: CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLEs...
		to: self. chrom, pos, rs, ref, alt, qual, filt, info, fmt, samples"""
		#if self.line.startswith('#CH'):
		#	self.samples = self.line.strip().split('\t')[9:]
		#	self.samples = ('|').join(self.samples)
		s = datetime.datetime.now()
		if self.line.startswith('#'):
			logging.error("The line starts with '#'!")
			sys.exit(1)
		self.chrom = self.infos[0]
		self.pos = self.infos[1]
		self.rs = self.infos[2]
		self.ref = self.infos[3]
		self.alt = self.infos[4]
		self.qual = self.infos[5]
		self.filt = self.infos[6]
		self.info = self.infos[7]
		self.fmt = self.infos[8]
		self.detail = self.infos[9:]

		# SnpEff
		self.ann = Snpeff(self.info)

		self.gt, self.ad, self.dp, self.zygous = [], [], [], []  # set
		
		# only one alternative allele. No mutil-alt.
		self.len_var = abs(len(self.alt) - len(self.ref))  # variant length. 0: SNP; >0: indel
		for single in self.detail:
			#print single
			more = single.split(':')
			if re.match('GT:AD:DP', self.fmt):  # here, just consider one sample, need improve to multiple samples
				# GT:AD:DP:GQ:PL  1/1:0,2155:2155:99:85781,6535,0
				self.gt.append(more[0])
				try:
				    self.ad.append(more[1].split(',')[1])
				except:
					self.ad.append('.')
				self.dp.append(more[2])  # alternative allelic depth
			else:
				self.gt.append(more[0])  # GT:GQ:PL
				self.ad.append('.')
				self.dp.append('.')

			if re.match('\d', more[0]):
				gts = re.split('[|/]', more[0])
				if gts[0] == gts[1]:
					self.zygous.append('homozygous')
				else:
					self.zygous.append('heterozygous')
		self.gt = ('|').join(self.gt)
		self.ad = ('|').join(self.ad)
		self.dp = ('|').join(self.dp)
		self.zygous = ('|').join(self.zygous)
		e = datetime.datetime.now()
		#logging.error("Read vcf line time: %s" % (e-s))
		return

class Readvcf(object):
	"read vcf file"
	def __init__(self, vcf):
		self.vcf = vcf
		self.get_samples()
		self.extract = self.extract()
		#super(Vcfline)

	def get_samples(self):
		"""get samples from #CHROM line"""
		with gzip.open(self.vcf, 'r') as infile:
			for line in infile:
				line = line.decode()
				if line.startswith("#C"):
					self.samples = line.strip().split('\t')[9:]
					self.samples = ('|').join(self.samples)
					break
		return

	def extract(self):
		"""@yield: chrom, pos, rs, ref, alt, qual, filt, info, fmt, samples"""
		with gzip.open(self.vcf, 'r') as infile:
			for line in infile:
				line = line.decode()
				if line.startswith("#"): continue  # release #CHROM
				yield(Vcfline(line))


class Tabix(object):
	"""read database info, return maybe multilines
	每行比subprocess节约0.05s"""
	def __init__(self, database, chrom, start, end, prefix=False):
		"""
		prefix:  chromosome prefix, True - 'chr1', False - '1'
		@param database:             combined database file
		@param chrom:                chromosome
		@param start:                start pos for tabix
		@param end:                  end pos for tabix
		@param prefix:               chromosome prefix.  'chrX': True; 'X': False
		"""
		self.database = database
		if prefix:
			self.chrom = chrom
		elif chrom == 'chrM':
			self.chrom = 'MT'
		else:
			self.chrom = chrom.replace('chr','')
		self.start = int(start)
		self.end = int(end)
		#print(self.chrom, self.start, self.end)
		self.extract = self.extract()

	def extract(self):
		"""feedbacks is str(maybe contain '\\n', multilines) or None
		@return feedbacks:     database result
		"""
		s = datetime.datetime.now()
		try:
			feedbacks = pysam.TabixFile(self.database).fetch(self.chrom, self.start-1, self.end)
		except:
			raise
		#print(feedbacks)  # '\n'
		#feedbacks = feedbacks.decode('UTF-8').strip() if feedbacks else ''
		e = datetime.datetime.now()
		#logging.error("tabix time: %s" % (e-s))
		return feedbacks

	def extract_discard(self):
		"""feedbacks is str(maybe contain '\\n', multilines) or None
		@return feedbacks:     database result
		"""
		s = datetime.datetime.now()
		feedbacks = subprocess.check_output(['tabix', self.database, '{}:{}-{}'.format(self.chrom, self.start, self.end)], shell=False)
		#print(feedbacks)  # '\n'
		feedbacks = feedbacks.decode('UTF-8').strip() if feedbacks else ''
		e = datetime.datetime.now()
		logging.error("tabix time: %s" % (e-s))
		return feedbacks
