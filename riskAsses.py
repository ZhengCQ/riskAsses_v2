#!/usr/bin/env  python
# -*- coding:UTF-8 -*-

import argparse
import sys
import os
from threading import Thread
from bin.riskAssesing import *


ARGS = argparse.ArgumentParser(description="危害度评估")
ARGS.add_argument(
	'-v', '--vcf', dest='vcf', required=True, help='snpEff注释后的vcf')
ARGS.add_argument(
	'-w', '--work_dir', dest='work_dir', default='.', help='工作目录，默认：当前目录。')
ARGS.add_argument(
	'-A', dest='A_population', required=True, help='A population的名称')
ARGS.add_argument(
	'-A_samples',  dest='A_samples', required=True, 
	help='A population的具体样本名称，多个样本以逗号分隔，如sample1,sample2')
ARGS.add_argument(
	'-B', dest='B_population', required=True, help='B population的名称')
ARGS.add_argument(
	'-B_samples',  dest='B_samples', required=True, 
	help='B population的具体样本名称，多个样本以逗号分隔，如sample1,sample2')
ARGS.add_argument(
	'-C_samples',  dest='C_samples', required=True, 
	help='C population（外群）的具体样本名称，多个样本以逗号分隔，如sample1,sample2')

def get_all_sites(vcf, work_dir):
	if vcf.endswith('gz'):
		os.system("vcftools --gzvcf %s --out %s/vcfstat"%(vcf, work_dir))
	else:
		os.system("vcftools --vcf %s --out %s/vcfstat"%(vcf, work_dir))
	f = open ('%s/vcfstat.log'%(work_dir))
	all_sites = re.findall(r'possible (\d+) Sites', f.read())[0] #用findall,不用search
	#os.system('rm ./vcfstat.log')
	return int(all_sites)

def main():
	args = ARGS.parse_args()
	if not args.vcf:
		print('Use --help for command line help')
		return
	A_samples_lst = args.A_samples.split(',')
	B_samples_lst = args.B_samples.split(',')
	C_samples_lst = args.C_samples.split(',')

	#all_sites = int(os.popen("less -S %s |wc -l" % (args.vcf)).read())
	#all_sites = get_all_sites(args.vcf, args.work_dir)
	all_sites = int(168773)
	#min_sites = int(10000)
	fix_sites = all_sites/2
	print 'There are %s sites'%(all_sites)
	outf = open('%s/block_riks.info'%args.work_dir, 'w')
	work = """
	"""
	for i in range(0, 10):
		B_samples_lst = random.sample(B_samples_lst, 4)
		risk = FuncRisk(args.vcf, args.A_population, args.B_population,
				A_samples_lst, B_samples_lst, C_samples_lst, fix_sites, all_sites)
		outf.write("{}\t{}\t{}\t{}\t{}\n".format(i, risk.sites, risk.risk_missense, risk.risk_synonymous, risk.risk_lof))
	outf.close()

if __name__ == '__main__':
	main()
