#!/usr/bin/env  python
# -*- coding:UTF-8 -*-
# @Author: Zheng ChenQing
# @Date: 2019.04.30
# @E-mail: zhengchenqing@qq.com



import argparse
import sys
import os
from bin.riskAssesing import *
import time
import multiprocessing as mp

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
ARGS.add_argument(
    '--fixlength', dest='fixlength', type=int, 
    help='随机刀切时，刀切的固定长度, 数值大于总位点数时候，用总位点数')
ARGS.add_argument(
    '--A_samples_num', dest='A_samples_num', type=int, 
    help='A组随机样本数, 默认不随机，num大于或等于实际样本数时，按照A组实际样本输入')
ARGS.add_argument(
    '--B_samples_num', dest='B_samples_num', type=int,
    help='B组随机样本数，默认不随机，num大于或等于实际样本数时，按照B组实际样本输入')
ARGS.add_argument(
    '--n_core', dest='n_core', type=int,default=1,
    help='多进程数目, 默认为1')
ARGS.add_argument(
    '--n_permutation', dest='n_permutation', type=int, default=10,
    help='排列组合次数, 默认为10')

def get_all_sites(vcf, work_dir):
	if vcf.endswith('gz'):
		os.system("vcftools --gzvcf %s --out %s/vcfstat"%(vcf, work_dir))
	else:
		os.system("vcftools --vcf %s --out %s/vcfstat"%(vcf, work_dir))
	f = open ('%s/vcfstat.log'%(work_dir),'r')
	all_sites = re.findall(r'possible (\d+) Sites', f.read())[0] #用findall,不用search
	#os.system('rm ./vcfstat.log')
	return int(all_sites)

def write_table(func_dict, outfile):
	out_order = ['functional', 'splice_donor_variant', 'stop_lost', 'stop_gained', 'start_lost','splice_acceptor_variant',\
				 'lof','missense_variant','missense_deleterious','missense_benign','intron_variant','synonymous_variant', 'stop_retained_variant', \
	 			'intergenic_region', 'splice_region_variant', 'upstream_gene_variant','downstream_gene_variant',\
	  			'dn_ds', 'total']
	  			#'lof_intergenic', 'dn_intergenic', 'lof_ds',
	outf = open(outfile , 'w')
	for k in out_order:
		try:
			outf.write(k + '\t' +'\t'.join([ str(i) for i in func_dict[k]]) + '\n')
		except:
			print ('%s not exists in func_dict for %s'%(k, outfile))

def run_funcrisk(i, vcf, A_population, B_population,
				A_samples_lst, B_samples_lst, C_samples_lst, fix_sites, all_sites, work_dir):
	start = time.time()
	print('Run funcRisk task %s (%s)...' % (i, os.getpid()))
	funcstat = FuncRisk(i, vcf, A_population, B_population, A_samples_lst, B_samples_lst, 
				C_samples_lst, fix_sites, all_sites, work_dir)
	
	write_json('%s/Func_HomHet_Count_%s.json'%(work_dir,i),funcstat.func_homhet_count_dict)
	write_json('%s/Func_di_HomHet_Count_%s.json'%(work_dir,i),funcstat.func_di_homhet_count_dict)
	write_json('%s/Func_Count_%s.json'%(work_dir,i),funcstat.func_count_dict)
	write_json('%s/Func_di_Count_%s.json'%(work_dir,i),funcstat.func_di_count_dict)
	
	#write_table(funcstat.func_count_dict, work_dir + "/func_stat_info_%s.xls"%(i))
	#write_table(funcstat.func_di_count_dict, work_dir + "/func_derived_info_%s.xls"%(i))
	end = time.time()

	print('Task %s runs %0.2f seconds.' % (i, (end - start)))

def write_json(outf,dict_info):
	with open(outf, 'w') as fi:
		fi.write(json.dumps(dict_info))

def main():
	args = ARGS.parse_args()
	if not args.vcf:
		print('Use --help for command line help')
		return
	try:
		os.makedirs(args.work_dir)
	except:
		print ('%s exists' %(args.work_dir))
	A_samples_lst = args.A_samples.split(',')
	B_samples_lst = args.B_samples.split(',')
	C_samples_lst = args.C_samples.split(',')
	#all_sites = int(os.popen("less -S %s |wc -l" % (args.vcf)).read())
	#all_sites = get_all_sites(args.vcf, args.work_dir)
	all_sites = int(68775)
	fix_sites = int(args.fixlength) if args.fixlength else all_sites/2
	if all_sites -1000 > fix_sites:
		print ('There are %s sites, cut %s sites to analysis'%(all_sites, fix_sites))
	else:
		print ('There are %s sites, cut %s sites great than all sites, and turn off cut'%(all_sites, fix_sites))

	pool = mp.Pool(int(args.n_core)) #启动多线程池
	for i in range(0, args.n_permutation):
		i = i + 1
		if args.A_samples_num and args.A_samples_num < len(A_samples_lst):
			A_samples_lst = random.sample(A_samples_lst, args.A_samples_num) #随机A样本
			print ("启动A组样本随机，随机样本为%s"%(','.join(A_samples_lst)))
		if args.B_samples_num and args.B_samples_num < len(B_samples_lst):
			B_samples_lst = random.sample(B_samples_lst, args.B_samples_num) #随机B样本
			print ("启动B组样本随机，随机样本为%s"%(','.join(B_samples_lst)))
		run_funcrisk(i, args.vcf, args.A_population, args.B_population,
				A_samples_lst, B_samples_lst, C_samples_lst, fix_sites, all_sites, args.work_dir)	
		#pool.apply_async(run_funcrisk, args=(i, args.vcf, args.A_population, args.B_population,
		#		A_samples_lst, B_samples_lst, C_samples_lst, fix_sites, all_sites, args.work_dir)) #函数写入到多线程池
	print('Waiting for all subprocesses done...')
	pool.close()
	pool.join()
	print('All subprocesses done.')
	pool.terminate()

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print ("时间总计%s"%(end - start))
