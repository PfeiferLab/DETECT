import subprocess
import argparse
import numpy as np
import statistics as st
import pandas as pd

parser = argparse.ArgumentParser(description="Creates a config file, which is then passed to the DETECT workflow")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-mv","--mv-vcf",dest="mv_vcf",help="MV VCF output from pipeline.",required=True)
required_args.add_argument("-v","--mutation-vcf",dest="mutation_vcf",help="original input mutation VCF.",required=True)
required_args.add_argument("-m","--mutation-table",dest="mutation_table",help="mutation table output by pipeline.",required = True)
required_args.add_argument("-p","--polymorphism-table",dest="poly_table",help="polymorphism table outout by pipeline.",required=True)
required_args.add_argument("-e","--error-table",dest="error_table",help="error table output by pipeline.",required=True)
required_args.add_argument("-rm","--reassembly-parent1-table",dest="reassembly_muts_table",help="mutation reassembly table.",required=True)
required_args.add_argument("-rp","--reassembly-parent2-table",dest="reassembly_poly_table",help="polymorphism reassembly table.",required=True)
required_args.add_argument("-re","--reassembly-child-table",dest="reassembly_error_table",help="error reassembly table.",required=True)
required_args.add_argument("-o","--output",dest="output_file",help="name of the output file, specified by the config file.",required = True)
args = parser.parse_args()

mv_vcf = args.mv_vcf
mutation_vcf = args.mutation_vcf
mutation_table = pd.read_csv(args.mutation_table,sep="\t")
poly_table = pd.read_csv(args.poly_table,sep='\t')
error_table = pd.read_csv(args.error_table,sep='\t')
reassembly_muts_table = pd.read_csv(args.reassembly_muts_table,sep=' ')
reassembly_poly_table = pd.read_csv(args.reassembly_poly_table,sep=' ')
reassembly_error_table = pd.read_csv(args.reassembly_error_table,sep=' ')
output_file = open(args.output_file,'w')

output_file.write('\t'.join(['Filter','min/max','Average','original_mutations','total_sites','total_mutations','total_mutation_mut_recall','total_mutation_precision','total_polymorphisms','total_other_sites','percentile_rec_5_95','filter_mutations','filter_mutation_recall','filter_mutation_precision','filter_polymorphisms','filter_other_sites'])+'\n')

def calc_precision(num_tp,num_fp):
    return num_tp/(num_tp+num_fp)

def calc_recall(num_tp,total_muts):
    return num_tp/total_muts

def filter_count(query_series,mut_rec,min_max):
    if min_max == 'min':
        filtered = query_series[query_series >= mut_rec]
    elif min_max == 'max':
        filtered = query_series[query_series <= mut_rec]
    else:
        input('not min or max')
    
    count = len(filtered)
    return count

original_mut_count = subprocess.run('grep -v "#" '+mutation_vcf+' | wc -l',shell=True,capture_output=True,text=True,check=True).stdout.strip()
original_mut_count = int(original_mut_count)

mv_total_sites = subprocess.run('grep -v "#" '+mv_vcf+' | wc -l',shell=True,capture_output=True,text=True,check=True).stdout.strip()
mv_total_sites = int(mv_total_sites)
mv_mut_count = len(mutation_table)
mv_poly_count = len(poly_table)
mv_error_count = len(error_table)
mv_mut_recall = calc_recall(mv_mut_count,original_mut_count)
mv_mut_precision = calc_precision(mv_mut_count,mv_total_sites)


#QUAL
ss='QUAL'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.05)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'min') 

filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_precision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#GQ
for name in ['child','parent_1','parent_2']:

    ss=name+'.GQ'
    print(ss)
    mut_stats = mutation_table[ss]
    mut_avg = mut_stats.mean()
    mut_rec = mut_stats.quantile(0.05)

    poly_stats = poly_table[ss]
    filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
    error_stats = error_table[ss]
    filtered_error_count = filter_count(error_stats,mut_rec,'min')

    filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
    total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

    filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
    filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

    output_file.write('\t'.join([ss,'min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#DPmin
for name in ['child','parent_1','parent_2']:
    ss=name+'.DP'
    print(ss)
    mut_stats = mutation_table[ss]
    mut_avg = mut_stats.mean()
    mut_rec = mut_stats.quantile(0.05)

    poly_stats = poly_table[ss]
    filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
    error_stats = error_table[ss]
    filtered_error_count = filter_count(error_stats,mut_rec,'min')

    filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
    total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

    filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
    filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

    output_file.write('\t'.join([ss,'min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#DPmax
for name in ['child','parent_1','parent_2']:
    ss=name+'.DP'
    print(ss)
    mut_stats = mutation_table[ss]
    mut_avg = mut_stats.mean()
    mut_rec = mut_stats.quantile(0.95)

    poly_stats = poly_table[ss]
    filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
    error_stats = error_table[ss]
    filtered_error_count = filter_count(error_stats,mut_rec,'max')

    filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
    total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

    filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
    filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

    output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')
 
#AD
for name in ['parent_2']: #['parent_1','parent_2']:
    ss=name+'.AD'
    print(ss)
    mut_stats = pd.Series([int(x.split(',')[1]) for x in mutation_table[ss]])
    mut_avg = mut_stats.mean()
    mut_rec = mut_stats.quantile(0.95)
    #print(type(poly_table[ss][0].split(',')[1]))
    bad_rows = poly_table[~poly_table[ss].astype(str).str.contains(",", na=False)]
    print(bad_rows)
    #stop

    poly_stats = pd.Series([int(x.split(',')[1]) for x in poly_table[ss]])
    filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
    error_stats = pd.Series([int(x.split(',')[1]) for x in error_table[ss]])
    filtered_error_count = filter_count(error_stats,mut_rec,'max')

    filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
    total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

    filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
    filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

    output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#ABmin
name='child'
ss=name+'.AD'
print(ss)
mut_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in mutation_table[ss]])
print(mut_stats)
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.95)

#ab_list=[]
#for x in poly_table[ss]:
#    ref,alt=x.split(',')
#    ref=float(ref)
#    alt=float(alt)
    
#    print(float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))))


#print([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in poly_table[ss]])
poly_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in poly_table[ss]])
#poly_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in poly_table[ss]])
filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
error_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in error_table[ss]])
filtered_error_count = filter_count(error_stats,mut_rec,'max')

filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([name+'.AB','max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#ABmax
name='child'
ss=name+'.AD'
print(ss)
mut_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in mutation_table[ss]])
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.05)

poly_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in poly_table[ss]])
filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
error_stats = pd.Series([float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in error_table[ss]])
filtered_error_count = filter_count(error_stats,mut_rec,'min')

filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([name+'.AB','min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#QD
ss='QD'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.05)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'min')

filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#MQRankSum min
ss='MQRankSum'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.05)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'min')

filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#MQRankSum max
ss='MQRankSum'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.95)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'max')

filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#ReadPosRankSum min
ss='ReadPosRankSum'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.05)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'min')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'min')

filtered_mut_count = filter_count(mut_stats,mut_rec,'min')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'min',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#ReadPosRankSum min
ss='ReadPosRankSum'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.95)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'max')

filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#FS
ss='FS'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.95)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'max')

filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#SOR
ss='SOR'
print(ss)
mut_stats = mutation_table[ss]
mut_avg = mut_stats.mean()
mut_rec = mut_stats.quantile(0.95)

poly_stats = poly_table[ss]
filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
error_stats = error_table[ss]
filtered_error_count = filter_count(error_stats,mut_rec,'max')

filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')

#Reassembly
for name in ['p1','p2','child']:
    ss=name+'_reassembled'
    print(ss)
    mut_stats = reassembly_muts_table[ss]
    mut_avg = mut_stats.mean()
    mut_rec = mut_stats.quantile(0.95)

    poly_stats = reassembly_poly_table[ss]
    filtered_poly_count = filter_count(poly_stats,mut_rec,'max')
    error_stats = reassembly_error_table[ss]
    filtered_error_count = filter_count(error_stats,mut_rec,'max')

    filtered_mut_count = filter_count(mut_stats,mut_rec,'max')
    total_filtered_count = filtered_mut_count + filtered_poly_count + filtered_error_count

    filtered_mut_recall = calc_recall(filtered_mut_count,original_mut_count)
    filtered_mut_pmut_recision = calc_precision(filtered_mut_count,total_filtered_count)

    output_file.write('\t'.join([ss,'max',str(mut_avg),str(original_mut_count),str(mv_total_sites),str(mv_mut_count),str(mv_mut_recall),str(mv_mut_precision),str(mv_poly_count),str(mv_error_count),str(mut_rec),str(filtered_mut_count),str(filtered_mut_recall),str(filtered_mut_precision),str(filtered_poly_count),str(filtered_error_count)])+'\n')
