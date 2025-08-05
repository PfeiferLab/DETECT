import pandas as pd
import argparse
import numpy as np
import sys
import os

#-m -v -e

parser = argparse.ArgumentParser(description="Creates a config file, which is then passed to the DETECT workflow")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-n","--num-iter",dest="num_iter",help="Iteration Number",required = True)
#required_args.add_argument("-m","--mut-file",dest="mut_file",help="Mutation table containing mutations",required = True)
#required_args.add_argument("-v","--variant-file",dest="variant_file",help="Mutation table containing polymorphisms", required = True)
#required_args.add_argument("-e","--error-file",dest="error_file",help="Mutation table containing errors (not mutations or polymorphisms)", required = True)
required_args.add_argument("-o","--out-file",dest="out_file",help="Output file", required = True)
args = parser.parse_args()

#Set argparse Variables
#mut_table = pd.read_table(args.mut_file)
#var_table = pd.read_table(args.variant_file)
#err_table = pd.read_table(args.error_file)
output_file = open(args.out_file,'w')

#>= is min, <= is max

if not os.path.exists('pipeline/run_outputs'):
    os.makedirs('pipeline/run_outputs')

run = int(args.num_iter)
rec_filter_output_dict = {}

def return_counts(run,filtername,bound):
    mut_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.mutations.table')
    var_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.polymorphisms.table')
    err_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.errors.table')
    if bound == 'max':
        val = max(mut_table[filtername]) 
        mut_counts = len(mut_table[mut_table[filtername] <= val])
        var_counts = len(var_table[var_table[filtername] <= val])
        err_counts = len(err_table[err_table[filtername] <= val])
    elif bound == 'min':
        val = min(mut_table[filtername])
        mut_counts = len(mut_table[mut_table[filtername] >= val])
        var_counts = len(var_table[var_table[filtername] >= val])
        err_counts = len(err_table[err_table[filtername] >= val])
    else:
        print("ERROR")
        sys.exit()
    total_counts = mut_counts + var_counts + err_counts
    return [run,filtername,bound,val,total_counts,mut_counts,var_counts,err_counts]

def return_AD(run,filtername,bound):
    mut_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.mutations.table')
    var_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.polymorphisms.table')
    err_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.errors.table')

    mut_list = []
    var_list = []
    err_list = []

    for site in mut_table[filtername]:
        ref,alt = site.split(',')
        alt = int(alt)
        mut_list.append(alt)

    for site in var_table[filtername]:
        ref,alt = site.split(',')
        alt = int(alt)
        var_list.append(alt)

    for site in err_table[filtername]:
        ref,alt = site.split(',')
        alt = int(alt)
        err_list.append(alt)

    best_ad = max(mut_list)
    mut_counts = len([x for x in mut_list if x <= best_ad])
    var_counts = len([x for x in var_list if x <= best_ad])
    err_counts = len([x for x in err_list if x <= best_ad])
    total_counts = mut_counts + var_counts + err_counts
    return [run,filtername,bound,best_ad,total_counts,mut_counts,var_counts,err_counts]

def return_AB(run,bound):
    mut_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.mutations.table')
    var_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.polymorphisms.table')
    err_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.errors.table')
    
    mut_list = []
    var_list = []
    err_list = []
    
    for site in mut_table['child.AD']:
        ref,alt = site.split(',')
        ref = int(ref)
        alt = int(alt)
        ab = round(alt/(ref+alt),2)
        mut_list.append(ab)
    
    for site in var_table['child.AD']:
        ref,alt = site.split(',')
        ref = int(ref)
        alt = int(alt)
        ab = alt/(ref+alt)
        var_list.append(alt/(ref+alt))
    
    for site in err_table['child.AD']:
        ref,alt = site.split(',')
        ref = int(ref)
        alt = int(alt)
        ab = alt/(ref+alt)
        err_list.append(alt/(ref+alt))

    if bound == 'max':
        best_ab = max(mut_list)
        mut_counts = len([x for x in mut_list if x <= best_ab])
        var_counts = len([x for x in var_list if x <= best_ab])
        err_counts = len([x for x in err_list if x <= best_ab])
    
    elif bound == 'min':
        best_ab = min(mut_list)
        mut_counts = len([x for x in mut_list if x >= best_ab])
        var_counts = len([x for x in var_list if x >= best_ab])
        err_counts = len([x for x in err_list if x >= best_ab])

    else:
        print("ERROR")
        sys.exit()
    total_counts = mut_counts + var_counts + err_counts
    return [run,'child.AB',bound,best_ab,total_counts,mut_counts,var_counts,err_counts]

total_output = []

run_output = [] 
num_muts = os.popen('grep -v "#" pipeline/mutations/input_mutations.'+str(run)+'.phased.vcf | wc -l').read()
mut_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.mutations.table')
#var_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.polymorphisms.table')
#err_table = pd.read_table('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.errors.table')
print('pipeline/MV/all_chr_trio.sorted.mark_dups.MV.'+str(run)+'.mutations.table')
print(mut_table.columns)
#DPMAX
max_dp_p1 = return_counts(run,'parent_1.DP','max') #max(mut_table['parent_1.DP'])
#print(mut_table['parent_1.DP'])
run_output.append(max_dp_p1)

max_dp_p2 = return_counts(run,'parent_2.DP','max') #max(mut_table['parent_2.DP'])
run_output.append(max_dp_p2)

max_dp_ch = return_counts(run,'child.DP','max') #max(mut_table['child.DP'])
run_output.append(max_dp_ch)

#DPMIN
min_dp_p1 = return_counts(run,'parent_1.DP','min')
run_output.append(min_dp_p1)

min_dp_p2 = return_counts(run,'parent_2.DP','min')
run_output.append(min_dp_p2)

min_dp_ch = return_counts(run,'child.DP','min')
run_output.append(min_dp_ch)

#GQ
gq_p1 = return_counts(run,'parent_1.GQ','min') #min(mut_table['parent_1.GQ'])
run_output.append(gq_p1)

gq_p2 = return_counts(run,'parent_2.GQ','min') #min(mut_table['parent_2.GQ'])
run_output.append(gq_p2)

gq_ch = return_counts(run,'child.GQ','min') #min(mut_table['child.GQ'])
run_output.append(gq_ch)

#QUAL
qual = return_counts(run,'QUAL','min') #min(mut_table['QUAL'])
run_output.append(qual)

#AD
ad_p1 = return_AD(run,'parent_1.AD','max') #mut_table['parent_1.AD']
run_output.append(ad_p1)

ad_p2 = return_AD(run,'parent_2.AD','max')
run_output.append(ad_p2)

ad_ch = list(mut_table['child.AD'])

#p1_ad = max([int(x.split(',')[1]) for x in ad_p1])
#p2_ad = max([int(x.split(',')[1]) for x in ad_p2])

#AB
ab_list = []
for item in ad_ch:
    ref,alt = item.strip().split(',')
    ref = float(ref)
    alt = float(alt)
    ab = alt/(ref+alt)
    ab_list.append(ab)

#ABMAX
best_abmax = return_AB(run,'max')
run_output.append(best_abmax)

#ABMIN
best_abmin = return_AB(run,'min')
run_output.append(best_abmin)

#QD
best_qd=return_counts(run,'QD','min')
run_output.append(best_qd)

#MQRSMAX
best_mqrsmax=return_counts(run,'MQRankSum','max')
run_output.append(best_mqrsmax)

#MQRSMIN
best_mqrsmin=return_counts(run,'MQRankSum','min')
run_output.append(best_mqrsmin)

#RPRSMAX
best_rprsmax=return_counts(run,'ReadPosRankSum','max')
run_output.append(best_rprsmax)

#RPRSMIN
best_rprsmin=return_counts(run,'ReadPosRankSum','min')
run_output.append(best_rprsmin)

#FS
best_fs = return_counts(run,'FS','max')
run_output.append(best_fs)

#SOR
best_sor = return_counts(run,'SOR','max')
run_output.append(best_sor)
#[run,filtername,bound,total_counts,mut_counts,var_counts,err_counts]
df = pd.DataFrame(run_output,columns = ['Run','FilterName','min/max','FilterValue','TotalSites','MutCounts','VarCounts','ErrCounts'])
print(df)
outfile = open('pipeline/run_outputs/run'+str(run)+'.statistics.txt','w')
df.to_csv(outfile,sep='\t',index=False)
#Conglomerate Runs
outfile.close()
