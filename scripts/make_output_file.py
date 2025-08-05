import argparse
import numpy as np
import statistics as st
import pandas as pd

parser = argparse.ArgumentParser(description="Creates a config file, which is then passed to the DETECT workflow")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-i","--input-file",dest="input_file",help="input mutation table.",required = True)
required_args.add_argument("-r","--reassembly-file",dest="reassembly_file",help="input mutation table.",required=True)
required_args.add_argument("-o","--output",dest="output_file",help="name of the output file, specified by the config file.",required = True)
args = parser.parse_args()

sd_list = np.arange(0.5,5.1,0.5)

output_file = open(args.output_file,'w')

input_file = pd.read_csv(args.input_file,sep='\t')
header_line = 'Filter'+'\t'+'min/max'+'\t'+'Average'+'\t'+'SD'+'\t'+'\t'.join([str(x) for x in sd_list])
output_file.write(header_line+'\n')
#QUAL
ss='QUAL'
input_file_stats = input_file[ss]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write('QUAL\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')

#GQ
for name in ['child','parent_1','parent_2']:

    ss=name+'.GQ'
    print(ss)
    input_file_stats = input_file[ss]
    mean = st.mean(input_file_stats)
    sd = st.stdev(input_file_stats)

    output_file.write(name+'.GQ\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')

#DP
for name in ['child','parent_1','parent_2']:
    ss=name+'.DP'
    print(ss)
    input_file_stats = input_file[ss]
    mean = st.mean(input_file_stats)
    sd = st.stdev(input_file_stats)
    output_file.write(name+'.DP\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')
    output_file.write(name+'.DP\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')

#AD
for name in ['parent_1','parent_2']:
    ss=name+'.AD'
    print(ss)
    input_file_stats = [int(x.split(',')[1]) for x in input_file[ss]]
    mean = st.mean(input_file_stats)
    sd = st.stdev(input_file_stats)
    output_file.write(name+'.AD\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')

#AB
name='child'
ss=name+'.AD'
print(ss)
input_file_stats = [float(int(x.split(',')[1])/(int(x.split(',')[0])+int(x.split(',')[1]))) for x in input_file[ss]]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write(name+'.AB\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')
output_file.write(name+'.AB\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')

#QD
ss='QD'
print(ss)
input_file_stats = input_file[ss]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write('QD\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')

#MQRankSum
ss='MQRankSum'
print(ss)
input_file_stats = input_file[ss]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write('MQRankSum\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')
output_file.write('MQRankSum\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')

#ReadPosRankSum
ss='ReadPosRankSum'
print(ss)
input_file_stats = input_file[ss]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write('ReadPosRankSum\tmin\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean-x) for x in sd*sd_list])+'\n')
output_file.write('ReadPosRankSum\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')


#FS
ss='FS'
print(ss)
input_file_stats = input_file[ss]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write('FS\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')

#SOR
ss='SOR'
input_file_stats = input_file[ss]
mean = st.mean(input_file_stats)
sd = st.stdev(input_file_stats)
output_file.write('SOR\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')

#Reassembly
reassembly_file = pd.read_csv(args.reassembly_file,sep='\t')
for name in ['parent_1','parent_2','child']:
    ss=name+'_reassembled'
    reassembly_file_stats = reassembly_file[ss]
    mean = st.mean(reassembly_file_stats)
    sd = st.stdev(reassembly_file_stats)
    output_file.write(name+'.reassembly\tmax\t'+str(mean)+'\t'+str(sd)+'\t'+'\t'.join([str(mean+x) for x in sd*sd_list])+'\n')
