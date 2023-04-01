import argparse
import numpy as np
import statistics

parser = argparse.ArgumentParser(description="Creates a config file, which is then passed to the DETECT workflow")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-n","--num-iterations",dest="num_iterations",help="Number of iterations specified by the user.",required = True)
required_args.add_argument("-o","--output",dest="output_file",help="name of the output file, specified by the config file.",required = True)
args = parser.parse_args()

iter_list = [str(x) for x in np.arange(1,int(args.num_iterations)+1)]

output_file = open(args.output_file,'w')


header_line = 'Filter'+'\t'+'Average'+'\t'+'Variance'+'\t'+'\t'.join(iter_list)
output_file.write(header_line+'\n')

iter_dict = {}

for item in iter_list:
    input_file = open('pipeline/outputs/output_'+item+'.txt')
    for line in input_file:
        if "Threshold" not in line and "original" not in line and "MVs" not in line and "Multihit" not in line:
            fields = line.strip().split()
            filter_name=fields[0]
            recommended_value = float(fields[1])
            if filter_name not in iter_dict.keys():
                iter_dict[filter_name] =[]
            iter_dict[filter_name].append(recommended_value)
        
        if "Multihit" in line:
            fields = line.strip().split(":")
            if fields[0] not in iter_dict.keys():
                iter_dict[fields[0]] =[]
            iter_dict[fields[0]].append(float(fields[1]))
for key in iter_dict.keys():
    print(iter_dict[key])
    mean=statistics.mean(iter_dict[key])
    if args.num_iterations != "1":
        variance=statistics.variance(iter_dict[key])
    else:
        variance=0
    output_file.write(key+'\t'+str(mean)+'\t'+str(variance)+'\t'+'\t'.join([str(x) for x in iter_dict[key]])+'\n')
