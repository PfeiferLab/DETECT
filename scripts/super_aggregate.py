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


header_line = 'Filter'+'\t'+'min/max'+'\t'+'Average'+'\t'+'Variance'+'\t'+'\t'.join(iter_list)
output_file.write(header_line+'\n')

iter_dict = {}

for run in iter_list:
    input_file = open('pipeline/run_outputs/run'+run+'.statistics.txt')
    for line in input_file:
        if "FilterName" not in line:
            fields = line.strip().split()
            filter_name=fields[1]
            bound=fields[2]
            print(iter_dict)
            recommended_value = float(fields[3])
            if filter_name not in iter_dict.keys():
                iter_dict[filter_name] = {}
            if bound not in iter_dict[filter_name].keys():
                iter_dict[filter_name][bound] = []

            iter_dict[filter_name][bound].append(recommended_value)
for filtername in iter_dict.keys():
    for bound in iter_dict[filtername].keys():
        print(iter_dict[filtername][bound])
        mean=statistics.mean(iter_dict[filtername][bound])
        if args.num_iterations != "1":
            variance=statistics.variance(iter_dict[filtername][bound])
        else:
            variance=0
        output_file.write(filtername+'\t'+bound+'\t'+str(round(mean,2))+'\t'+str(round(variance,2))+'\t'+'\t'.join([str(x) for x in iter_dict[filtername][bound]])+'\n')
