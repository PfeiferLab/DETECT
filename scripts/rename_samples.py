from sys import argv
import sys
import gzip
if argv[1].endswith('.gz'):
    input_file = gzip.open(argv[1],'rt')
else:
    input_file = open(argv[1],'r')
old_names = argv[2].split(",")

parent_1 = old_names[0]
parent_2 = old_names[1]
child = old_names[2]

for line in input_file:
        if "#CHROM" in line:
            names = line.strip().split()[-3:]
            for i in range(0,len(names)):
                if names[i] == parent_1:
                    names[i] = "parent_1"
                elif names[i] == parent_2:
                    names[i] = "parent_2"
                elif names[i] == child:
                    names[i] = "child"
                else:
                    print("Names do not match samples in the VCF!")
                    sys.exit()
            print("\t".join(line.strip().split()[:-3]) + '\t'+ '\t'.join(names))
        else:
            print(line.strip()) 
