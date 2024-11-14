from sys import argv

input_file = open(argv[1])

for line in input_file:
    if "#" in line:
        if "CHROM" in line:
            print('\t'.join(line.strip().split()[:8]))
        else:
            print(line.strip())
    else:
        fields=line.strip().split()
        print(fields[0]+'\t'+fields[1]+'\t.\t'+fields[3]+'\t'+fields[4]+'\t.\t.\tAF=0.5')
