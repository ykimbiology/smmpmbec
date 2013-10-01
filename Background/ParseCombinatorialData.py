import math

# The input text file JSCombinatorialCompilation.txt contains raw IC50 values
# from completed combinatorial libraries that 'worked', i.e. that showed at
# least some binding for something, were repeated etc. 

lines = open("RawData/JSCombinatorialCompilationIC50.txt","rt").readlines()
header = lines[0].strip().split("\t")
residues = "ACDEFGHIKLMNPQRSTVWY"

# This script simply reformats that file, removes the inequalities
# and converts to log10(ic50), to make it easier for Mathematica to handle

vals = {}
count =0
for allele in range(5, len(header)):
    print header[allele], 
    for pos in range(9):
        print pos+1, 
        log_ic_list = []
        for line in range(20):
            f = lines[1+20*pos + line].strip().split("\t")
            assert(int(f[4])==pos+1)
            assert(residues[len(log_ic_list)]==f[3])
            meas = f[allele]
            if meas[0]=="<":
                ineq = True
                ic50 = math.log10(float(meas[1:]))
            else:
                ineq = False
                ic50 = math.log10(float(meas))
            log_ic_list.append(ic50)
        vals[header[allele],pos+1] = log_ic_list
        count+=1
    print ""
print count   

out = open("FormattedCombinatorials_log_ic50.txt","wt")
out.write("Allele\tPosition")
for i in range(20):
    out.write("\t%s" % residues[i])
for (allele, pos) in sorted(vals.keys()):
    out.write("\n%s\t%s" % (allele, pos))
    for v in vals[allele,pos]:
        out.write("\t%s" % v)
out.close()
