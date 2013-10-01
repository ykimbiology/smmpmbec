f = open("InverseCovariance.txt");
out = open("Code for InverseCovariance.txt", "wt");
row = 0
col = 0
for row in range(20):
    vals = f.readline().split("\t");
    out.write("\n\t");
    for col in range(20):
        out.write("m(%d,%d) = %s; " % (row, col, vals[col].strip()))
f.close()
f = open("InverseCovariance3.txt");

out.write("\n")

row = 0
col = 0
for row in range(20):
    vals = f.readline().split("\t");
    out.write("\n\t");
    for col in range(20):
        out.write("m(%d,%d) = %s; " % (row, col, vals[col].strip()))
out.close()
    
