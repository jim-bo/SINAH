'''
translates clusters from opera to bundles
'''

# imports
import os
import sys

# parameters.
cluster_file = sys.argv[1]
bundle_file = sys.argv[2]
out_file = sys.argv[3]

# load valid bundles.
valid = set()
fin = open(cluster_file, "rb")
for line in fin:
    
    # tokenize.
    line = line.strip().split()
    ctga = line[0]
    ctgb = line[2]
    
    # add to set.
    valid.add( (ctga, ctgb) )
    valid.add( (ctgb, ctga) )
fin.close()
'''
First Contig    Orientation     Second Contig   Orientation     Distance        Standard Deviation      Size
2885    -       3035    +       -4554   250     1
2884    -       2986    -       -16760  250     1
2873    -       2888    -       -3398   250     1
2856    -       2924    -       -10843  250     1

212     2718    11.357419       0.000000        0.000000        0.000000
1695    2076    32.997731       0.000000        0.000000        0.000000
298     1777    0.000000        0.000000        81.000000       0.000000
1022    1713    0.000000        0.000000        20.966074       0.000000

'''

# filter the bundles.
fout = open(out_file, "wb")
fin = open(bundle_file, "rb")
for line in fin:
    
    # tokenize.
    tmp = line.strip().split()
    ctga = tmp[0]
    ctgb = tmp[1]
    
    # check.
    if (ctga, ctgb) in valid or (ctgb, ctga) in valid:
        fout.write(line)
fin.close()
fout.close()


