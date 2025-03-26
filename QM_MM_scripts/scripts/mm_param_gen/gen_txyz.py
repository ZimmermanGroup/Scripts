import os
import sys

old_txyz = sys.argv[1] #old_txyz file with MM2 types
matched_prm = sys.argv[2] #matched_atoms.prm
new_txyz = sys.argv[3]  #new txyz file name

with open(old_txyz, 'r+') as f:
    lines_txyz = f.readlines()
f.close()

with open(matched_prm, 'r+') as g:
    lines_prm = g.readlines()
g.close()

with open(new_txyz, 'w+') as h:
    h.write(lines_txyz[-1].split()[0])
    h.write("\n\n")
    for idx, line in enumerate(lines_txyz[1:]):
        line_list = line.split()
        prm_line_list = lines_prm[idx].split()
        line_list[5] = prm_line_list[2]
        h.write('  '.join(line_list))
        h.write("\n")

h.close() 