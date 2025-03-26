import os
import sys

rtf_file = sys.argv[1] # rtf file generated with cgenff

with open(rtf_file, 'r', errors='ignore') as f:
    lines = f.readlines()
f.close()

string_list = []
with open('tmpAtoms.prm', 'w+') as g:
    for line in lines:
        line_list = line.split()
        if len(line_list) >= 1 and line_list[0] == 'ATOM':
            str_write = line_list[2] + " " + line_list[3] + '\n'
            if str_write not in string_list:
                string_list.append(str_write)
                g.write(str_write)

g.close()

string_list2 = []
with open('lgAtoms.prm', 'w+') as h:
    for line in lines:
        line_list = line.split()
        if len(line_list) >= 1 and line_list[0] == 'ATOM':
            str_write2 = line_list[2] + '\n'
            if str_write2 not in string_list2:
                string_list2.append(str_write2)
                h.write(str_write2)

h.close()  

with open('atoms.prm', 'w+') as k:
    for line in lines:
        line_list = line.split()
        if len(line_list) >= 1 and line_list[0] == 'ATOM':
            str_write3 = line_list[1] + "  " + line_list[2] + "  " + line_list[3] + "\n"
            k.write(str_write3)