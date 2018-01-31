#!/usr/bin/python
import sys
import glob
script, filename = sys.argv
script = sys.argv[0]
filename = sys.argv[1]


# Author: Jerelle A Joseph jaj52@cam.ac.uk
# If you want to colour the free energy (FE) disconnectivity graph based on some order parameter 
# you need to get the potential energy (PE)  minimum that is representative of each FE group. 
#
# ***USAGE***
# python get_FreeEnergy_GroupRep.py minima_groups.<TEMPERATURE>
#
# ***FILE NEEDED***
# minima_groups.<TEMPERATURE> - output file from PATHSAMPLE specifying the composition of each FE group 
#
# This script will create:
# 1. FEgroup_PErep.out - specifies the FE group index (in min.data.regrouped.<TEMEPERATURE>) 
#                        and the representative PE minimum index (in min.data)
# 2. PErep.out - list the PE minimum indices (one per line) for each FE group
#
# 3. The script can also dump the free energy groups into individual files and also produce a file that gives 
# a summary of number of minima in each FE group ("GROUP_COUNT"). Uncomment the lines at the end of the script if you
# want this. 
#
# ***BEYOND THE SCRIPT***
#
# You will then need to compute the order parameter ONLY for the minima specified in "PErep.out"
# Then, as for PE disconnectivity graphs, provide the values for your order parameter in a file listed the same order
# as the PE mininma appear in "PErep.out" and include the necessary key words in your dinfo file. 
# The values will then match the FE groups once you DO NOT change the order! 
#
# ***CAUTION***
#
# The order in which the groups are lised in minima_groups.<TEMEPERATURE> match the order in which they appear
# in min.data.regrouped.<TEMPERATURE> but NOT necessarily how they are listed in min.data.regrouped.resorted!!!
# Hence, using min.data.regrouped.<TEMPERATURE> for your analyses is highly recommended!!


def get_rep(fileofgroups):
    '''Get a PE minimum representative for each FE group'''
    output1 = open('FEgroup_PErep.out', 'w')
    output1.write('FE group      PE rep')
    output1.write('\n')

    with open(fileofgroups) as f:
        lines = f.readlines()

    this_group = []
    for index, line in enumerate(lines):
        if not 'group' in line: 
            if line.strip() != '':
                this_group.append(line)
        elif 'group' in line:
            header = line.split()
            groupnum = header[1]
            rep = this_group[0]   
            PE_rep.append(rep)     
            output1.write('%s      %s' %(groupnum,rep))
            output1.write('\n')
            this_group = []    
    return 


def dump_groups(minima_groups):
    '''Dumps each free energy group in a separate file. Useful for colouring the 
        PE disconnectivity graph based on the composition of specific FE groups'''
    with open(minima_groups) as f:
        lines = f.readlines()
	
    currset=[]
    for index, line in enumerate(lines):
        if not 'group' in line:
            currset.append(line)
        elif 'group' in line:
            header = line.split()
            groupnum = header[1]
            newgroup = open('group_%s' % groupnum, 'w')
            newgroup.write(line)
            
            for min_ in currset:
                newgroup.write(min_)
            newgroup.close()
            currset=[]
    return 



def group_info(singlegroup):
    '''Gets the number of minima present in each group and the group number'''
    with open(singlegroup) as f:
        first_line = f.readline()
    data = first_line.split()
    group_number = data [1]
    num_of_min = data[7]
    
    return [group_number, num_of_min]



def group_count():
    '''Produces a file with the FE group number and the number of minima in each'''
    output = open('GROUP_COUNT', 'w')

    for file_ in glob.glob('group_*'):
        group_data = group_info(file_)
        num_min = group_data[1]
        group_num = group_data[0]
        output.write(str(group_num))
        output.write(' ')
        output.write(str(num_min))
        output.write('\n')
    return 

# CODE TO RUN
PE_rep = []

get_rep(filename)

output2 = open('PErep.out', 'w')

for r in PE_rep:
    output2.write(r)

# Uncomment this line if you want to dump the free energy groups into separate files 
############## dump_groups(filename)

# Uncommment this line if you want to get a summary of ONLY the number of PE minima
# in each free group. NOTE you must uncomment the <dump_groups(fileofgroups)> line above
# if you want to use this
############# #group_count()
