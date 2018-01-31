import sys

f_in = sys.argv[1]
n_sp = int(sys.argv[2])
dic_frqs = {}

for i in xrange(n_sp):
    dic_frqs[i+1] = []


with open(f_in , "r") as f:
    for line in f:
        sp_id = int(line.split()[0])
        dic_frqs[sp_id].append(float(line.split()[1]))

out_f = open('frqs.sorted' , "w")
removed = open('sp.remove' , "w")
remove = 0
string = ''
for entry in dic_frqs.keys():
    if len(dic_frqs[entry]) >= 1:
        frq = str("%20.10f") % dic_frqs[entry][0]
        id_sp = str("%8i") % entry
        out_f.write(id_sp + '      ' + frq + "\n")
    else:
        remove += 1 
        string += str("%8i") % entry
        string += "\n"
        print 'No frequency found for stationary point: ',entry

out_f.close()

if remove > 0:
    removed = open('sp.remove' , "w")
    removed.write(str(remove) + "\n" + string)
    removed.close()


