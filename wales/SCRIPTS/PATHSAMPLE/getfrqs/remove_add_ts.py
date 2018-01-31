ts_f = "remove_add.ts"
ts_remove = []
with open(ts_f , "r") as _:
   for line in _:
      try:
         ts_remove.append(int(line.strip()))
      except ValueError:
         pass

input_f = "frqs.sorted.ts"
out_f = open("frqs.sorted.ts.new","w")


with open(input_f , "r") as f:
   for line in f:
      id = int(line.split()[0])
      if id in ts_remove:
         pass
      else:
         out_f.write(line)
out_f.close()
