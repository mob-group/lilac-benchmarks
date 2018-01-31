#!/usr/bin/env python
import networkx as nx
import numpy as np
import sys
import math

def duplicate_test(graph,min1,min2):
    return min1 in G.neighbors(min2)
if len(sys.argv) < 4:
   print 'dijkstra_test.py minA minBi npairs ([True/False] for debugging)'
   exit()
   
A = int(sys.argv[1])
B = int(sys.argv[2])
npairs = int(sys.argv[3])
if len(sys.argv) == 5: 
   if sys.argv[4] == 'True':
      debug = True
      print 'Debugging is on'
   else:
      debug = False
else:
   debug = False

#default for input files
min_f = "min.data"
ts_f = "ts.data"
paird_f = "pairdist"
pairl_f = "pairlist"

print
print '##############################################################################'
print '|Script to test whether A and B are connected by pairlist and ts.data        |'
print '|If there is no connection, PATHSAMPLE will not be able attempt connections  |'
print '|other than connecting A and B directly.                                     |'
print '|If there is a connection, the weight given is the sum of the distances on   |' 
print '|the shortest dijkstra path without a weighting function.                    |'
print '|Written by Konstantin Roeder (kr366)                                        |'
print '##############################################################################'
print

#read files into arrays
mindata = np.genfromtxt(min_f , dtype="float,float,int,float,float,float")
tsdata = np.genfromtxt(ts_f , dtype="float,float,int,int,int,float,float,float")
paird = np.genfromtxt(paird_f , dtype = float)
pairl = np.genfromtxt(pairl_f , dtype = int)

#reshape paird and pairl to 100 entries per minima as it ought to be (genfromtxt reads it as 10x10 for each minima)
paird.shape = paird.shape[0]/(npairs/10),npairs
pairl.shape = pairl.shape[0]/(npairs/10),npairs

if (paird.shape[0] != pairl.shape[0]) or (paird.shape[0] != mindata.shape[0]):
   raise Exception('Number of minima and number of pair entries are not equal.')
   

nmin = mindata.shape[0]
nts = tsdata.shape[0]

#print 'There are',nmin,'minima and',nts,'transition states.'

#create graph and add in all minima and transition states as edges with weight 0.0
G = nx.Graph()

for minima in xrange(nmin):
    G.add_node(minima + 1,emin=mindata[minima][0])
    if debug:
       print 'Read minimum' , minima + 1


for ts in xrange(nts):
    if (tsdata[ts][3] != tsdata[ts][4]):
       if not duplicate_test(G,tsdata[ts][3],tsdata[ts][4]):
          G.add_edge(tsdata[ts][3],tsdata[ts][4],weight=0.0,ets=tsdata[ts][0])
          if debug:
             print 'Read TS between minima',tsdata[ts][3],'and',tsdata[ts][4]
   # else:
      # print 'TS',ts+1,'connects minimum',tsdata[ts][3],'with itself'

#add egdes for every pair in pairl and use paird as weight

for minima in xrange(nmin):
    for entry in xrange(len(pairl[minima])):
        if not duplicate_test(G , minima + 1 , pairl[minima][entry]):
           G.add_edge(minima + 1 , pairl[minima][entry] , weight = paird[minima][entry])
          # G.add_edge(minima + 1 , pairl[minima][entry] , weight = math.exp(paird[minima][entry]))
           if debug:
              print 'Added edge between minima',minima + 1,'and',pairl[minima][entry]
try:
   wtot = nx.dijkstra_path_length(G,A,B)
   print 'The shortest dijkstra path has weight',wtot
   path = nx.dijkstra_path(G,A,B)
   sum_ = 0 
   for node in xrange(len(path)-1):
      sum_ += math.exp(G[path[node]][path[node+1]]['weight'])
      print '%6d %6d %12.4f  %16.8g' %(path[node], path[node+1], G[path[node]][path[node+1]]['weight'], sum_)
except nx.exception.NetworkXNoPath:
   print 'The A and B sets are not connected by ts.data and pairlist.'
