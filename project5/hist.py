import numpy as np
from matplotlib import pyplot as plt
from sys import argv

filename = '/home/filiphl/Desktop/figs/hist1.txt'

infile = open(filename, "r")
firstline = infile.readline().split()
nrb = len(firstline)
nrp = int(firstline[0])
infile.close()

infile = open(filename, "r")
particles = []
i=0
for line in infile:
	w = line.split()
	if w[0] >0:
		particles.append([])
		for j in range(nrb):
			particles[i].append(float(w[j]))
	i+=1

Bar = plt.bar(range(nrb), particles[int(argv[1])])
plt.axis([0,nrb,0,particles[0][0]])
plt.show()



infile.close()
