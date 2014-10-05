from numpy import *
from matplotlib import *
from matplotlib import pyplot as plt

#Reading in the data.
infile = open("/home/filiphl/Desktop/project3/PlanetData.txt", "r")
n = int(infile.readline().split()[1])	#Number of objects.
names = infile.readline().split()	#Stores their names, in order to show in the plot.
for i in range(len(names)):
	names[i] = names[i][:-2]  	#Removes the "_x" and "_y" from the names, for the plot.
planets = []	
for i in range(2*n):
	planets.append([])		#Setting up arrays for the data

time = []
for line in infile:			#Stores the data in arrays.
	words = line.split()
	time.append(float(words[0]))	#Time. Not that I ever use it...LOL.
	for i in range(2*n):
		planets[i].append(float(words[i+1]))

#Plot of every planet/star/whatever 
for i in range(n):
	if i == 0:
		plt.plot(planets[i*2], planets[i*2 +1], '*')
	else:
		plt.plot(planets[i*2], planets[i*2 +1], '-')
	plt.hold('on')
plt.xlabel('x / AU')
plt.ylabel('y / AU')
plt.legend(names[1::2],numpoints=1)	#Skips the first element, and every other element.
plt.show()
	
