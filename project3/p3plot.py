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
		plt.plot(planets[i*2], planets[i*2 +1], '*', linewidth=2.0)
	else:
<<<<<<< HEAD
                plt.plot(planets[i*2], planets[i*2 +1], '-', linewidth=2.0)
=======
                plt.plot(planets[i*2], planets[i*2 +1], '-', linewidth=3.0)
>>>>>>> FETCH_HEAD
	plt.hold('on')
plt.xlabel('x / AU')
plt.ylabel('y / AU')
plt.legend(names[1::2],numpoints=1)	#Skips the first element, and every other element.
plt.axis('equal')
yrs = time[-1]/(24.0*365.0)
plt.title('Time period: %.1f years.' %yrs)
plt.grid()

#Plotting the energies...
#Reading in the data.
infile = open("/home/filiphl/Desktop/project3/PlanetEnergy.txt", "r")

time = []
Ek = []
Ep = []
E = []
infile.readline()
for line in infile:			#Stores the data in arrays.
	words = line.split()
<<<<<<< HEAD
	time.append(float(words[0])/(365.*24.))	#Time. Not that I ever use it...LOL.
=======
	time.append(float(words[0]))	#Time. Not that I ever use it...LOL.
>>>>>>> FETCH_HEAD
	Ek.append(float(words[1]))
	Ep.append(float(words[2]))
	E.append(float(words[3]))

plt.figure(2)
<<<<<<< HEAD
plt.plot(time,Ek, '-', linewidth=2.0)
plt.hold("on")
plt.plot(time,Ep, '-', linewidth=2.0)
plt.plot(time,E, '-', linewidth=2.0)	
=======
plt.plot(time,Ek)
plt.hold("on")
plt.plot(time,Ep)
plt.plot(time,E)	
>>>>>>> FETCH_HEAD
plt.xlabel('Years')
plt.ylabel('Energy')
plt.title('Energy of the system')
plt.legend(["Kinetic","Potential", "Total"])
<<<<<<< HEAD
plt.axis([min(time), max(time), min(Ep)*1.1, max(Ek)*1.1])
=======
>>>>>>> FETCH_HEAD

plt.figure(3)
plt.subplot(3,1,1)
plt.plot(time,Ek)
<<<<<<< HEAD
plt.legend(["Kinetic"])
plt.subplot(3,1,2)
plt.plot(time,Ep)
plt.ylabel('Energy')
plt.legend(["Potential"])
plt.subplot(3,1,3)
plt.plot(time,E)	
plt.xlabel('Years')
=======
plt.xlabel('Years')
plt.ylabel('Energy')
plt.legend(["Kinetic"])
plt.subplot(3,1,2)
plt.plot(time,Ep)
plt.legend(["Potential"])
plt.subplot(3,1,3)
plt.plot(time,E)	
>>>>>>> FETCH_HEAD
plt.legend(["Total"])

plt.show()
	
