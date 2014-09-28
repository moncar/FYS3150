from numpy import *
from matplotlib import *
from matplotlib import pyplot as plt
from os import system
infile = open("/home/filiphl/Desktop/v-er.txt", "r")
line = infile.readlines()
v = []
for word in line:	#Makes an array of the v-values.
	v.append(float(word))
system("rm /home/filiphl/Desktop/v-er.txt") #Deletes the .txt file.


N = int(v[-1])
v = v[:-1]

h = 1./(N+1)
u = []
x = []
for i in range(N):
	x.append(i*h)
	u.append(1-(1-exp(-10))*x[i]-exp(-10*x[i]))

plt.plot(x,u, "-r" )
plt.hold("on")
plt.plot(x,v)
plt.legend(["u(x)", "v"])
plt.title("N = %d" %N)
plt.xlabel("x")
plt.savefig("/home/filiphl/Desktop/fhl_p2_%d.jpeg" %N)	#Stores the figures.
plt.show()


