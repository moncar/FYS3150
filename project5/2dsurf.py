from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sys import argv

num = argv[1] #int(raw_input('which time step? : '))
filename = '/home/filiphl/Desktop/figs/mat'+str(num)+'.txt'
#Reading in the data.
infile = open(filename, "r")
n = len(infile.readline().split())
infile.close()
infile = open(filename, "r")

z = []
for i in range(n):
	z.append([])
	r = infile.readline()
	w = r.split()
	for j in range(n):
		z[i].append(float(w[j]))
Z = matrix(z)



X = linspace(0,1,n)
Y = linspace(0,1,n)
X, Y = meshgrid(X, Y)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,  linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
ax.set_zlabel("u(x,y)")
plt.axis([0,1,0,1])
plt.show()



