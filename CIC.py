
# coding: utf-8

##ojo esta hecho con 100 x100x100 enves de 1000x1000x1000, 
##ademas va a salir un error de complex casting, ignorenlo es -
##por tratar los complejos dela matriz transformada, a los cuales - 
##ya se les habia sacado la norma porque estamos tratando con un campo vectorial
import numpy as np
import matplotlib.pylab as plt


# In[2]:

n = 100
#espacio de densidades:
grid = [[[0 for x in xrange(n)] for y in xrange(n)] for z in xrange(n)]


# In[3]:

m = 1
x = np.genfromtxt('Serena-Venus.txt',delimiter=" ", usecols=(1))
y = np.genfromtxt('Serena-Venus.txt',delimiter=" ", usecols=(2))
z = np.genfromtxt('Serena-Venus.txt',delimiter=" ", usecols=(3))
##

x = x[0:100]
y = y[0:100]
z = z[0:100]

x = x - x.min()
y = y - y.min()
z = z - z.min()


# In[4]:

X = (x.max() - x.min())/n
Y = (y.max() - y.min())/n
Z = (z.max() - z.min())/n
xb = [0]
yb = [0]
zb = [0]
counterx = 0
countery = 0
counterz = 0
for i in range(99):
    counterx += X
    xb.append(counterx)
for j in range(99):
    countery += Y
    yb.append(countery)
for k in range(99):
    counterz += Z
    zb.append(counterz)


# In[5]:

#i,j,k sera el centro de las celdas :D y xp,yp,zp la posicion de la particula :D :D
for i in range(99-1):
    for j in range(99-1):
        for k in range(99-1):
            dx = x[i] - (xb[i]-xb[i-1])
            dy = y[j] - (yb[j]-yb[j-1])
            dz = z[k] - (zb[k]-zb[k-1])
            xd = 1 - dx
            yd = 1 - dy
            zd = 1 - dz
            grid[i][j][k] = np.abs(grid[i][j][k] + m*xd*yd*zd)
            grid[i+1][j][k] = np.abs(grid[i+1][j][k] + m*dx*yd*zd)     
            grid[i][j+1][k] = np.abs(grid[i][j+1][k] + m*xd*dy*zd)   
            grid[i+1][j+1][k] = np.abs(grid[i+1][j+1][k] + m*dx*dy*zd)  
            grid[i][j][k+1] = np.abs(grid[i][j][k+1] + m*xd*yd*dz)  
            grid[i+1][j][k+1] = np.abs(grid[i+1][j][k+1] + m*dx*yd*dz) 
            grid[i][j+1][k+1] = np.abs(grid[i][j+1][k+1] + m*xd*dy*dz)
            grid[i+1][j+1][k+1] = np.abs(grid[i+1][j+1][k+1] + m*dx*dy*dz)


# In[6]:

dgo = -1 * np.fft.fftn(grid)


# In[7]:

def G(x,y,z):
    return (-3/(8*((np.sin(2*np.pi*x/xb[1]))**2+(np.sin(2*np.pi*y/yb[1]))**2+(np.sin(2*np.pi*z/zb[1]))**2)))


# In[8]:

phigo = dgo
for i in range(99-1):
    for j in range(99-1):
        for k in range(99-1):
            if i == 0 and j == 0 and k == 0:
                phigo[i][j][k] = 0
            else:
                phigo[i][j][k] = G(i,j,k)*dgo[i][j][k]


# In[9]:

phi = np.fft.ifftn(phigo)


# In[10]:

for i in range(len(phi[:][0][0])):
    for j in range(len(phi[0][:][0])):
        for k in range(len(phi[0][0][:])):
            phi[i][j][k] = np.sqrt(phi[i][j][k].real**2 + phi[i][j][k].imag**2)




# In[12]:

fgrav = np.gradient(phi,X,Y,Z)


# In[13]:

fgrav = fgrav[2][:][:][:] # 1 layer problemas con la funcion, pero esta hace lo que queremos, el gradiente



# In[16]:

mini = np.abs(np.amin(fgrav))

for i in range(len(fgrav[:][0][0])):
    for j in range(len(fgrav[0][:][0])):
        for k in range(len(fgrav[0][0][:])):
            fgrav[i][j][k] = fgrav[i][j][k] + mini


# In[17]:

maxi = np.amax(fgrav)
for i in range(len(fgrav[:][0][0])):
    for j in range(len(fgrav[0][:][0])):
        for k in range(len(fgrav[0][0][:])):
            fgrav[i][j][k] = fgrav[i][j][k]/maxi


# In[18]:

for i in range(len(fgrav[:][0][0])):
    for j in range(len(fgrav[0][:][0])):
        for k in range(len(fgrav[0][0][:])):
            fgrav[i][j][k] = float(fgrav[i][j][k])


# In[19]:

max = []
min = []
# se normalizo los valores discretos del campo y se definio un umbral , para identificar los maximos y minimos.
for i in range(len(fgrav[:][0][0])):
    for j in range(len(fgrav[0][:][0])):
        for k in range(len(fgrav[0][0][:])):
            if fgrav[i][j][k] > 0.97:
                a = [i,j,k]
                max.append(a)
            elif fgrav[i][j][k] < 0.03:
                a = [i,j,k]
                min.append(a)
            else:
                pass
                


# In[32]:

mesh = fgrav[:][:][0]
X = np.genfromtxt('Serena-Venus.txt',delimiter=" ", usecols=(1))
Y = np.genfromtxt('Serena-Venus.txt',delimiter=" ", usecols=(2))
Z = np.genfromtxt('Serena-Venus.txt',delimiter=" ", usecols=(3))


# In[47]:

plt.figure(figsize=(15,10))
plt.contourf(mesh,cmap=plt.get_cmap('CMRmap'))
plt.colorbar()
plt.scatter(X,Y)
plt.xlim(33.8,36.5)
plt.ylim(33.7,36.8)
plt.title('Proyeccion en una cara de un cubo conteniendo el campo de fuerza gravitacional')
plt.savefig("Grafica1")
# un valor mas cercano a 1 indica un minimo de la forma contraria para valores cercanos a 0.


# In[48]:

from mpl_toolkits.mplot3d.axes3d import Axes3D
fig = plt.figure(figsize=(16,10))
ax = fig.add_subplot(1,1,1, projection='3d')
ax.scatter(X,Y,Z)
plt.savefig("GraficaExtra")






