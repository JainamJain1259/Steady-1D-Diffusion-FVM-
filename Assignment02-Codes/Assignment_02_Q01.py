

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

num_mesh = 11                   # number of mesh points
k=0.5                           #conductivity
L= 0.02                         #Length of the rod
q = 1000000                     #Source Term
del_x=L/((num_mesh-1.0))        # mesh size (Delta_x)
xmesh    = np.zeros(num_mesh)   # Mesh points

#defining specific coefficients
a_p = np.zeros(num_mesh)        #coefficient of the pth node
a_e = np.zeros(num_mesh)        #coefficient of east side of pth node
a_w = np.zeros(num_mesh)        #coefficient of west side of pth node

# Diagonal elements of system matrix
d    = np.zeros(num_mesh)        # main diagonal elements
u    = np.zeros(num_mesh)        # upper diagonal
l    = np.zeros(num_mesh)        # lower diagonal

sol   = np.zeros(num_mesh)        # solution at the previous time step

# RHs of the discretized linear algebraic system
rhs    = np.zeros(num_mesh)

# RHs of the discretized linear algebraic system
f    = np.zeros(num_mesh)# Compute the location of mesh points

for i in range(0, len(xmesh)):
    xmesh[i] = i * del_x
print(xmesh)

#applying boundary conditions#

T_0 = 100
T_l = 200

a_w[0] = 0.0
a_e[0] = 0.0
a_p[0] = 1.0
rhs[0] = 100
rhs[num_mesh-1] = 200
a_w[num_mesh-1] = 0.0
a_e[num_mesh-1] = 0.0
a_p[num_mesh-1] = 1.0
for i in range(1, num_mesh-1):
    a_p[i] = ((2.0)*k)/del_x
    a_w[i] = ((-1.0)*k)/del_x
    a_e[i] = ((-1.0)*k)/del_x
    rhs[i] = (q*del_x)   

#Exact solution
exact = np.zeros(num_mesh)
for i in range(0,num_mesh):
  exact[i] = ((T_l-T_0)/L + (q*(L-xmesh[i])/(2*k)))*xmesh[i] + T_0


#==============================THOMAS_ALGORITHM==========================#
dia1    = np.zeros(num_mesh)
rhs1    = np.zeros(num_mesh)
dia1[0] = a_p[0]
rhs1[0] = rhs[0]
for i in range(1,num_mesh-1):
        dia1[i] = a_p[i] - a_w[i]*a_e[i-1]/dia1[i-1]
        rhs1[i] = rhs[i] - rhs1[i-1]*a_w[i]/dia1[i-1]
       # print(dia1[i], u[i],l[i],rhs1[i])
sol[num_mesh-1] = rhs[num_mesh-1]/a_p[num_mesh-1]
for i in range(len(a_p)-2,-1,-1):
        sol[i] = (rhs1[i]-a_e[i]*sol[i+1])/dia1[i]

#=======================plot the converged results=========================#
plt.plot(xmesh,sol,'r-o',marker='v')
plt.plot(xmesh,exact,'g-o',marker='^')
plt.xlabel('Mesh Points')
plt.ylabel('Temperature')
plt.show()



