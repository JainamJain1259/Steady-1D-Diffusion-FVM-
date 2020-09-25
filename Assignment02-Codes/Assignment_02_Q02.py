
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

num_mesh = 11      # number of mesh points
k=1000             #conductivity

del_x=1.0/((num_mesh-1.0)*2)     # mesh size (Delta_x)(here the length of rod is givem to be 0.5m)
xmesh    = np.zeros(num_mesh)    # Mesh points

#defining specific coefficients
a_p = np.zeros(num_mesh)         #coefficient of the pth node 
a_e = np.zeros(num_mesh)         #coefficient of east side of pth node
a_w = np.zeros(num_mesh)         #coefficient of west side of pth node

sol   = np.zeros(num_mesh)        # solution at the previous time step

# RHs of the discretized linear algebraic system
rhs    = np.zeros(num_mesh)

for i in range(0, len(xmesh)):
    xmesh[i] = i * del_x

a_w[0] = 0.0
a_e[0] = 0.0
a_p[0] = 1.0
rhs[0] = 100
rhs[num_mesh-1] = 500
a_w[num_mesh-1] = 0.0
a_e[num_mesh-1] = 0.0
a_p[num_mesh-1] = 1.0
for i in range(1, num_mesh-1):
    a_p[i] = ((-2.0)*k)/del_x
    a_w[i] = ((1.0)*k)/del_x
    a_e[i] = ((1.0)*k)/del_x


#Exact solution
exact = np.zeros(num_mesh)
for i in range(0,num_mesh):
  exact[i] = 800*(xmesh[i])+100


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



#======================plot the converged results=========================#
plt.plot(xmesh,sol,'r-o',marker='v')
plt.plot(xmesh,exact,'g-o',marker='^')
plt.xlabel('Mesh Points')
plt.ylabel('Temperature')
plt.show()



