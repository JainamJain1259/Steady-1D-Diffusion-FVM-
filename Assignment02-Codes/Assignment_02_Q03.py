

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

num_mesh = 21                  #number of mesh points    
L= 0.6                          #Length of the rod
q = 1000000                     #Source Term
del_x=L/((num_mesh-1.0))        # mesh size (Delta_x)
xmesh    = np.zeros(num_mesh)   # Mesh points
start = 0.0                     #start of the rod
end = 0.6                       #end of the rod
interface_1 = 0.3               #Position of 1st interface(end of first composite rod)
interface_2 = 0.45              #Position of 2nd interface(end of 2nd composite rod)
temperature_end = 293           #Temperature at the left end
temperature_fluid = 1073        #Temperature of the fluid at left end 
h_coeff = 25                    #fluid heat transfer coefficient value
source = 0                      #source term


conductivity_a = np.zeros(num_mesh)   #Arithmetic mean of conductivity
conductivity_h = np.zeros(num_mesh)   #Harmonic mean of conductivity
conductivity_1 = 20                   #Condutivity in first sub-domain
conductivity_2 = 1.5                  #Conducitivty in second sub-domain
conductivity_3 = 50                   #Conductivity in third sub-domain

sol_A = np.zeros(num_mesh)            #Solution of Arithmetic mean
sol_H = np.zeros(num_mesh)            #Solution of Harmonic mean

#===================Defining matrix elements============================#
dia_A= np.zeros(num_mesh)             #arithmetic coefficient of the pth node
low_A = np.zeros(num_mesh)            #arithmetic coefficient of east side of pth node
upp_A = np.zeros(num_mesh)            #arithmetic coefficient of west side of pth node
func_A = np.zeros(num_mesh)           #right hand side elements, for arithmetic mean

dia_H= np.zeros(num_mesh)             #harmonic coefficient of the pth node
low_H = np.zeros(num_mesh)            #harmonic coefficient of east side of pth node
upp_H = np.zeros(num_mesh)            #harmonic coefficient of west side of pth node
func_H = np.zeros(num_mesh)           #right hand side elements, for harmonic mean

Temp_exact = np.zeros(num_mesh)       #Exact solution of the given problem

#initialising an array of mesh points to locate each node
for i in range(0, len(xmesh)):
    xmesh[i] = i * del_x

#Idea: I have took a conductivity array and updated the value in the array at each and every node. 
#If there is an interface between any two nodes then Arithmetic as well as Harmonic mean is calculated 
#If there is no interface between two nodes then the constant conductivity value of respective sub-domain is feeded.
for i in range(0,num_mesh):
  if(xmesh[i]<=interface_1):
    if(xmesh[i+1]>interface_1):
      fe = (xmesh[i+1]-interface_1)/del_x
      conductivity_a[i] = fe*(conductivity_2)+ (1-fe)*conductivity_1 #Arithmetic interpolation at first interface
      conductivity_h[i] = 1/((fe/conductivity_2)+((1-fe)/(conductivity_1))) #Harmonic interpolation at first interface
    else:
      conductivity_a[i] = conductivity_1
      conductivity_h[i]  = conductivity_1
  elif(xmesh[i]<=interface_2 and xmesh[i]>interface_1):
      if(xmesh[i+1]>interface_2):
       fe = (xmesh[i+1]-interface_2)/del_x 
       conductivity_a[i] = fe*(conductivity_3)+ (1-fe)*conductivity_2 #Arithmetic interpolation at first interface
       conductivity_h[i] = 1/((fe/conductivity_3)+((1-fe)/(conductivity_2))) #Harmonic interpolation at first interface
      else:
        conductivity_a[i] = conductivity_2
        conductivity_h[i] = conductivity_2
  elif(xmesh[i]>interface_2 and xmesh[i]<= end):
      conductivity_a[i] = conductivity_3
      conductivity_h[i] = conductivity_3

#-----------------------------EXACT VALUE CALCULATION--------------------------------------------------
for i in range(0, num_mesh):
    if (xmesh[i]<=interface_1):
        m1 = ((801.48-875.53)/(interface_1-start))  #exact solution in the first sub-domain
        Temp_exact[i] = 875.53 + (m1*xmesh[i])
    elif (xmesh[i]<=interface_2) and (xmesh[i]>interface_1):
        m2 = ((307.81-801.48)/(interface_2-interface_1))
        Temp_exact[i] = 801.48 + (m2*(xmesh[i] - interface_1))    #exact solution in the second sub-domain
    elif (xmesh[i]<=end) and (xmesh[i]>interface_2):
        m3 = ((293-307.81)/(end-interface_2))
        Temp_exact[i] = 307.81 + (m3*(xmesh[i] - interface_2))    #exact solution in the thrid sub-domain
#================Diagonal Entries============#
for i in range(1,num_mesh-1):
  low_A[i] = ((-1.0)*conductivity_a[i-1])/del_x
  upp_A[i] = ((-1.0)*conductivity_a[i])/del_x
  dia_A[i] = (conductivity_a[i]/del_x) + (conductivity_a[i-1]/del_x)
  func_A[i] = source*del_x
  low_H[i] = ((-1.0)*conductivity_h[i-1])/del_x
  upp_H[i] = ((-1.0)*conductivity_h[i])/del_x
  dia_H[i] = (conductivity_h[i-1]/del_x) + (conductivity_h[i]/del_x)
  func_H[i] = (source*del_x)/2

#==================Boundary Condition=================#

low_A[0] = 0.0
upp_A[0] = 0.0
dia_A[0] = 1.0                                     #Boundary conditions for left hand side(Arithmetic)
func_A[0] = 875.531                                #This temperature is obtained by equating the fluxes at interface 1,2,3 and 4


low_H[0] = 0.0
upp_H[0] = 0.0
dia_H[0] = 1.0                                     #Boundary conditions for left Hand Side(Harmonic)
func_H[0] = 875.531                                #This temperature is obtained by equating the fluxes at interface 1,2,3 and 4


low_A[num_mesh-1] = 0.0
upp_A[num_mesh-1] = 0.0
dia_A[num_mesh-1] = 1.0                            #Boundary conditions for right hand side(Arithmetic)
func_A[num_mesh-1] = temperature_end


low_H[num_mesh-1] = 0.0
upp_H[num_mesh-1] = 0.0
dia_H[num_mesh-1] = 1.0                            #Boundary conditions for right hand side(Harmonic)
func_H[num_mesh-1] = temperature_end


#==============================THOMAS_ALGORITHM(For Arithmetic)==========================#
dia1    = np.zeros(num_mesh)
rhs1    = np.zeros(num_mesh)
dia1[0] = dia_A[0]
rhs1[0] = func_A[0]
for i in range(1,num_mesh-1):
        dia1[i] = dia_A[i] - low_A[i]*upp_A[i-1]/dia1[i-1]
        rhs1[i] = func_A[i] - rhs1[i-1]*low_A[i]/dia1[i-1]
sol_A[num_mesh-1] = func_A[num_mesh-1]/dia_A[num_mesh-1]
for i in range(len(dia_A)-2,-1,-1):
        sol_A[i] = (rhs1[i]-upp_A[i]*sol_A[i+1])/dia1[i]

#==============================THOMAS_ALGORITHM(For Arithmetic)==========================#
dia2    = np.zeros(num_mesh)
rhs2    = np.zeros(num_mesh)
dia2[0] = dia_H[0]
rhs2[0] = func_H[0]
for i in range(1,num_mesh-1):
        dia2[i] = dia_H[i] - low_H[i]*upp_H[i-1]/dia2[i-1]
        rhs2[i] = func_H[i] - rhs2[i-1]*low_H[i]/dia2[i-1]
       # print(dia1[i], u[i],l[i],rhs1[i])
sol_H[num_mesh-1] = func_H[num_mesh-1]/dia_H[num_mesh-1]
for i in range(len(dia_H)-2,-1,-1):
        sol_H[i] = (rhs2[i]-upp_H[i]*sol_H[i+1])/dia2[i]

#=======================plot the converged results=========================#
plt.plot(xmesh,sol_A,'r-o',marker='v')
plt.plot(xmesh,Temp_exact,'g-o',marker='^')
plt.plot(xmesh,sol_H,'b-o',marker='*')
plt.xlabel('Mesh Points')
plt.ylabel('Temperature')
plt.show()
