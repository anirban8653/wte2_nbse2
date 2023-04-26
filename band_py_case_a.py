import numpy as np
from scipy.linalg import eigh
# import sympy as sp
# from sympy import DiracDelta

#import math
import matplotlib.pyplot as plt
# import time
# import numba
from numba import jit


#from mpi4py import MPI 
#comm = MPI.COMM_WORLD
#myRank = comm.Get_rank()
##print("my rank is ", myRank)
#coreN = comm.Get_size()
#if myRank ==0:
#    print("this is the total number of cores", coreN)


i = 0
Nx=200
Ny=40
Nbandw=8
muw=0.0




ypos = np.arange(Ny)

@jit()
def DiracDelta(x,y):
    if y == -x:
        return 1
    elif y == x:
        return 0
    elif y != x:
        return 0
@jit()    
def gso(kxs,rxyp):
    return np.array([[0.74*DiracDelta(0.0,rxyp)-0.82*np.cos(1. *kxs)*DiracDelta(0.0,rxyp)-0.016*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),
             (0.-0.012j)*DiracDelta(0.61,-rxyp)-(0.+0.012j)*DiracDelta(0.39,rxyp)+0.14  *np.exp((0.-1j) *kxs) 
            *DiracDelta(0.39,rxyp)-0.14  *np.exp((0.+1j) *kxs)*DiracDelta(0.39,rxyp),0.51  *np.exp((0.-0.5j) *kxs) 
            *DiracDelta(0.36,-rxyp)+0.51  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.36,-rxyp),0.39  *np.exp((0.-0.5j) *kxs) 
            *DiracDelta(0.25,rxyp)-0.39  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp)+0.29  *np.exp((0.-1.5j) *kxs) 
            *DiracDelta(0.25,rxyp)-0.29  *np.exp((0.+1.5j) *kxs)*DiracDelta(0.25,rxyp),(0.+0.062j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),
             -0.05*DiracDelta(0.61,-rxyp)-0.051*DiracDelta(0.39,rxyp),0,-0.011  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,rxyp)-0.011 
             *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp)],[(0.+0.012j)*DiracDelta(0.39,-rxyp)-0.14  *np.exp((0.-1j) *kxs)*DiracDelta(0.39,-rxyp)+0.14 
                                                           *np.exp((0.+1j) *kxs)*DiracDelta(0.39,-rxyp)+(0.+0.012j)*DiracDelta(0.61,rxyp),0.13*DiracDelta(1.,-rxyp)-1.75 
                                                          *DiracDelta(0.0,rxyp)+(1.13-0.01j)  *np.exp((0.-1j) *kxs)*DiracDelta(0.0,rxyp)+(1.13+0.01j) 
                                                           *np.exp((0.+1j) *kxs)*DiracDelta(0.0,rxyp)+0.13*DiracDelta(1.,rxyp),-0.39  *np.exp((0.-0.5j) *kxs) 
                                                          *DiracDelta(0.25,rxyp)+0.39  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp)-0.29 
                                                           *np.exp((0.-1.5j) *kxs)*DiracDelta(0.25,rxyp)+0.29  *np.exp((0.+1.5j) *kxs)
                                                          *DiracDelta(0.25,rxyp),0.4  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.14,-rxyp)+0.4 
                                                           *np.exp((0.+0.5j) *kxs)*DiracDelta(0.14,-rxyp),0.051*DiracDelta(0.39,-rxyp)+0.05
                                                          *DiracDelta(0.61,rxyp),(0.+0.08j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),-0.011  
                                                           *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,rxyp)-0.011  *np.exp((0.+0.5j) *kxs) 
                                                          *DiracDelta(0.25,rxyp),0],[0.51  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.36,rxyp)+0.51 
                                                                                     *np.exp((0.+0.5j) *kxs)*DiracDelta(0.36,rxyp),0.39  *np.exp((0.-0.5j) *kxs) 
                                                                                    *DiracDelta(0.25,-rxyp)-0.39  *np.exp((0.+0.5j) *kxs) 
                                                                                    *DiracDelta(0.25,-rxyp)+0.29  *np.exp((0.-1.5j) *kxs)
                                                                                    *DiracDelta(0.25,-rxyp)-0.29  *np.exp((0.+1.5j) *kxs)
                                                                                    *DiracDelta(0.25,-rxyp),0.74*DiracDelta(0.0,rxyp)-0.82*np.cos(1. *kxs) 
                                                                                    *DiracDelta(0.0,rxyp)+0.016*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),(0.+0.012j)
                                                                                    *DiracDelta(0.39,-rxyp)+0.14  *np.exp((0.-1j) *kxs)*DiracDelta(0.39,-rxyp)-0.14  
                                                                                     *np.exp((0.+1j) *kxs)*DiracDelta(0.39,-rxyp)+(0.+0.012j)*DiracDelta(0.61,rxyp),0,
                                                                                     0.011  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.011  *np.exp((0.+0.5j) *kxs)
                                                                                    *DiracDelta(0.25,-rxyp),(0.-0.062j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),0.051*DiracDelta(0.39,-rxyp)+0.05
                                                                                    *DiracDelta(0.61,rxyp)],[-0.39  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.39  *np.exp((0.+0.5j) *kxs)
                                                                                                            *DiracDelta(0.25,-rxyp)-0.29  *np.exp((0.-1.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.29  
                                                                                                             *np.exp((0.+1.5j) *kxs)*DiracDelta(0.25,-rxyp),0.4  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.14,rxyp)+0.4
                                                                                                             *np.exp((0.+0.5j) *kxs)*DiracDelta(0.14,rxyp),(0.-0.012j)*DiracDelta(0.61,-rxyp)-(0.+0.012j)
                                                                                                            *DiracDelta(0.39,rxyp)-0.14  *np.exp((0.-1j) *kxs)*DiracDelta(0.39,rxyp)+0.14  *np.exp((0.+1j) *kxs)
                                                                                                            *DiracDelta(0.39,rxyp),0.13*DiracDelta(1.,-rxyp)-1.75*DiracDelta(0.0,rxyp)+(1.13+0.01j)  *np.exp((0.-1j) *kxs)
                                                                                                            *DiracDelta(0.0,rxyp)+(1.13-0.01j)  *np.exp((0.+1j) *kxs)*DiracDelta(0.0,rxyp)+0.13*DiracDelta(1.,rxyp),0.011  
                                                                                                             *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.011  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,-rxyp),0,-0.05
                                                                                                            *DiracDelta(0.61,-rxyp)-0.051*DiracDelta(0.39,rxyp),(0.-0.08j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs)],[(0.-0.062j)
                                                                                                                                                                                                          *DiracDelta(0.0,rxyp)*np.sin(1. *kxs),0.05
                                                                                                                                                                                                          *DiracDelta(0.61,-rxyp)+0.051*DiracDelta(0.39,rxyp),0,0.011 
                                                                                                                                                                                                           *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,rxyp)+0.011 
                                                                                                                                                                                                           *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp),0.74
                                                                                                                                                                                                          *DiracDelta(0.0,rxyp)-0.82*np.cos(1. *kxs)*DiracDelta(0.0,rxyp)+0.016
                                                                                                                                                                                                          *DiracDelta(0.0,rxyp)*np.sin(1. *kxs),(0.+0.012j)*DiracDelta(0.61,-rxyp)+(0.+0.012j)
                                                                                                                                                                                                     *DiracDelta(0.39,rxyp)+0.14  *np.exp((0.-1j) *kxs)*DiracDelta(0.39,rxyp)-0.14 
                                                                                                                                                                                                           *np.exp((0.+1j) *kxs)*DiracDelta(0.39,rxyp),0.51  *np.exp((0.-0.5j) *kxs)
                                                                                                                                                                                                          *DiracDelta(0.36,-rxyp)+0.51  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.36,-rxyp),0.39 
                                                                                                                                                                                                           *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,rxyp)-0.39  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp)+0.29
                                                                                                                                                                                                           *np.exp((0.-1.5j) *kxs)*DiracDelta(0.25,rxyp)-0.29  *np.exp((0.+1.5j) *kxs)*DiracDelta(0.25,rxyp)],[-0.051*DiracDelta(0.39,-rxyp)-0.05*DiracDelta(0.61,rxyp),(0.-0.08j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),0.011  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,rxyp)+0.011  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp),0,(0.-0.012j)*DiracDelta(0.39,-rxyp)-0.14  *np.exp((0.-1j) *kxs)*DiracDelta(0.39,-rxyp)+0.14  *np.exp((0.+1j) *kxs)*DiracDelta(0.39,-rxyp)-(0.+0.012j)*DiracDelta(0.61,rxyp),0.13*DiracDelta(1.,-rxyp)-1.75*DiracDelta(0.0,rxyp)+(1.13+0.01j)  *np.exp((0.-1j) *kxs)*DiracDelta(0.0,rxyp)+(1.13-0.01j)  *np.exp((0.+1j) *kxs)*DiracDelta(0.0,rxyp)+0.13*DiracDelta(1.,rxyp),-0.39  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,rxyp)+0.39  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,rxyp)-0.29  *np.exp((0.-1.5j) *kxs)*DiracDelta(0.25,rxyp)+0.29  *np.exp((0.+1.5j) *kxs)*DiracDelta(0.25,rxyp),0.4  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.14,-rxyp)+0.4  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.14,-rxyp)],[0,-0.011  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)-0.011  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,-rxyp),(0.+0.062j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),-0.051*DiracDelta(0.39,-rxyp)-0.05*DiracDelta(0.61,rxyp),0.51  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.36,rxyp)+0.51  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.36,rxyp),0.39  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)-0.39  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.29  *np.exp((0.-1.5j) *kxs)*DiracDelta(0.25,-rxyp)-0.29  *np.exp((0.+1.5j) *kxs)*DiracDelta(0.25,-rxyp),0.74*DiracDelta(0.0,rxyp)-0.82*np.cos(1. *kxs)*DiracDelta(0.0,rxyp)-0.016*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),(0.-0.012j)*DiracDelta(0.39,-rxyp)+0.14 *np.exp((0.-1j) *kxs)*DiracDelta(0.39,-rxyp)-0.14  *np.exp((0.+1j) *kxs)*DiracDelta(0.39,-rxyp)-(0.+0.012j)*DiracDelta(0.61,rxyp)],[-0.011  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)-0.011  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,-rxyp),0,0.05*DiracDelta(0.61,-rxyp)+0.051*DiracDelta(0.39,rxyp),(0.+0.08j)*DiracDelta(0.0,rxyp)*np.sin(1. *kxs),-0.39  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.39  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.25,-rxyp)-0.29  *np.exp((0.-1.5j) *kxs)*DiracDelta(0.25,-rxyp)+0.29  *np.exp((0.+1.5j) *kxs)*DiracDelta(0.25,-rxyp),0.4  *np.exp((0.-0.5j) *kxs)*DiracDelta(0.14,rxyp)+0.4  *np.exp((0.+0.5j) *kxs)*DiracDelta(0.14,rxyp),(0.+0.012j)*DiracDelta(0.61,-rxyp)+(0.+0.012j)*DiracDelta(0.39,rxyp)-0.14  *np.exp((0.-1j) *kxs)*DiracDelta(0.39,rxyp)+0.14  *np.exp((0.+1j) *kxs)*DiracDelta(0.39,rxyp),0.13*DiracDelta(1.,-rxyp)-1.75*DiracDelta(0.0,rxyp)+(1.13-0.01j)  *np.exp((0.-1j) *kxs)*DiracDelta(0.0,rxyp)+(1.13+0.01j)  *np.exp((0.+1j) *kxs)*DiracDelta(0.0,rxyp)+0.13*DiracDelta(1.,rxyp)]])

# print(gso(1,1)[1,1])

cp = np.diag(np.array([-muw for i in range(Ny)]))

r  = np.array([[0.,0.],[0.,-0.39],[0.5,-0.64],[0.5,-0.25],[0.,0.],[0.,-0.39],[0.5,-0.64],[0.5,-0.25]])



eigenval = []


@jit()
def t(jx,y,m,n):
    if jx>Nx-1 or jx<0 or abs(y) > 2  or m > 7 or m < 0 or n > 7 or n < 0:
        return None
    return gso((2*jx/Nx-1)*np.pi,y-r[m][1]+r[n][1])[m,n]


@jit()
def tijmn(jx,lx,ly):
                if lx <= Ny*Nbandw/2:
                    m = round(((lx-1/10) % int(Nbandw/2)) + 1/10)
                    i = round(((lx-1/10) // int(Nbandw/2)) + 1)
                else:
                    m = round(((lx-1/10) % int(Nbandw/2)) + 1/10 + 4)
                    i = round(((lx-1/10) // int(Nbandw/2)) + 1 - Ny)
                if ly <= Ny*Nbandw/2:
                    n = round(((ly-1/10) % int(Nbandw/2)) + 1/10)
                    j = round(((ly-1/10) // int(Nbandw/2)) + 1)
                else:
                    n = round(((ly-1/10) % int(Nbandw/2)) + 1/10 + 4)
                    j = round(((ly-1/10) // int(Nbandw/2)) + 1 - Ny)
                # print(jx,i,j,(ypos[i-1]-ypos[j-1]),m,n, t(jx,(ypos[i-1]-ypos[j-1]),m-1,n-1))
                # return t(jx,(ypos[i-1]-ypos[j-1]),m-1,n-1)
                
            
                if t(jx,ypos[i-1]-ypos[j-1],m-1,n-1) is not None:
                    return t(jx,ypos[i-1]-ypos[j-1],m-1,n-1)
    #            elif t(jx,ypos[i-1]-ypos[j-1] -Ny,m-1,n-1) is not None:
    #               return t(jx,ypos[i-1]-ypos[j-1]-Ny,m-1,n-1)
    #          elif t(jx,ypos[i-1]-ypos[j-1] + Ny,m-1,n-1) is not None:
    #             return t(jx,ypos[i-1]-ypos[j-1] + Ny,m-1,n-1)

               
    
for jx in range(Nx):
    
    a = np.array([[tijmn(jx,lx+1, ly+1) for lx in range(Ny*Nbandw)] for ly in              range(Ny*Nbandw)],dtype=np.dtype('complex128'))  
    a[np.isnan(a)] = 0
    # print(a)
    
    v,w = eigh(a)
    eigenval.append(np.real(v))
    print(jx)

    
eval_list = np.array(eigenval)

#%%
###########################################################################
# plotting data

plot_data = []
for j in range(Nx):
    for i in range(len(eval_list[0])):
        a = [(2*j/Nx-1)*np.pi,eval_list[j,i]] 
        plot_data.append(a)
# print(eval_list)
# print(plot_data)

#print(plot_data.flatten())
x = np.array([plot_data[i][0] for i in range(len(plot_data))])
y = np.array([plot_data[i][1] for i in range(len(plot_data))])


plt.rcParams["figure.figsize"] = [5, 5]
plt.rcParams["figure.autolayout"] = True
plt.plot(x,y, 'bo', markersize=0.5)
plt.ylim(-0.4,0.4)
#plt.legend(loc=(0.62,0.72),frameon=False,fontsize=15)
plt.xlabel('$ak_x$',fontsize=15)
plt.ylabel('$\omega$ (eV)',fontsize=15)
plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                top=True, left=True, right=True)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig('case_a_surface_band.png', dpi=600)
plt.show()

