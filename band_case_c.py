import numpy as np
from numpy.linalg import eig
# *1jmport sympy as sp
# from sympy*1jmport DiracDelta

# import math
import matplotlib.pyplot as plt
# *1jmport time
import numba
from numba import jit

i = 0
Nx = 40
Ny = 500
Nbandw = 8
muw = 0.0

xpos = np.arange(Nx)


@jit()
def DiracDelta(x, y):
    if y == -x:
        return 1
    elif y == x:
        return 0
    elif y != x:
        return 0


@jit()
def gso(rxsp, kys):
    return np.array([[(-0.163566-0.00319154*1j)*DiracDelta(1.,-rxsp)+0.295217*DiracDelta(1,rxsp-1)-(0.163566-0.00319154*1j)*DiracDelta(1.,+rxsp),0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,-rxsp)-
(0.+0.00478731*1j) *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)-(0.+0.00478731*1j) *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1)-0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,+rxsp),
0.203461 *np.exp((0.-0.36*1j)*kys)*DiracDelta(0.5,-rxsp)+0.203461 *np.exp((0.-0.36*1j)*kys)*DiracDelta(0.5,+rxsp),0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+
0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,-rxsp)-0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp)-0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,+rxsp),
-0.0123672*DiracDelta(1.,-rxsp)+0.0123672*DiracDelta(1.,+rxsp),-0.0203461 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)-0.0199471 *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1),
0,-0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp)],[-0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,-rxsp)+
(0.+0.00478731*1j) *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)+(0.+0.00478731*1j) *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1)+0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,+rxsp),
(0.450805-0.00398942*1j)*DiracDelta(1.,-rxsp)-0.698149*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.-1.*1j)*kys)*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.+1.*1j)*kys)*DiracDelta(1,rxsp-1)+
(0.450805+0.00398942*1j)*DiracDelta(1.,+rxsp),-0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,-rxsp)+
0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp)+0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,+rxsp),0.159577 *np.exp((0.-0.14*1j)*kys)*DiracDelta(0.5,-rxsp)+
0.159577 *np.exp((0.-0.14*1j)*kys)*DiracDelta(0.5,+rxsp),0.0203461 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)+0.0199471 *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1),
-0.0159577*DiracDelta(1.,-rxsp)+0.0159577*DiracDelta(1.,+rxsp),-0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp),0]
,[0.203461 *np.exp((0.+0.36*1j)*kys)*DiracDelta(0.5,-rxsp)+0.203461 *np.exp((0.+0.36*1j)*kys)*DiracDelta(0.5,+rxsp),0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+
0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,-rxsp)-0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp)-0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,+rxsp),
(-0.163566+0.00319154*1j)*DiracDelta(1.,-rxsp)+0.295217*DiracDelta(1,rxsp-1)-(0.163566+0.00319154*1j)*DiracDelta(1.,+rxsp),0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,-rxsp)
+(0.+0.00478731*1j) *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)+(0.+0.00478731*1j) *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1)-0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,+rxsp)
,0,0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp),0.0123672*DiracDelta(1.,-rxsp)-0.0123672*DiracDelta(1.,+rxsp),
0.0203461 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)+0.0199471 *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1)],[-0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-
0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,-rxsp)+0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp)+0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,+rxsp),
0.159577 *np.exp((0.+0.14*1j)*kys)*DiracDelta(0.5,-rxsp)+0.159577 *np.exp((0.+0.14*1j)*kys)*DiracDelta(0.5,+rxsp),-0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,-rxsp)-
(0.+0.00478731*1j) *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)-(0.+0.00478731*1j) *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1)+0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,+rxsp),
(0.450805+0.00398942*1j)*DiracDelta(1.,-rxsp)-0.698149*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.-1.*1j)*kys)*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.+1.*1j)*kys)*DiracDelta(1,rxsp-1)+
(0.450805-0.00398942*1j)*DiracDelta(1.,+rxsp),0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp),0,
-0.0203461 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)-0.0199471 *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1),0.0159577*DiracDelta(1.,-rxsp)-0.0159577*DiracDelta(1.,+rxsp)],
[0.0123672*DiracDelta(1.,-rxsp)-0.0123672*DiracDelta(1.,+rxsp),0.0203461 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)+0.0199471 *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1),
0,0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp),(-0.163566+0.00319154*1j)*DiracDelta(1.,-rxsp)+
0.295217*DiracDelta(1,rxsp-1)-(0.163566+0.00319154*1j)*DiracDelta(1.,+rxsp),0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,-rxsp)+(0.+0.00478731*1j) *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)+
(0.+0.00478731*1j) *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1)-0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,+rxsp),0.203461 *np.exp((0.-0.36*1j)*kys)*DiracDelta(0.5,-rxsp)+
0.203461 *np.exp((0.-0.36*1j)*kys)*DiracDelta(0.5,+rxsp),0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,-rxsp)-
0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp)-0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,+rxsp)],[-0.0203461 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)-
0.0199471 *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1),0.0159577*DiracDelta(1.,-rxsp)-0.0159577*DiracDelta(1.,+rxsp),0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+
0.00438837 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp),0,-0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,-rxsp)-(0.+0.00478731*1j) *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)-
(0.+0.00478731*1j) *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1)+0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,+rxsp),(0.450805+0.00398942*1j)*DiracDelta(1.,-rxsp)-
0.698149*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.-1.*1j)*kys)*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.+1.*1j)*kys)*DiracDelta(1,rxsp-1)+(0.450805-0.00398942*1j)*DiracDelta(1.,+rxsp),
-0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,-rxsp)+0.155587 *np.exp((0.+0.25*1j)*kys)*DiracDelta(0.5,+rxsp)+
0.115693 *np.exp((0.+0.25*1j)*kys)*DiracDelta(1.5,+rxsp),0.159577 *np.exp((0.-0.14*1j)*kys)*DiracDelta(0.5,-rxsp)+0.159577 *np.exp((0.-0.14*1j)*kys)*DiracDelta(0.5,+rxsp)],
[0,-0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp),-0.0123672*DiracDelta(1.,-rxsp)+
0.0123672*DiracDelta(1.,+rxsp),-0.0203461 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)-0.0199471 *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1),0.203461 *np.exp((0.+0.36*1j)*kys)*DiracDelta(0.5,-rxsp)
+0.203461 *np.exp((0.+0.36*1j)*kys)*DiracDelta(0.5,+rxsp),0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)+0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,-rxsp)-
0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp)-0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,+rxsp),(-0.163566-0.00319154*1j)*DiracDelta(1.,-rxsp)+
0.295217*DiracDelta(1,rxsp-1)-(0.163566-0.00319154*1j)*DiracDelta(1.,+rxsp),0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,-rxsp)-(0.+0.00478731*1j) *np.exp((0.-0.39*1j)*kys)*DiracDelta(1,rxsp-1)-
(0.+0.00478731*1j) *np.exp((0.+0.61*1j)*kys)*DiracDelta(1,rxsp-1)-0.0558519 *np.exp((0.-0.39*1j)*kys)*DiracDelta(1.,+rxsp)],[-0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-
0.00438837 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp),0,0.0203461 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)+0.0199471 *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1),
-0.0159577*DiracDelta(1.,-rxsp)+0.0159577*DiracDelta(1.,+rxsp),-0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,-rxsp)-0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,-rxsp)+
0.155587 *np.exp((0.-0.25*1j)*kys)*DiracDelta(0.5,+rxsp)+0.115693 *np.exp((0.-0.25*1j)*kys)*DiracDelta(1.5,+rxsp),0.159577 *np.exp((0.+0.14*1j)*kys)*DiracDelta(0.5,-rxsp)+
0.159577 *np.exp((0.+0.14*1j)*kys)*DiracDelta(0.5,+rxsp),-0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,-rxsp)+(0.+0.00478731*1j) *np.exp((0.+0.39*1j)*kys)*DiracDelta(1,rxsp-1)+
(0.+0.00478731*1j) *np.exp((0.-0.61*1j)*kys)*DiracDelta(1,rxsp-1)+0.0558519 *np.exp((0.+0.39*1j)*kys)*DiracDelta(1.,+rxsp),(0.450805-0.00398942*1j)*DiracDelta(1.,-rxsp)-
0.698149*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.-1.*1j)*kys)*DiracDelta(1,rxsp-1)+0.0518625 *np.exp((0.+1.*1j)*kys)*DiracDelta(1,rxsp-1)+(0.450805+0.00398942*1j)*DiracDelta(1.,+rxsp)]])


# print(gso(0, 0)[0, 0])

cp = np.diag(np.array([-muw for i in range(Nx)]))

r = np.array([[0., 0.], [0., -0.39], [0.5, -0.64], [0.5, -0.25], [0., 0.], [0., -0.39], [0.5, -0.64], [0.5, -0.25]])

# tw = np.zeros((Nx,5,8,8))

# tijmn = np.zeros((Nx*Nbandw,Nx*Nbandw))


eigenval = []

@jit()
def t(x, jy, m, n):
    if jy > Ny - 1 or jy < 0 or abs(x) > 2 or m > 7 or m < 0 or n > 7 or n < 0:
        return None
    return np.sqrt(2 * np.pi) * gso(x + r[m][0] - r[n][0], (2 * jy / Ny - 1) * np.pi)[m, n]
    # print(t(1,-3,0,1))

@jit()  
def tijmn(jy,lx, ly):
    if lx <= Nx * Nbandw / 2:
        m = round(((lx - 1 / 10) % int(Nbandw / 2)) + 1 / 10)
        i = round(((lx - 1 / 10) // int(Nbandw / 2)) + 1)
    else:
        m = round(((lx - 1 / 10) % int(Nbandw / 2)) + 1 / 10 + 4)
        i = round(((lx - 1 / 10) // int(Nbandw / 2)) + 1 - Nx)
    if ly <= Nx * Nbandw / 2:
        n = round(((ly - 1 / 10) % int(Nbandw / 2)) + 1 / 10)
        j = round(((ly - 1 / 10) // int(Nbandw / 2)) + 1)
    else:
        n = round(((ly - 1 / 10) % int(Nbandw / 2)) + 1 / 10 + 4)
        j = round(((ly - 1 / 10) // int(Nbandw / 2)) + 1 - Nx)
    # print(jx,i,j,(ypos[i-1]-ypos[j-1]),m,n, t(jx,(ypos[i-1]-ypos[j-1]),m-1,n-1))
    # return t(jx,(ypos[i-1]-ypos[j-1]),m-1,n-1)

    if t(xpos[i - 1] - xpos[j - 1], jy, m - 1, n - 1) is not None:
        return t(xpos[i - 1] - xpos[j - 1], jy, m - 1, n - 1)
    #elif t(xpos[i-1]-xpos[j-1]-Nx,jy,m-1,n-1) is not None:
      #  return t(xpos[i-1]-xpos[j-1]-Nx,jy,m-1,n-1)
    #elif t((xpos[i-1]-xpos[j-1]+Nx),jy,m-1,n-1) is not None:
     #   return t(xpos[i-1]-xpos[j-1]+Nx,jy,m-1,n-1)
for jy in range(Ny):
    a = np.array([[tijmn(jy,lx + 1, ly + 1) for lx in range(Nx * Nbandw)] for ly in range(Nx * Nbandw)],
                 dtype=np.dtype('complex128'))
    a[np.isnan(a)] = 0
    # print(a)

    v, w = eig(a)
    eigenval.append(np.real(v))
    print(jy)

eval_list = np.array(eigenval)
# print(eval_list)
# %%
###########################################################################
# plotting data

plot_data = []
for j in range(Ny):
    for i in range(len(eval_list[0])):
        a = [(2 * j / Ny - 1), eval_list[j, i]]
        plot_data.append(a)
# print(eval_list)
# print(plot_data)

# print(plot_data.flatten())
x = np.array([plot_data[i][0] for i in range(len(plot_data))])
y = np.array([plot_data[i][1] for i in range(len(plot_data))])

plt.rcParams["figure.figsize"] = [5, 5]
plt.rcParams["figure.autolayout"] = True
plt.plot(x, y, 'bo', markersize=0.5)
plt.ylim(-0.4, 0.1)
# plt.legend(loc=(0.62,0.72),frameon=False,fontsize=15)
plt.xlabel('$ak_x$', fontsize=15)
plt.ylabel('$\omega$ (eV)', fontsize=15)
plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                top=True, left=True, right=True)
plt.yticks(fontsize=15)
plt.xticks((-1, -0.5, 0, 0.5, 1), fontsize=15)
plt.savefig('case_c_surface_band.png', dpi=600)
plt.show()
