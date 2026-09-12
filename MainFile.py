from __future__ import division
import numpy as np
import os
from shutil import copy, rmtree

#Output folders: create them if missing, empty them if they already exist
for folder in ['./outputsVertices', './outputsJunctions', './outputsT1']:
    if os.path.exists(folder):
        rmtree(folder)
    os.makedirs(folder)

#Functions and classes: vertices-cells-tissue
exec(open("classes_and_functions.py").read())

factor = np.sqrt(2/(3*np.sqrt(3)))
mu  = 0.2
Gamma_act = 0.03
tau_fluct = 0.4
sigma_fluct = 0.02
difusion = (sigma_fluct**2)/tau_fluct
tau = 10

KA = 1
KP = 0
L_x = 19* np.sqrt(3)*factor
L_y = 39 *factor
L_z = 0
l_T1 = 0.05
u0 = 0
Kl = 1
Kc = 0.1
Ke = 0.2
e_critico = 0.1
T0 = 0.1

dt = 0.004
d_pasos = 300
d_pasos_T1 = 10
total_txt = 120


#Initial configuration of vertices and cells
celdas = np.loadtxt('./config/celda.txt')
data = np.loadtxt('./config/vertices.txt')
area0_ec_mov = np.loadtxt('./config/celulas.txt')
vtipo1 = np.loadtxt('./config/topology.txt')[:,0].astype(int)
top_vtipo1 = np.loadtxt('./config/topology.txt', dtype =float).astype(int)

celda = celdas[:,0]
v_celda = celdas[:,1]
xx = data[:,0]
yy = data[:,1]
tipo_v = np.ones(len(data))
adys_tipo1 = list(zip(top_vtipo1[:,1].astype(int),top_vtipo1[:,2].astype(int),top_vtipo1[:,3].astype(int),
                      top_vtipo1[:,4].astype(int),top_vtipo1[:,5].astype(int),top_vtipo1[:,6].astype(int)))
cells_tipo1 =list(zip(top_vtipo1[:,7].astype(int),top_vtipo1[:,8].astype(int),top_vtipo1[:,9].astype(int)))


x = []
y = []
for i in range(len(xx)):
    r = [xx[i]*1, yy[i]*1]
    r_bien = resta_inicial(r)
    x.append(r_bien[0])
    y.append(r_bien[1])

n_celulas = int(celda[len(celda)-1]+1)
n_vertices = int(len(x))


#Creation of the tissue
exec(open("tissue_creation.py").read())

for v in T1.vs:
    resta_inicial(v.r)


T1.cal_area()
T1.data()
copy('./outputsVertices/data_vertices.txt','./outputsVertices/0_vertices.txt')

lista_i = []
lista_j = []
flats = []
for i in range(n_vertices):
    flat = [T1.vs[i].ady[0][0],T1.vs[i].ady[0][1],T1.vs[i].ady[1][0]]
    flats.append(flat)
for i in range(n_vertices):
    for j in flats[i]:
        if j>i and j in np.array(T1.vs[i].ady).flatten():
            lista_i.append(i)
            lista_j.append(j)

np.savetxt('./outputsJunctions/0.txt', np.c_[lista_i, lista_j], fmt='%1.5f')

for i in range(int(d_pasos * total_txt)-1):
    T1.evol_vertex((i+1)*dt)
    c = (i + 1) % d_pasos
    c2 = (i + 1) % d_pasos_T1
    if c2 == 0 :
        T1.try_T1_tejido((i+1)*dt)
    d = int((i + 1) / d_pasos)
    if c == 0 :
        T1.data()
        copy('./outputsVertices/data_vertices.txt', './outputsVertices/'+ str(d) + '_vertices.txt')

        n_i3 = []
        n_i4= []
        lista_i = []
        lista_j = []
        flats3 = []
        flats4 = []
        for i in range(n_vertices):
            if len(T1.vs[i].cells)==3:
                n_i3.append(i)
                flat = [T1.vs[i].ady[0][0],T1.vs[i].ady[0][1],T1.vs[i].ady[1][0]]
                flats3.append(flat)
            if len(T1.vs[i].cells)==4:
                n_i4.append(i)
                flat = [T1.vs[i].ady[0][0],T1.vs[i].ady[0][1],T1.vs[i].ady[1][0],T1.vs[i].ady[2][0]]
                flats4.append(flat)
        for ii in range(len(n_i3)):
            i = n_i3[ii]
            for j in flats3[ii]:
                if j>i and len(T1.vs[j].cells)==3:
                    lista_i.append(i)
                    lista_j.append(j)
        for ii in range(len(n_i4)):
            i = n_i4[ii]
            for j in flats4[ii]:
                lista_i.append(i)
                lista_j.append(j)


        np.savetxt('./outputsJunctions/'+str(d)+'.txt', np.c_[lista_i, lista_j], fmt='%1.5f')


