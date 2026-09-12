from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from shutil import copy
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib
import random
from matplotlib import rcParams


plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 16

def suma(lista1, lista2):
    lista = []
    for i in range(len(lista1)):
        lista.append(lista1[i]+lista2[i])
    return lista

def resta(lista1, lista2):
    lista = []
    for i in range(len(lista1)):
        lista.append(lista1[i]-lista2[i])
    return lista

def div(lista1, lista2):
    lista = []
    for i in range(len(lista1)):
        lista.append(lista1[i]/lista2[i])
    return lista

def mult(lista1, lista2):
    lista = []
    for i in range(len(lista1)):
        lista.append(lista1[i]*lista2[i])
    return lista
time_SUNITS = (300*120-1)*0.004
texto0 = './outputsJunctions/0.txt'
data0 = np.loadtxt(texto0)


rcParams['font.family'] = 'sans-serif'
rcParams['mathtext.default'] = 'regular'

data4fold = np.loadtxt('./outputsT1/4-fold.txt')
ti_4fold = data4fold[:,0]
vi_4fold = data4fold[:,1]
vj_4fold = data4fold[:,2]
lij_pre = data4fold[:,3]
tension_pre = data4fold[:,4]
gamma_pre = data4fold[:,5]
fluct_pre = data4fold[:,6]

num_dataT1 = sum(1 for line in open('./outputsT1/T1_original.txt'))
dataT1 = np.loadtxt('./outputsT1/T1_original.txt')
if num_dataT1>1:
    T1ori = np.copy(dataT1)
    tf_ori =  T1ori[:,0]
    vi_ori = T1ori[:,1]
    vj_ori = T1ori[:,2]
elif num_dataT1==1:
    T1ori = np.copy(dataT1)
    tf_ori =  [T1ori[0]]
    vi_ori = [T1ori[1]]
    vj_ori = [T1ori[2]]
else:
    tf_ori =  []
    vi_ori = []
    vj_ori = []

num_dataRev = sum(1 for line in open('./outputsT1/T1_perpendicular.txt'))
dataRev = np.loadtxt('./outputsT1/T1_perpendicular.txt')
if num_dataRev>1:
    T1perp = np.copy(dataRev)
    tf_perp =  T1perp[:,0]
    vi_perp = T1perp[:,1]
    vj_perp = T1perp[:,2]
elif num_dataRev==1:
    T1perp = np.copy(dataRev)
    tf_perp =  [T1perp[0]]
    vi_perp = [T1perp[1]]
    vj_perp = [T1perp[2]]
else:
    tf_perp = []
    vi_perp = []
    vj_perp = []

vi = []
vj = []
ti = []
tf = []
tipo = []
tension_ori=[]
tension_perp=[]
for i in range(len(ti_4fold)):
    vi.append(vi_4fold[i])
    vj.append(vj_4fold[i])
    ti.append(ti_4fold[i])
    t1 = ti_4fold[i]
    t2 =  time_SUNITS
    q = 0

    tf_probable_ori=t2
    tf_probable_perp=t2
    if len(tf_ori)>0:
        for j in range(len(tf_ori)):
            if vi_ori[j] == vi[i] and vj_ori[j] == vj[i] and tf_ori[j]>=t1 and tf_ori[j]<t2 and q==0:
                tf_probable_ori = tf_ori[j]
                q=1
            else:
                if q==0:
                    tf_probable_ori=t2

    q=0
    if len(tf_perp)>0:
        for j in range(len(tf_perp)):
            if vi_perp[j] == vi[i] and vj_perp[j] == vj[i] and tf_perp[j]>=t1 and tf_perp[j]<t2 and q==0:
                tf_probable_perp = tf_perp[j]
                q=1
            else:
                if q==0:
                    tf_probable_perp=t2

    if tf_probable_ori==t2 and tf_probable_perp==t2:
        tf.append( time_SUNITS)
        tipo.append(0)
    elif tf_probable_ori<tf_probable_perp:
        tf.append(tf_probable_ori)
        tipo.append(1)
        tension_ori.append(tension_pre[i])
    elif tf_probable_ori>tf_probable_perp:
        tf.append(tf_probable_perp)
        tipo.append(-1)
        tension_perp.append(tension_pre[i])

fig, (ax) = plt.subplots(figsize=(3.8,4))

pparam = dict( xlabel='$t$ [min]', ylabel='# delayed event $\\left[10^3\\right]$')
h = 5*(28/60)
delta_original=[]
delta_perp=[]

delta_unresolved = []
delta_unresolved_vi = []
nni=1
for i in range(len(vi)):
    if tipo[i]==0:
        ax.hlines(y=nni*0.001, xmin=h*ti[i], xmax=h*tf[i], linewidth=.4, color='red', alpha=0.8,zorder=30)
        ax.scatter(h*ti[i], nni*0.001,marker='o', s=.0001,color="gray",zorder=10)
        ax.scatter(h*tf[i], nni*0.001,marker='o', s=3,color="white",edgecolor='red',linewidth=.1)
        delta_unresolved.append(h*tf[i]-h*ti[i])
        delta_unresolved_vi.append(vi[i])
        nni=nni+1
    if tipo[i]==1:
        if h*tf[i]-h*ti[i]<0.84:
            delta_original.append(0)
        else:
            ax.hlines(y=nni*0.001, xmin=h*ti[i], xmax=h*tf[i], linewidth=.2, color='k')
            ax.scatter(h*tf[i], nni*0.001,marker='o', s=.5,color='k', zorder=20)
            ax.scatter(h*ti[i], nni*0.001,marker='o', s=.0001,color="gray",zorder=10)
            delta_original.append(h*tf[i]-h*ti[i])
            nni=nni+1
    if tipo[i]==-1:
        if h*tf[i]-h*ti[i]<0.84:
            delta_perp.append(0)
        else:
            ax.hlines(y=nni*0.001, xmin=h*ti[i], xmax=h*tf[i], linewidth=.2, color='k')
            ax.scatter(h*tf[i], nni*0.001,marker='o', s=.5,color='k', zorder=20)
            ax.scatter(h*ti[i], nni*0.001,marker='o', s=.0001,color="gray",zorder=10)
            delta_perp.append(h*tf[i]-h*ti[i])
            nni=nni+1


np.savetxt('./outputsT1/lifetime_unresolved.txt',np.c_[range(len(delta_unresolved)),delta_unresolved,delta_unresolved_vi], fmt='%1.5f')

ax.set(**pparam)
plt.xlim(0,350)
ax.set_xticks([0,100,200,300])

plt.tight_layout()
fig.savefig('FigureDelayedEvents.png',dpi=500)
