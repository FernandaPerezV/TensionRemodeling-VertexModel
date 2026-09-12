import numpy as np
import matplotlib.pyplot as plt

class Vertex:
    def __init__(self, pos):
        x, y = pos
        if x < 0.0:
            self.x = L_x + x
        elif x > L_x:
            self.x = -L_x + x
        else:
            self.x = x

        if y < 0.0:
            self.y = L_y + y
        elif y > L_y:
            self.y = -L_y + y
        else:
            self.y = y

class Tissue:
    def __init__(self, vertices, topology_i, topology_j):
        self.vs = vertices
        self.ti = topology_i
        self.tj = topology_j

    def tensiones(self):
        for i in range(len(self.ti)):
            v1 = int(self.ti[i])
            v2 = int(self.tj[i])
            if v1 in unresolved:
                plt.scatter(self.vs[v1].x, self.vs[v1].y, s=20, c='r', zorder=50,alpha=1)
            if v2 in unresolved:
                plt.scatter(self.vs[v1].x, self.vs[v1].y, s=20, c='r', zorder=50,alpha=1)

            x1,y1=self.vs[v1].x,self.vs[v1].y
            x2,y2=self.vs[v2].x,self.vs[v2].y

            if np.abs(x2-x1)>L_x/2 or np.abs(y2-y1)>L_y/2:
                if np.abs(y2-y1)<L_y/2: 
                    if x2-x1>L_x/2:
                        x2left = x2-L_x
                        x1right = x1+L_x
                        plt.plot([x1right,x2],
                                 [y1,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2left],
                                 [y1,y2], '-', color='k',zorder=10)
                    else:
                        x2left = x2+L_x
                        x1right = x1-L_x
                        plt.plot([x1right,x2],
                                 [y1,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2left],
                                 [y1,y2], '-', color='k',zorder=10)

                if np.abs(x2-x1)<L_x/2:
                    if y1-y2>L_y/2:
                        y2left = y2+L_y
                        y1right = y1-L_y
                        plt.plot([x1,x2],
                                 [y1right,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2],
                                 [y1,y2left], '-', color='k',zorder=10)

                    else:
                        y2left = y2-L_y
                        y1right = y1+L_y
                        plt.plot([x1,x2],
                                 [y1right,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2],
                                 [y1,y2left], '-', color='k',zorder=10)

                else: 
                    if x2-x1>L_x/2 and y2-y1>L_y/2:
                        x2left = x2-L_x
                        x1right = x1+L_x
                        y2left = y2-L_y
                        y1right = y1+L_y
                        plt.plot([x1right,x2],
                                 [y1right,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2left],
                                 [y1,y2left], '-', color='k',zorder=10)
                        plt.plot([x1right,x2],
                                 [y1right-L_y,y2-L_y], '-', color='k',zorder=10)

                        plt.plot([x1,x2left],
                                 [y1+L_y,y2left+L_y], '-', color='k',zorder=10)

                    elif x2-x1>L_x/2 and y1-y2>L_y/2:
                        x2left = x2-L_x
                        x1right = x1+L_x
                        y2left = y2+L_y
                        y1right = y1-L_y
                        plt.plot([x1right,x2],
                                 [y1right,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2left],
                                 [y1,y2left], '-', color='k',zorder=10)

                    elif x1-x2>L_x/2 and y2-y1>L_y/2:
                        x2left = x2+L_x
                        x1right = x1-L_x
                        y2left = y2-L_y
                        y1right = y1+L_y
                        plt.plot([x1right,x2],
                                 [y1right,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2left],
                                 [y1,y2left], '-', color='k',zorder=10)

                    elif x1-x2>L_x/2 and y1-y2>L_y/2:
                        x2left = x2+L_x
                        x1right = x1-L_x
                        y2left = y2+L_y
                        y1right = y1-L_y
                        plt.plot([x1right,x2],
                                 [y1right,y2], '-', color='k',zorder=10)
                        plt.plot([x1,x2left],
                                 [y1,y2left], '-', color='k',zorder=10)

                        plt.plot([x1right,x2],
                                 [y1right-L_y,y2-L_y], '-', color='k',zorder=10)

                        plt.plot([x1,x2left],
                                 [y1+L_y,y2left+L_y], '-', color='k',zorder=10)

            else:
                plt.plot([x1,x2],
                         [y1,y2], '-', color='k',zorder=10)
                plt.plot([x1,x2],
                         [y1,y2], '-', color='k',zorder=10)


def crear_T(texto_vertices,texto_topologia):
    data = np.loadtxt(texto_vertices)
    topology = np.loadtxt(texto_topologia)
    nros_c = data[:,0]

    t = []
    for i in range(int(nros_c[len(nros_c)-1])+1):
        t.append(Vertex((data[i,1],data[i,2])))

    T1 = Tissue(t,topology[:,0],topology[:,1])
    T1.tensiones()


plt.close('all')

factor = np.sqrt(2/(3*np.sqrt(3)))

L_x = 19* np.sqrt(3)*factor
L_y = 13.*3 *factor

unresolved_all = np.loadtxt('./outputsT1/lifetime_unresolved.txt')[:,2]
unresolved_all_time = np.loadtxt('./outputsT1/lifetime_unresolved.txt')[:,1]
unresolved = []

for i in range(len(unresolved_all)):
    if unresolved_all_time[i]>100:
        unresolved.append(unresolved_all[i])


for ii in [119]: #119 used in PRX
    fig, ax = plt.subplots(figsize=(5,3.8))

    i = int(ii)
    crear_T('./outputsVertices/'+str(i)+'_vertices.txt','./outputsJunctions/'+str(i)+'.txt')
    plt.tick_params(
    axis='both',
    which='both',
    bottom=False,
    top=False,
    left=False,
    labelleft=False,
    labelbottom=False)

    plt.axis('scaled')
    plt.xlim([0,L_x])
    plt.ylim([0,L_y])
    plt.savefig('LastSnapshot.png', dpi=400)

    plt.close('all')
