from __future__ import division
import numpy as np

np.random.seed(111)

class Vertex:
    def __init__(self,n, vert_cell, velocidad, celulas,adyacentes, dt, eje, par_4fold):
        self.id = n
        self.x = vert_cell[0]
        self.y = vert_cell[1]
        self.z = vert_cell[2]
        self.r = (self.x, self.y, self.z)
        self.vx = velocidad[0]
        self.vy = velocidad[1]
        self.vz = velocidad[2]
        self.cells = celulas
        self.ady = adyacentes
        self.dt = dt
        self.eje = eje
        self.par = par_4fold

    def v_nuevas_tipo12(self):
        vx = 0
        vy = 0
        vz = 0

        for j in range(len(self.cells)):
            k_a =  celulas[self.cells[j]].KA
            c = self.cells[j]
            vi = int(self.ady[j][0])
            vd = int(self.ady[j][1])

            l_ii = modulo(resta_vectores(t[vi].r,self.r))
            l_id = modulo(resta_vectores(t[vd].r,self.r))
            if k_a == 0:
                T_ii = 0
                T_id = 0
            else:
                T_ii = (Tensiones[self.id][vi]+ TensionesFluct[self.id][vi]+ Tensiones4fold[self.id][vi] + Tact[self.id][vi]* l_ii) * 0.5
                T_id = (Tensiones[self.id][vd]+ TensionesFluct[self.id][vd]+ Tensiones4fold[self.id][vd] + Tact[self.id][vd]* l_id) * 0.5

            a = resta_vectores(t[vi].r,self.r)
            b = resta_vectores(t[vd].r,self.r)

            vx = vx +  T_ii * a[0]/ l_ii
            vy = vy +  T_ii * a[1]/ l_ii

            vx = vx +  T_id * b[0]/ l_id
            vy = vy +  T_id * b[1]/ l_id

            #Parte del area
            k = (0.0, 0.0, 1.0)
            rest = resta_vectores(t[vd].r, t[vi].r)
            v_tot = np.cross(rest,k)
            area0_c = area0_ec_mov[c]
            AREA = 0.5* k_a * (celulas[c].area - area0_c)
            vx = vx +  AREA * v_tot[0]
            vy = vy +  AREA * v_tot[1]
            vz = 0

        return [vx/mu, vy/mu, vz/mu]

class Cell(Vertex):
    def __init__(self, numero, indices,  dt, A, KA, KP):
        self.n = numero
        self.ind = indices
        self.dt = dt
        self.area = A
        self.KA = KA
        self.KP = KP

    def cal_area(self):
        v = self.ind
        area = 0
        vects = []
        for i in range(len(v)-2):
            vects.append(resta_vectores(t[int(v[i+1])].r,t[int(v[0])].r))
        for j in range(len(vects)-1):
            k = (0.0, 0.0, 1.0)
            area = area - 0.5 * np.dot(np.cross(vects[j],vects[j+1]), k)
        self.area = area

    def cal_centro(self):
        c = (0,0,0)
        v = self.ind
        for i in range(len(v)-1):
            add = resta_vectores(t[int(v[i])].r,t[int(v[0])].r)
            c = suma(c,add)
        posible_centro = suma(t[int(v[0])].r,div(c,len(v)-1))
        centroencaja = resta_inicial([posible_centro[0],posible_centro[1]])
        return [centroencaja[0],centroencaja[1],0]

class Tissue(Cell):
    def __init__(self, vertices, celulas, dt):
        self.vs = vertices
        self.cs = celulas
        self.dt = dt

    def T1_tejido(self, now):
        flats = []
        for i in range(n_vertices):
            if len(self.vs[i].ady) == 0:
                flat = []
            else:
                flat = [self.vs[i].ady[0][0],self.vs[i].ady[0][1],self.vs[i].ady[1][0]]
            flats.append(flat)

        for i in range(n_vertices):
            for j in flats[i]:

                if j>i and j in np.array(self.vs[i].ady).flatten():
                    lij = modulo(resta_vectores(self.vs[i].r, self.vs[j].r))

                    if lij < l_T1:
                        if len(self.vs[i].cells)==3 and len(self.vs[j].cells)==3:
                            if previo_T1(i,j)==1:
                                t4f_saved= Tensiones[i][j]+TensionesFluct[i][j]+Tact[i][j]*lij
                                l4f_saved = lij
                                self.vs[i].par[1] = t4f_saved
                                self.vs[i].par[2] = l4f_saved

                                delete_v1_original(i,j,0,0) 

                                file = open('./outputsT1/4-fold.txt', 'a')
                                n1,n2,n3,n4=round(now,3), int(i), int(j), round(lij,4)
                                n5,n6,n7 = round(Tensiones[i][j],4), round(Tact[i][j],4), round(TensionesFluct[i][j],4)
                                file.write(str(n1)+' '+str(n2)+' '+str(n3)+' '+str(n4)+' '+str(n5)+' '+str(n6)+' '+str(n7)+ '\n')

    def try_T1_tejido(self, now):
        n = []
        for i in range(n_vertices):
            if len(self.vs[i].cells)==4:
                n.append(i)
        for i in range(len(n)):
            vi1 = n[i]
            vi2 = self.vs[n[i]].par[0]
            eje = self.vs[vi1].eje


            nT12 = np.random.normal(0.1,0.01)
            nl120 = np.random.normal(l_T1,l_T1*0.1)
            if nl120 <0:
                nl120=0
            if nT12 <0:
                nT12=0
            producto1= create_v1_original(vi1,nT12,nl120)
            Textra = nT12+Gamma_act*1.5*l_T1
            delete_v1_original(vi1,vi2,Textra,self.vs[n[i]].par[1])
            producto2= create_v1_perpendicular(vi1,nT12,nl120)

            if producto1<0 and producto2<0:
                delete_v1_perpendicular(vi1,vi2,Textra,self.vs[n[i]].par[1])
                self.vs[vi1].eje = eje

            else:
                if producto1>producto2:
                    delete_v1_perpendicular(vi1,vi2,Textra,self.vs[n[i]].par[1])
                    self.vs[vi1].eje = eje
                    create_v1_original(vi1,nT12,nl120)
                    self.vs[vi1].par = [-1,0,0]
                    self.vs[vi1].eje = [-1,-1]
                    file = open('./outputsT1/T1_original.txt', 'a')
                    n1,n2,n3=round(now,3), int(vi1), int(vi2)
                    file.write(str(n1)+' '+str(n2)+' '+str(n3)+ '\n')

                else:
                    file = open('./outputsT1/T1_perpendicular.txt', 'a')
                    self.vs[vi1].par = [-1,0,0]
                    self.vs[vi1].eje = [-1,-1]
                    n1,n2,n3=round(now,3), int(vi1), int(vi2)
                    file.write(str(n1)+' '+str(n2)+' '+str(n3)+ '\n')

    def evol_vertex(self,now):
        self.pos_nuevas_vertex()
        self.cal_area()
        self.T1_tejido(now)
        self.nuevos_LargosNat_Tensiones()
        self.cal_area()

    def pos_nuevas_vertex(self):
        velocidades_vm1 = []
        for n in range(len(self.vs)):
            if len(self.vs[n].cells)==0:
                velocidades_vm1.append([0,0,0])
            else:
                velocidades_vm1.append(self.vs[n].v_nuevas_tipo12())

        q1 = 0
        for n in range(len(self.vs)):
            vels1 = velocidades_vm1[q1]
            q1 = q1 + 1
            self.vs[n].vx = vels1[0]
            self.vs[n].vy = vels1[1]
            self.vs[n].vz = vels1[2]

            self.vs[n].x = self.vs[n].x + self.vs[n].vx*self.vs[n].dt
            self.vs[n].y = self.vs[n].y + self.vs[n].vy*self.vs[n].dt
            self.vs[n].z = self.vs[n].z + self.vs[n].vz*self.vs[n].dt
            self.vs[n].r = (self.vs[n].x, self.vs[n].y, self.vs[n].z)

            self.vs[n].r = resta_inicial(self.vs[n].r)
            self.vs[n].x = self.vs[n].r[0]
            self.vs[n].y = self.vs[n].r[1]
            self.vs[n].z = self.vs[n].r[2]

    def data(self):
        pos = []
        new01 = []
        for i in range(len(self.vs)):
            pos.append((self.vs[i].x,self.vs[i].y))
            new01.append(i)

        np.savetxt('./outputsVertices/data_vertices.txt',  np.c_[new01, pos], fmt='%1.10f')

    def cal_area(self):
        for i in cs_moviles:
            self.cs[i].cal_area()

    def nuevos_LargosNat_Tensiones(self):
        for i in range(len(self.vs)):
            if len(self.vs[i].cells)==3:
                for j in [self.vs[i].ady[0][0],self.vs[i].ady[0][1],self.vs[i].ady[1][0]]:
                    if j>i and len(self.vs[j].cells)==3:
                        lij = modulo(resta_vectores(self.vs[i].r, self.vs[j].r))
                        e_ij = (lij - LargosNat[i][j])/LargosNat[i][j]
                        if np.abs(e_ij) < e_critico:
                            alpha = 0
                        else:
                            if e_ij>=0:
                                alpha = Ke
                            else:
                                if lij<l_T1:
                                    alpha=0
                                else:
                                    alpha = Kc

                        v_Lij = -Kl*(LargosNat[i][j]-lij)
                        v_Tij = -alpha*(lij-LargosNat[i][j]) - (Tensiones[i][j]-T0)/tau
                        new_Lij = LargosNat[i][j] + self.dt * v_Lij
                        new_Tij = Tensiones[i][j] + self.dt * v_Tij
                        LargosNat[i][j] = new_Lij
                        LargosNat[j][i] = new_Lij

                        new2_Tij = TensionesFluct[i][j]*(1-(dt/tau_fluct)) + np.sqrt(2*difusion*dt)*np.random.normal(0,1)
                        TensionesFluct[i][j] = new2_Tij
                        TensionesFluct[j][i] = new2_Tij

                        if new_Tij+new2_Tij+Tact[i][j]*lij<0:
                            new_Tij = -1*(Tact[i][j]*lij + new2_Tij)

                        Tensiones[i][j] = new_Tij
                        Tensiones[j][i] = new_Tij


            if len(self.vs[i].cells)==4:
                for j in [self.vs[i].ady[0][0],self.vs[i].ady[0][1],self.vs[i].ady[1][0],self.vs[i].ady[2][0]]:
                    if len(self.vs[j].cells)==3:
                        lij = modulo(resta_vectores(self.vs[i].r, self.vs[j].r))
                        e_ij = (lij - LargosNat[i][j])/LargosNat[i][j]
                        if np.abs(e_ij) < e_critico:
                            alpha = 0
                        else:
                            if e_ij>=0:
                                alpha = Ke
                            else:
                                if lij<l_T1:
                                    alpha=0
                                else:
                                    alpha = Kc

                        v_Lij = -Kl*(LargosNat[i][j]-lij)
                        v_Tij = -alpha*(lij-LargosNat[i][j])- (Tensiones[i][j]-T0)/tau

                        new_Lij = LargosNat[i][j] + self.dt * v_Lij
                        new_Tij = Tensiones[i][j] + self.dt * v_Tij
                        LargosNat[i][j] = new_Lij
                        LargosNat[j][i] = new_Lij

                        new2_Tij = TensionesFluct[i][j]*(1-(dt/tau_fluct)) + np.sqrt(2*difusion*dt)*np.random.normal(0,1)
                        TensionesFluct[i][j] = new2_Tij
                        TensionesFluct[j][i] = new2_Tij

                        if new_Tij+new2_Tij+Tact[i][j]*lij<0:
                            new_Tij = -1*(Tact[i][j]*lij + new2_Tij)

                        Tensiones[i][j] = new_Tij
                        Tensiones[j][i] = new_Tij

                    if len(self.vs[j].cells)==4 and j>i:
                        lij = modulo(resta_vectores(self.vs[i].r, self.vs[j].r))
                        e_ij = (lij - LargosNat[i][j])/LargosNat[i][j]
                        if np.abs(e_ij) < e_critico:
                            alpha = 0
                        else:
                            if e_ij>=0:
                                alpha = Ke
                            else:
                                if lij<l_T1:
                                    alpha=0
                                else:
                                    alpha = Kc

                        v_Lij = -Kl*(LargosNat[i][j]-lij)
                        v_Tij = -alpha*(lij-LargosNat[i][j])- (Tensiones[i][j]-T0)/tau

                        new_Lij = LargosNat[i][j] + self.dt * v_Lij
                        new_Tij = Tensiones[i][j] + self.dt * v_Tij
                        LargosNat[i][j] = new_Lij
                        LargosNat[j][i] = new_Lij

                        new2_Tij = TensionesFluct[i][j]*(1-(dt/tau_fluct)) + np.sqrt(2*difusion*dt)*np.random.normal(0,1)
                        TensionesFluct[i][j] = new2_Tij
                        TensionesFluct[j][i] = new2_Tij

                        if new_Tij+new2_Tij+Tact[i][j]*lij <0:
                            new_Tij = -1*(Tact[i][j]*lij + new2_Tij )

                        Tensiones[i][j] = new_Tij
                        Tensiones[j][i] = new_Tij

def create_v1_perpendicular(vi1,nT12,nl120):
    cells_i1 = T1.vs[vi1].cells 
    cellA = cells_i1[0]
    cellB = cells_i1[1]
    cellC = cells_i1[2]
    cellD = cells_i1[3]
    vert_A = T1.cs[cellA].ind
    vert_B = T1.cs[cellB].ind
    vert_C = T1.cs[cellC].ind
    vert_D = T1.cs[cellD].ind


    for i in vert_A:
        if i in vert_B and i not in vert_C and i not in vert_D:
            vert2 = i
    for i in vert_A:
        if i in vert_D and i not in vert_B and i not in vert_C:
            vert3 = i
    for i in vert_C:
        if i in vert_B and i not in vert_A and i not in vert_D:
            vert6 = i
    for i in vert_C:
        if i in vert_D and i not in vert_A and i not in vert_B:
            vert5 = i
    vert4 = vi1
    vi2 = T1.vs[vi1].par[0]
    vert1 = T1.vs[vi1].par[0]
    t4f = T1.vs[vi1].par[1]


    c1y2 = []
    for c in [cellA,cellB,cellC,cellD]:
        if c not in T1.vs[vert4].eje:
            c1y2.append(c)
    c1 = c1y2[0]
    c2 = c1y2[1]
    r_T1 = resta_vectores(celulas[c1].cal_centro(),celulas[c2].cal_centro())
    rnorm_T1 = div(r_T1,modulo(r_T1))

    #MODIFY ADY
    for i in range(len(T1.vs[vert3].ady)):
        for j in range(len(T1.vs[vert3].ady[i])):
            if T1.vs[vert3].ady[i][j] == vert4:
                T1.vs[vert3].ady[i][j] = vert1
    for i in range(len(T1.vs[vert5].ady)):
        for j in range(len(T1.vs[vert5].ady[i])):
            if T1.vs[vert5].ady[i][j] == vert4:
                T1.vs[vert5].ady[i][j] = vert1

    T1.vs[vert4].cells = [cellA, cellB, cellC]
    T1.vs[vert4].ady = [[vert2,vert1],[vert6,vert2],[vert1,vert6]]
    T1.vs[vert1].cells = [cellC, cellD, cellA]
    T1.vs[vert1].ady = [[vert5,vert4],[vert3,vert5],[vert4,vert3]]

    #MODIFY TENSIONES AND NATURAL LENGTHS
    Tensiones[vert3][vert1] = Tensiones[vert3][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert1][vert3] = Tensiones[vert3][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert5][vert1] = Tensiones[vert5][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert1][vert5] = Tensiones[vert5][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert3][vert4]= 0
    Tensiones[vert4][vert3]= 0
    Tensiones[vert5][vert4]= 0
    Tensiones[vert4][vert5]= 0
    Tensiones[vert2][vert4]= Tensiones[vert2][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert4][vert2]= Tensiones[vert4][vert2]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert6][vert4]= Tensiones[vert6][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert4][vert6]= Tensiones[vert4][vert6]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    TensionesFluct[vert3][vert1] = TensionesFluct[vert3][vert4]
    TensionesFluct[vert1][vert3] = TensionesFluct[vert3][vert4]
    TensionesFluct[vert5][vert1] = TensionesFluct[vert5][vert4]
    TensionesFluct[vert1][vert5] = TensionesFluct[vert5][vert4]
    TensionesFluct[vert3][vert4]= 0
    TensionesFluct[vert4][vert3]= 0
    TensionesFluct[vert5][vert4]= 0
    TensionesFluct[vert4][vert5]= 0

    Tensiones4fold[vert6][vert4]= 0
    Tensiones4fold[vert4][vert6]= 0
    Tensiones4fold[vert5][vert4]= 0
    Tensiones4fold[vert4][vert5]= 0
    Tensiones4fold[vert2][vert4]= 0
    Tensiones4fold[vert4][vert2]= 0
    Tensiones4fold[vert3][vert4]= 0
    Tensiones4fold[vert4][vert3]= 0



    LargosNat[vert3][vert1]= LargosNat[vert3][vert4]
    LargosNat[vert1][vert3]= LargosNat[vert3][vert4]
    LargosNat[vert5][vert1]= LargosNat[vert5][vert4]
    LargosNat[vert1][vert5]= LargosNat[vert5][vert4]
    LargosNat[vert3][vert4]= 0
    LargosNat[vert4][vert3]= 0
    LargosNat[vert5][vert4]= 0
    LargosNat[vert4][vert5]= 0

    #MODIFY IND cell A
    add_A = []
    q = 0
    for i in range(len(T1.cs[cellA].ind)-1):
        if T1.cs[cellA].ind[i] == vert4 and q==0:
            add_A.append(vert4)
            add_A.append(vert1)
            q=1
        else:
            add_A.append(T1.cs[cellA].ind[i])
    T1.cs[cellA].ind = add_A
    if T1.cs[cellA].ind[0] != T1.cs[cellA].ind[len(T1.cs[cellA].ind)-1]:
        T1.cs[cellA].ind.append(T1.cs[cellA].ind[0])

    #MODIFY IND cell C
    add_C = []
    q = 0
    for i in range(len(T1.cs[cellC].ind)-1):
        if T1.cs[cellC].ind[i] == vert5 and q==0:
            add_C.append(vert5)
            add_C.append(vert1)
            q=1
        else:
            add_C.append(T1.cs[cellC].ind[i])
    T1.cs[cellC].ind = add_C
    if T1.cs[cellC].ind[0] != T1.cs[cellC].ind[len(T1.cs[cellC].ind)-1]:
        T1.cs[cellC].ind.append(T1.cs[cellC].ind[0])

    #MODIFY IND cell D
    add_D = []
    q = 0
    for i in range(len(T1.cs[cellD].ind)-1):
        if T1.cs[cellD].ind[i] == vert4 and q==0:
            add_D.append(vert1)
            q=1
        else:
            add_D.append(T1.cs[cellD].ind[i])
    T1.cs[cellD].ind = add_D
    if T1.cs[cellD].ind[0] != T1.cs[cellD].ind[len(T1.cs[cellD].ind)-1]:
        T1.cs[cellD].ind.append(T1.cs[cellD].ind[0])

    T1.vs[vert4].eje = [-1,-1]
    T1.vs[vert4].par[0] = -1


    pos_vi1p = suma(T1.vs[vert4].r,mult(rnorm_T1,(3/4)*l_T1))
    pos_vi2p = suma(T1.vs[vert4].r,mult(rnorm_T1,(-3/4)*l_T1))

    t[vi1].r = pos_vi1p
    t[vi1].x = pos_vi1p[0]
    t[vi1].y = pos_vi1p[1]
    t[vi1].z = pos_vi1p[2]

    t[vi2].r = pos_vi2p
    t[vi2].x = pos_vi2p[0]
    t[vi2].y = pos_vi2p[1]
    t[vi2].z = pos_vi2p[2]

    for c in [cellA,cellB,cellC,cellD]:
        T1.cs[c].cal_area()

    Tensiones[vert4][vert1] = nT12
    Tensiones[vert1][vert4] = nT12
    TensionesFluct[vert4][vert1] = 0
    TensionesFluct[vert1][vert4] = 0
    LargosNat[vert4][vert1] = nl120
    LargosNat[vert1][vert4] = nl120

    vel_i = T1.vs[vert4].v_nuevas_tipo12()
    vel_j = T1.vs[vert1].v_nuevas_tipo12()
    producto = np.dot(resta_vectores(T1.vs[vert4].r,T1.vs[vert1].r), resta(vel_i,vel_j))

    return producto

def create_v1_original(vi1,nT12,nl120):
    cells_i1 = T1.vs[vi1].cells 
    cellA = cells_i1[0]
    cellB = cells_i1[1]
    cellC = cells_i1[2]
    cellD = cells_i1[3]

    vert_A = T1.cs[cellA].ind
    vert_B = T1.cs[cellB].ind
    vert_C = T1.cs[cellC].ind
    vert_D = T1.cs[cellD].ind


    for i in vert_A:
        if i in vert_B and i not in vert_C and i not in vert_D:
            vert2 = i
    for i in vert_A:
        if i in vert_D and i not in vert_B and i not in vert_C:
            vert3 = i
    for i in vert_C:
        if i in vert_B and i not in vert_A and i not in vert_D:
            vert6 = i
    for i in vert_C:
        if i in vert_D and i not in vert_A and i not in vert_B:
            vert5 = i
    vert4 = vi1
    vi2 = T1.vs[vi1].par[0]
    vert1 = T1.vs[vi1].par[0]
    t4f = T1.vs[vi1].par[1]


    c1 = T1.vs[vert4].eje[0]
    c2 = T1.vs[vert4].eje[1]
    r_T1 = resta_vectores(celulas[c1].cal_centro(),celulas[c2].cal_centro())
    rnorm_T1 = div(r_T1,modulo(r_T1))

    #MODIFY ADY
    for i in range(len(T1.vs[vert6].ady)):
        for j in range(len(T1.vs[vert6].ady[i])):
            if T1.vs[vert6].ady[i][j] == vert4:
                T1.vs[vert6].ady[i][j] = vert1
    for i in range(len(T1.vs[vert5].ady)):
        for j in range(len(T1.vs[vert5].ady[i])):
            if T1.vs[vert5].ady[i][j] == vert4:
                T1.vs[vert5].ady[i][j] = vert1

    T1.vs[vert4].cells = [cellA, cellB, cellD]
    T1.vs[vert4].ady = [[vert2,vert3],[vert1,vert2],[vert3,vert1]]
    T1.vs[vert1].cells = [cellC, cellD, cellB]
    T1.vs[vert1].ady = [[vert5,vert6],[vert4,vert5],[vert6,vert4]]

    #MODIFY TENSIONES AND NATURAL LENGTHS
    Tensiones[vert6][vert1] = Tensiones[vert6][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert1][vert6] = Tensiones[vert6][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert5][vert1] = Tensiones[vert5][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert1][vert5] = Tensiones[vert5][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert6][vert4]= 0
    Tensiones[vert4][vert6]= 0
    Tensiones[vert5][vert4]= 0
    Tensiones[vert4][vert5]= 0
    Tensiones[vert2][vert4]= Tensiones[vert2][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert4][vert2]= Tensiones[vert4][vert2]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert3][vert4]= Tensiones[vert3][vert4]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    Tensiones[vert4][vert3]= Tensiones[vert4][vert3]-0.25*(nT12+Gamma_act*1.5*l_T1)+0.25*t4f
    TensionesFluct[vert6][vert1] = TensionesFluct[vert6][vert4]
    TensionesFluct[vert1][vert6] = TensionesFluct[vert6][vert4]
    TensionesFluct[vert5][vert1] = TensionesFluct[vert5][vert4]
    TensionesFluct[vert1][vert5] = TensionesFluct[vert5][vert4]
    TensionesFluct[vert6][vert4]= 0
    TensionesFluct[vert4][vert6]= 0
    TensionesFluct[vert5][vert4]= 0
    TensionesFluct[vert4][vert5]= 0

    Tensiones4fold[vert6][vert4]= 0
    Tensiones4fold[vert4][vert6]= 0
    Tensiones4fold[vert5][vert4]= 0
    Tensiones4fold[vert4][vert5]= 0
    Tensiones4fold[vert2][vert4]= 0
    Tensiones4fold[vert4][vert2]= 0
    Tensiones4fold[vert3][vert4]= 0
    Tensiones4fold[vert4][vert3]= 0

    LargosNat[vert6][vert1]= LargosNat[vert6][vert4]
    LargosNat[vert1][vert6]= LargosNat[vert6][vert4]
    LargosNat[vert5][vert1]= LargosNat[vert5][vert4]
    LargosNat[vert1][vert5]= LargosNat[vert5][vert4]
    LargosNat[vert6][vert4]= 0
    LargosNat[vert4][vert6]= 0
    LargosNat[vert5][vert4]= 0
    LargosNat[vert4][vert5]= 0

    #MODIFY IND cell B
    add_B = []
    q = 0
    for i in range(len(T1.cs[cellB].ind)-1):
        if T1.cs[cellB].ind[i] == vert6 and q==0:
            add_B.append(vert6)
            add_B.append(vert1)
            q=1
        else:
            add_B.append(T1.cs[cellB].ind[i])
    T1.cs[cellB].ind = add_B
    if T1.cs[cellB].ind[0] != T1.cs[cellB].ind[len(T1.cs[cellB].ind)-1]:
        T1.cs[cellB].ind.append(T1.cs[cellB].ind[0])

    #MODIFY IND cell C
    add_C = []
    q = 0
    for i in range(len(T1.cs[cellC].ind)-1):
        if T1.cs[cellC].ind[i] == vi1 and q==0:
            add_C.append(vert1)
            q=1
        else:
            add_C.append(T1.cs[cellC].ind[i])
    T1.cs[cellC].ind = add_C
    if T1.cs[cellC].ind[0] != T1.cs[cellC].ind[len(T1.cs[cellC].ind)-1]:
        T1.cs[cellC].ind.append(T1.cs[cellC].ind[0])

    #MODIFY IND cell D
    add_D = []
    q = 0
    for i in range(len(T1.cs[cellD].ind)-1):
        if T1.cs[cellD].ind[i] == vi1 and q==0:
            add_D.append(vi1)
            add_D.append(vert1)
            q=1
        else:
            add_D.append(T1.cs[cellD].ind[i])
    T1.cs[cellD].ind = add_D
    if T1.cs[cellD].ind[0] != T1.cs[cellD].ind[len(T1.cs[cellD].ind)-1]:
        T1.cs[cellD].ind.append(T1.cs[cellD].ind[0])

    T1.vs[vert4].eje = [-1,-1]
    T1.vs[vert4].par[0] = -1

    pos_vi1p = suma(T1.vs[vert4].r,mult(rnorm_T1,(3/4)*l_T1))
    pos_vi2p = suma(T1.vs[vert4].r,mult(rnorm_T1,(-3/4)*l_T1))


    t[vi1].r = pos_vi1p
    t[vi1].x = pos_vi1p[0]
    t[vi1].y = pos_vi1p[1]
    t[vi1].z = pos_vi1p[2]

    t[vi2].r = pos_vi2p
    t[vi2].x = pos_vi2p[0]
    t[vi2].y = pos_vi2p[1]
    t[vi2].z = pos_vi2p[2]

    for c in [cellA,cellB,cellC,cellD]:
        T1.cs[c].cal_area()

    Tensiones[vert4][vert1] = nT12
    Tensiones[vert1][vert4] = nT12
    TensionesFluct[vert4][vert1] = 0
    TensionesFluct[vert1][vert4] = 0
    LargosNat[vert4][vert1] = nl120
    LargosNat[vert1][vert4] = nl120

    vel_i = T1.vs[vert4].v_nuevas_tipo12()
    vel_j = T1.vs[vert1].v_nuevas_tipo12()
    producto = np.dot(resta_vectores(T1.vs[vert4].r,T1.vs[vert1].r), resta(vel_i,vel_j))

    return producto

def delete_v1_original(vi1,vi2,Textra,Tmenos): 
    t4f = T1.vs[vi1].par[1]
    not_repeated = []
    cells_i1 = T1.vs[vi1].cells 
    cells_i2 = T1.vs[vi2].cells 
    for i in range(len(cells_i1)):
        if cells_i1[i] not in cells_i2:
            not_repeated.append(cells_i1[i]) 
    for i in range(len(cells_i2)):
        if cells_i2[i] not in cells_i1:
            not_repeated.append(cells_i2[i]) 
    if cells_i1[0]== not_repeated[0]: 
        cells_i1_new = [cells_i1[0],cells_i1[1],not_repeated[1]]
    elif cells_i1[1]== not_repeated[0]: 
        cells_i1_new = [cells_i1[1],cells_i1[2],not_repeated[1]]
    else:  
        cells_i1_new = [cells_i1[2],cells_i1[0],not_repeated[1]]
    if cells_i2[0]== not_repeated[1]: 
        cells_i2_new = [cells_i2[0],cells_i2[1],not_repeated[0]]
    elif cells_i2[1]== not_repeated[1]: 
        cells_i2_new = [cells_i2[1],cells_i2[2],not_repeated[0]]
    else: 
        cells_i2_new = [cells_i2[2],cells_i2[0],not_repeated[0]]

    cellA = cells_i1_new[0]
    cellB = cells_i1_new[1]
    cellC = cells_i1_new[2]
    cellD = cells_i2_new[1]
    vert_A = T1.cs[cellA].ind[:len(T1.cs[cellA].ind)-1]
    vert_B = T1.cs[cellB].ind[:len(T1.cs[cellB].ind)-1]
    vert_C = T1.cs[cellC].ind[:len(T1.cs[cellC].ind)-1]
    vert_D = T1.cs[cellD].ind[:len(T1.cs[cellD].ind)-1]
  

    for i in vert_A:
        if i in vert_B and i not in vert_C and i not in vert_D:
            vert2 = i
    for i in vert_A:
        if i in vert_D and i not in vert_B and i not in vert_C:
            vert3 = i
    for i in vert_C:
        if i in vert_B and i not in vert_A and i not in vert_D:
            vert6 = i
    for i in vert_C:
        if i in vert_D and i not in vert_A and i not in vert_B:
            vert5 = i
    vert4 = vi1
    vert1 = vi2

    #MODIFY ADY
    for i in range(len(T1.vs[vert6].ady)):
        for j in range(len(T1.vs[vert6].ady[i])):
            if T1.vs[vert6].ady[i][j] == vert1:
                T1.vs[vert6].ady[i][j] = vert4

    for i in range(len(T1.vs[vert5].ady)):
        for j in range(len(T1.vs[vert5].ady[i])):
            if T1.vs[vert5].ady[i][j] == vert1:
                T1.vs[vert5].ady[i][j] = vert4
    T1.vs[vert4].cells = [cellA, cellB, cellC, cellD]
    T1.vs[vert4].ady = [[vert2,vert3],[vert6,vert2],[vert5,vert6],[vert3,vert5]]
    T1.vs[vert1].cells = []
    T1.vs[vert1].ady = []

    #MODIFY TENSIONES AND NATURAL LENGTHS
    Tensiones[vert6][vert4]= Tensiones[vert6][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert6]= Tensiones[vert6][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert5][vert4]= Tensiones[vert5][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert5]= Tensiones[vert5][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert6][vert1] = 0
    Tensiones[vert1][vert6] = 0
    Tensiones[vert5][vert1] = 0
    Tensiones[vert1][vert5] = 0
    Tensiones[vert2][vert4]= Tensiones[vert2][vert4]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert2]= Tensiones[vert4][vert2]+0.25*Textra-0.25*Tmenos
    Tensiones[vert3][vert4]= Tensiones[vert3][vert4]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert3]= Tensiones[vert4][vert3]+0.25*Textra-0.25*Tmenos
    TensionesFluct[vert6][vert4] = TensionesFluct[vert6][vert1]
    TensionesFluct[vert4][vert6] = TensionesFluct[vert6][vert1]
    TensionesFluct[vert5][vert4] = TensionesFluct[vert5][vert1]
    TensionesFluct[vert4][vert5] = TensionesFluct[vert5][vert1]
    TensionesFluct[vert6][vert1] = 0
    TensionesFluct[vert1][vert6] = 0
    TensionesFluct[vert5][vert1] = 0
    TensionesFluct[vert1][vert5] = 0

    Tensiones4fold[vert6][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert6]= 0.25*t4f
    Tensiones4fold[vert5][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert5]= 0.25*t4f
    Tensiones4fold[vert2][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert2]= 0.25*t4f
    Tensiones4fold[vert3][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert3]= 0.25*t4f


    LargosNat[vert6][vert4]= LargosNat[vert6][vert1]
    LargosNat[vert4][vert6]= LargosNat[vert6][vert1]
    LargosNat[vert5][vert4]= LargosNat[vert5][vert1]
    LargosNat[vert4][vert5]= LargosNat[vert5][vert1]
    LargosNat[vert6][vert1] = 0
    LargosNat[vert1][vert6] = 0
    LargosNat[vert5][vert1] = 0
    LargosNat[vert1][vert5] = 0

    #MODIFY IND cell B
    for i in T1.cs[cellB].ind:
        if i == vi2:
            T1.cs[cellB].ind.remove(vi2)
    if T1.cs[cellB].ind[0] != T1.cs[cellB].ind[len(T1.cs[cellB].ind)-1]:
        T1.cs[cellB].ind.append(T1.cs[cellB].ind[0])

    #MODIFY IND cell C
    add_C = []
    q = 0
    for i in range(len(T1.cs[cellC].ind)-1):
        if T1.cs[cellC].ind[i] == vi2 and q==0:
            add_C.append(vi1)
            q=1
        else:
            add_C.append(T1.cs[cellC].ind[i])
    T1.cs[cellC].ind = add_C
    if T1.cs[cellC].ind[0] != T1.cs[cellC].ind[len(T1.cs[cellC].ind)-1]:
        T1.cs[cellC].ind.append(T1.cs[cellC].ind[0])

    #MODIFY IND cell D
    for i in T1.cs[cellD].ind:
        if i == vi2:
            T1.cs[cellD].ind.remove(vi2)
    if T1.cs[cellD].ind[0] != T1.cs[cellD].ind[len(T1.cs[cellD].ind)-1]:
        T1.cs[cellD].ind.append(T1.cs[cellD].ind[0])

    T1.vs[vert4].eje = [cellA,cellC]
    T1.vs[vert4].par[0] = vert1

    mid = suma(T1.vs[vert4].r,div(resta_vectores(T1.vs[vert1].r,T1.vs[vert4].r),2))

    T1.vs[vi1].r = mid
    T1.vs[vi1].x = mid[0]
    T1.vs[vi1].y = mid[1]
    T1.vs[vi1].z = mid[2]

    T1.vs[vi2].r = mid
    T1.vs[vi2].x = mid[0]
    T1.vs[vi2].y = mid[1]
    T1.vs[vi2].z = mid[2]

def delete_v1_perpendicular(vi1,vi2,Textra,Tmenos): 
    t4f = T1.vs[vi1].par[1]
    not_repeated = []
    cells_i1 = T1.vs[vi1].cells 
    cells_i2 = T1.vs[vi2].cells 
    for i in range(len(cells_i1)):
        if cells_i1[i] not in cells_i2:
            not_repeated.append(cells_i1[i]) 
    for i in range(len(cells_i2)):
        if cells_i2[i] not in cells_i1:
            not_repeated.append(cells_i2[i]) 
    if cells_i1[0]== not_repeated[0]: 
        cells_i1_new = [cells_i1[0],cells_i1[1],not_repeated[1]]
    elif cells_i1[1]== not_repeated[0]: 
        cells_i1_new = [cells_i1[1],cells_i1[2],not_repeated[1]]
    else:  
        cells_i1_new = [cells_i1[2],cells_i1[0],not_repeated[1]]
    if cells_i2[0]== not_repeated[1]: 
        cells_i2_new = [cells_i2[0],cells_i2[1],not_repeated[0]]
    elif cells_i2[1]== not_repeated[1]: 
        cells_i2_new = [cells_i2[1],cells_i2[2],not_repeated[0]]
    else: 
        cells_i2_new = [cells_i2[2],cells_i2[0],not_repeated[0]]

    cellB = cells_i1_new[0]
    cellC = cells_i1_new[1]
    cellD = cells_i1_new[2]
    cellA = cells_i2_new[1]
    vert_A = T1.cs[cellA].ind[:len(T1.cs[cellA].ind)-1]
    vert_B = T1.cs[cellB].ind[:len(T1.cs[cellB].ind)-1]
    vert_C = T1.cs[cellC].ind[:len(T1.cs[cellC].ind)-1]
    vert_D = T1.cs[cellD].ind[:len(T1.cs[cellD].ind)-1]

   
    for i in vert_A:
        if i in vert_B and i not in vert_C and i not in vert_D:
            vert2 = i
    for i in vert_A:
        if i in vert_D and i not in vert_B and i not in vert_C:
            vert3 = i
    for i in vert_C:
        if i in vert_B and i not in vert_A and i not in vert_D:
            vert6 = i
    for i in vert_C:
        if i in vert_D and i not in vert_A and i not in vert_B:
            vert5 = i
    vert4 = vi1
    vert1 = vi2

    #MODIFY ADY
    for i in range(len(T1.vs[vert3].ady)):
        for j in range(len(T1.vs[vert3].ady[i])):
            if T1.vs[vert3].ady[i][j] == vert1:
                T1.vs[vert3].ady[i][j] = vert4

    for i in range(len(T1.vs[vert5].ady)):
        for j in range(len(T1.vs[vert5].ady[i])):
            if T1.vs[vert5].ady[i][j] == vert1:
                T1.vs[vert5].ady[i][j] = vert4

    T1.vs[vert4].cells = [cellA, cellB, cellC, cellD]
    T1.vs[vert4].ady = [[vert2,vert3],[vert6,vert2],[vert5,vert6],[vert3,vert5]]
    T1.vs[vert1].cells = []
    T1.vs[vert1].ady = []


    #MODIFY TENSIONES AND NATURAL LENGTHS
    Tensiones[vert3][vert4]= Tensiones[vert3][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert3]= Tensiones[vert3][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert5][vert4]= Tensiones[vert5][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert5]= Tensiones[vert5][vert1]+0.25*Textra-0.25*Tmenos
    Tensiones[vert3][vert1] = 0
    Tensiones[vert1][vert3] = 0
    Tensiones[vert5][vert1] = 0
    Tensiones[vert1][vert5] = 0
    Tensiones[vert2][vert4]= Tensiones[vert2][vert4]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert2]= Tensiones[vert4][vert2]+0.25*Textra-0.25*Tmenos
    Tensiones[vert6][vert4]= Tensiones[vert6][vert4]+0.25*Textra-0.25*Tmenos
    Tensiones[vert4][vert6]= Tensiones[vert4][vert6]+0.25*Textra-0.25*Tmenos
    TensionesFluct[vert3][vert4]= TensionesFluct[vert3][vert1]
    TensionesFluct[vert4][vert3]= TensionesFluct[vert3][vert1]
    TensionesFluct[vert5][vert4]= TensionesFluct[vert5][vert1]
    TensionesFluct[vert4][vert5]= TensionesFluct[vert5][vert1]
    TensionesFluct[vert3][vert1] = 0
    TensionesFluct[vert1][vert3] = 0
    TensionesFluct[vert5][vert1] = 0
    TensionesFluct[vert1][vert5] = 0

    Tensiones4fold[vert6][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert6]= 0.25*t4f
    Tensiones4fold[vert5][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert5]= 0.25*t4f
    Tensiones4fold[vert2][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert2]= 0.25*t4f
    Tensiones4fold[vert3][vert4]= 0.25*t4f
    Tensiones4fold[vert4][vert3]= 0.25*t4f

    LargosNat[vert3][vert4]= LargosNat[vert3][vert1]
    LargosNat[vert4][vert3]= LargosNat[vert3][vert1]
    LargosNat[vert5][vert4]= LargosNat[vert5][vert1]
    LargosNat[vert4][vert5]= LargosNat[vert5][vert1]
    LargosNat[vert3][vert1] = 0
    LargosNat[vert1][vert3] = 0
    LargosNat[vert5][vert1] = 0
    LargosNat[vert1][vert5] = 0

    #MODIFY IND cell C
    for i in T1.cs[cellC].ind:
        if i == vi2:
            T1.cs[cellC].ind.remove(vi2)
    if T1.cs[cellC].ind[0] != T1.cs[cellC].ind[len(T1.cs[cellC].ind)-1]:
        T1.cs[cellC].ind.append(T1.cs[cellC].ind[0])

    #MODIFY IND cell D
    add_D = []
    q = 0
    for i in range(len(T1.cs[cellD].ind)-1):
        if T1.cs[cellD].ind[i] == vi2 and q==0:
            add_D.append(vi1)
            q=1
        else:
            add_D.append(T1.cs[cellD].ind[i])
    T1.cs[cellD].ind = add_D
    if T1.cs[cellD].ind[0] != T1.cs[cellD].ind[len(T1.cs[cellD].ind)-1]:
        T1.cs[cellD].ind.append(T1.cs[cellD].ind[0])

    #MODIFY IND cell A
    for i in T1.cs[cellA].ind:
        if i == vi2:
            T1.cs[cellA].ind.remove(vi2)
    if T1.cs[cellA].ind[0] != T1.cs[cellA].ind[len(T1.cs[cellA].ind)-1]:
        T1.cs[cellA].ind.append(T1.cs[cellA].ind[0])

    T1.vs[vert4].eje = [cellA,cellC]
    T1.vs[vert4].par[0] = vert1

    mid = suma(T1.vs[vert4].r,div(resta_vectores(T1.vs[vert1].r,T1.vs[vert4].r),2))
    mid = resta_inicial(mid)

    T1.vs[vi1].r = mid
    T1.vs[vi1].x = mid[0]
    T1.vs[vi1].y = mid[1]
    T1.vs[vi1].z = mid[2]

    T1.vs[vi2].r = mid
    T1.vs[vi2].x = mid[0]
    T1.vs[vi2].y = mid[1]
    T1.vs[vi2].z = mid[2]

def previo_T1(vi1,vi2):
    not_repeated = []
    cells_i1 = T1.vs[vi1].cells 
    cells_i2 = T1.vs[vi2].cells 
    for i in range(len(cells_i1)):
        if cells_i1[i] not in cells_i2:
            not_repeated.append(cells_i1[i]) 
    for i in range(len(cells_i2)):
        if cells_i2[i] not in cells_i1:
            not_repeated.append(cells_i2[i]) 
    if cells_i1[0]== not_repeated[0]: 
        cells_i1_new = [cells_i1[0],cells_i1[1],not_repeated[1]]
    elif cells_i1[1]== not_repeated[0]: 
        cells_i1_new = [cells_i1[1],cells_i1[2],not_repeated[1]]
    else:  
        cells_i1_new = [cells_i1[2],cells_i1[0],not_repeated[1]]
    if cells_i2[0]== not_repeated[1]: 
        cells_i2_new = [cells_i2[0],cells_i2[1],not_repeated[0]]
    elif cells_i2[1]== not_repeated[1]: 
        cells_i2_new = [cells_i2[1],cells_i2[2],not_repeated[0]]
    else: 
        cells_i2_new = [cells_i2[2],cells_i2[0],not_repeated[0]]

    cellA = cells_i1_new[0]
    cellB = cells_i1_new[1]
    cellC = cells_i1_new[2]
    cellD = cells_i2_new[1]
    vert_A = T1.cs[cellA].ind[:len(T1.cs[cellA].ind)-1]
    vert_B = T1.cs[cellB].ind[:len(T1.cs[cellB].ind)-1]
    vert_C = T1.cs[cellC].ind[:len(T1.cs[cellC].ind)-1]
    vert_D = T1.cs[cellD].ind[:len(T1.cs[cellD].ind)-1]

    if len(vert_A)>=3 and len(vert_C)>=3 and len(vert_B)>=4 and len(vert_D)>=4:
        return 1
    else:
        return -1

def r_caja_rect(x,y,u0):
    Minvers = [[1/(1-u0**2), -u0/(1-u0**2)],[-u0/(1-u0**2), 1/(1-u0**2)]]
    return (Minvers[0][0]*x+Minvers[0][1]*y, Minvers[1][0]*x+Minvers[1][1]*y, 0)

def resta_inicial(a):
    pos = r_caja_rect(a[0],a[1],u0)
    x,y = pos[0], pos[1]

    if pos[0] < 0.0:
        x = L_x + pos[0]
    if pos[0] > L_x:
        x = -L_x + pos[0]

    if pos[1] < 0.0:
        y = L_y + pos[1]
    if pos[1] > L_y:
        y = -L_y + pos[1]

    return [x,y,0]

def suma(a,b):
    c = []
    for i in range(len(a)):
        c.append(a[i] + b[i])
    return c

def resta(a,b):
    c = []
    for i in range(len(a)):
        c.append(a[i] - b[i])
    return c

def div(v,a):
    return (v[0]/a, v[1]/a, v[2]/a)

def mult(v,a):
    return (v[0]*a, v[1]*a, v[2]*a)

def resta_vectores(vec1,vec2):  
    a = r_caja_rect(vec1[0],vec1[1],u0)
    b = r_caja_rect(vec2[0],vec2[1],u0)
    #print a, b
    if abs(a[0]-b[0]) < L_x/2.:
        q = a[0]-b[0]
    if a[0]-b[0] >= L_x/2.:
        q = a[0] - L_x -b[0]
    if a[0]-b[0] <= -L_x/2.:
        q = a[0] + L_x -b[0]
    if abs(a[1]-b[1]) < L_y/2.:
        n = a[1]-b[1]
    if a[1]-b[1] >= L_y/2.:
        n = a[1] - L_y -b[1]
    if a[1]-b[1] <= -L_y/2.:
        n = a[1] + L_y -b[1]
    if abs(a[2]-b[2]) < L_z/2.:
        o = a[2]-b[2]
    if a[2]-b[2] >= L_z/2.:
        o = a[2] - L_z -b[2]
    if a[2]-b[2] <= -L_z/2.:
        o = a[2] + L_z -b[2]

    v_rect = [q,n,o] 
    return [v_rect[0],v_rect[1],v_rect[2]]

def modulo(vector):
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
