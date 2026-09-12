#Creation of vertices
pos0 = []
vel0 = []
cells0 = []
t = []

for i in range(n_vertices):
    pos0.append([x[i],y[i],0])
    vel0.append((0,0,0))

    if tipo_v[i]==1:
        l = 0
        for j in range(len(vtipo1)):
            if vtipo1[j] == i and l==0:
                v1 = j
                cells_i_horario = cells_tipo1[v1]
                ady_i_horario = [[adys_tipo1[v1][0],adys_tipo1[v1][1]],[adys_tipo1[v1][2],adys_tipo1[v1][3]],[adys_tipo1[v1][4],adys_tipo1[v1][5]]]
                t.append(Vertex(i, pos0[i],vel0[i],np.array(cells_i_horario), np.array(ady_i_horario), dt, [-1,-1], [-1,0,0]))
                l = l+1


#Creations of cells
celulas = []
for i in range(n_celulas):
    b = []
    for j in range(len(celda)):
        if celda[j] == i:
            b.append(int(v_celda[j]))
    celulas.append(Cell(i,b, dt, np.random.normal(1,0.1), KA, KP))  #area00[i], KA, KP))   #b tiene repetido b[inicial] b[final]

    b.append(b[0])


#Creation of tissue
T1 = Tissue(t, celulas, dt)

lista_i = []
lista_j = []
lij_save = []
flats = []
for i in range(n_vertices):
    flat = [T1.vs[i].ady[0][0],T1.vs[i].ady[0][1],T1.vs[i].ady[1][0]]
    flats.append(flat)
for i in range(n_vertices):
    for j in flats[i]:
        if j>i and j in np.array(T1.vs[i].ady).flatten():
            lista_i.append(i)
            lista_j.append(j)
            lij_save.append(modulo(resta_vectores(T1.vs[i].r,T1.vs[j].r)))
tensiones_distribution = np.random.normal(0.1,0.01,len(lista_i))


LargosNat = []
LargosNat00 = []
for i in range(n_vertices):
    LargosNat.append([])
    LargosNat00.append([])
    for j in range(n_vertices):
        LargosNat[i].append(0)
        LargosNat00[i].append(0)

for i in range(len(lista_i)):
    v1 = int(lista_i[i])
    v2 = int(lista_j[i])
    LargosNat[v1][v2] = lij_save[i]
    LargosNat[v2][v1] = lij_save[i]
    LargosNat00[v1][v2] = lij_save[i]
    LargosNat00[v2][v1] = lij_save[i]

Tensiones = []
Tensiones4fold = []
TensionesFluct = []
for i in range(n_vertices):
    Tensiones.append([])
    Tensiones4fold.append([])
    TensionesFluct.append([])
    for j in range(n_vertices):
        Tensiones[i].append(0)
        Tensiones4fold[i].append(0)
        TensionesFluct[i].append(0)

for i in range(len(lista_i)):
    v1 = int(lista_i[i])
    v2 = int(lista_j[i])
    Tensiones[v1][v2] = tensiones_distribution[i]
    Tensiones[v2][v1] = tensiones_distribution[i]

for i in range(n_vertices):
    for j in range(n_vertices):
        if Tensiones[i][j]<0:
            Tensiones[i][j]=0


Tact = []
for i in range(n_vertices):
    Tact.append([])
    for j in range(n_vertices):
        Tact[i].append(Gamma_act)

cs_moviles = range(n_celulas)