# Cálculo de Diagramas de interacción para verificación de la capacidad de las secciones NSR-10.

import math
import numpy as np
import matplotlib.pyplot as plt

## INPUT SECCIÓN ==============
# Dimensiones sección S3
b =  1000 #[mm]
h =  500 #[mm]

d4 = 50      # [mm]
d3 = 190     # [mm]
d2 = 500-190 # [mm]
d1 = 500-50  # [mm]


# Áreas de barras =====
N6 = 284 # [mm]^2
N7 = 387 # [mm]^2
N8 = 510 # [mm]^2 
# =====================
As_1 = 6*N7
As_2 = 2*N7
As_3 = 2*N7
As_4 = 6*N7


d = [d1, d2, d3, d4]
As= [As_1, As_2, As_3, As_4]

## INPUT MATERIAL ==============
fc = 35  # [MPa] f'c
fy = 420 # [MPa]

Es = 200e3 # [MPa]

phi = 0.75 # Estribos
#phi = 0.80 # Espiral

eps_c = 0.003  # C.10.2.3. Máxima deformación unitaria concreto en compresión.
eps_y = fy/Es  # C.10.3.2. Máxima unitaria acero en fluencia (approx 0.002).

# C.10.2.7.3 ==================
if fc <= 28:
    beta1 = 0.85
elif fc <= 35:
    beta1 = 0.80
elif fc <= 42:
    beta1 = 0.75
elif fc <= 48:
    beta1 = 0.70
else:
    beta1 = 0.65
    
#print(beta1)
# =============================


Ag  = b*h         #[mm]^2 Área bruta de la sección 
Ast = np.sum(As)  #[mm]^2 Área total acero 


# =============================
P0 = 0.85*fc*(Ag-Ast) + fy*Ast #[N] -> Capacidad axial máxima

A = [0,P0] #Pure axial compresion case. 



def neutralAx_depth(d1, Z, eps_y, eps_c): #Profundidad del eje neutro.
    return (eps_c/(eps_c-Z*eps_y))*d1   



# ---------------------------- Procedimiento iterativo -------------------------------------

Pn= []
Mn= []

findD = False # Auxiliar variable

for Z in np.arange(-10, 0.40, 0.01):
#for Z in [0.5, 0.25, 0, -0.5, -1, -2.5, -4, -6]:
    c = neutralAx_depth(d[0], Z, eps_y, eps_c)
    a = beta1*c # Profundidad del bloque de compresión de Whitney  

    #print(c)

    eps_s = [(( c -d[i])/c)*eps_c for i in range(len(d))]

    # Cálculo de tensiones en el acero
    fs = [np.sign(Es*eps_s[i])*min(abs(Es*eps_s[i]),fy) for i in range(len(eps_s))]

    # Cálculo de fuerza en el concreto
    Cc = (0.85*fc*a*b)

    # Cálculo de fuerzas en el acero
    Fs = [ fs[i]*As[i] if a<d[0] else (fs[i]-0.85*fc)*As[i] for i in range(len(fs))]

    # Par de fuerza y momento [Pn,Mn]
    Pn.append((Cc + np.sum(Fs)))
    Mn.append((Cc*(h-a)/2 + np.sum( [Fs[i]*(h/2 - d[i]) for i in range(len(Fs))])))

    if np.isclose(Z, 0, rtol=1e-08, atol=1e-08, equal_nan=False)  : B = [Mn[-1],Pn[-1]] # Zero stress in tension reinforcement
    if np.isclose(Z,-1, rtol=1e-08, atol=1e-08, equal_nan=False)  : C = [Mn[-1],Pn[-1]] # Balanced failure and compression-controlled limit

    if ~findD:
        if np.isclose(eps_s[0], -2.5*eps_y, rtol=1e-05, atol=1e-05, equal_nan=False) : 
            D = [Mn[-1],Pn[-1]] # Tension controlled limit
            findD = True
            #print("Found D")
        
F = [0, -Ast*fy] #Pure axial tension case. 

# ------------------------------------------------------------------------------------------
#for Z in [0.5, 0.25, 0, -0.5, -1, -2.5, -4, -6]:


# Diagrama sin penalizar ------------------------------
fig, ax = plt.subplots()


for i in range(len(Mn)-1):
    ax.plot([Mn[i],Mn[i+1]],[Pn[i],Pn[i+1]],'b') 

ax.plot([Mn[-1],0],[Pn[-1],A[1]],'b')
ax.plot([Mn[ 0],0],[Pn[ 0],F[1]],'b')


ax.plot([A[0],B[0],C[0],D[0],F[0]],[A[1],B[1],C[1],D[1],F[1]], 'ob')

# -----------------------------------------------------

# Diagrama de diseño ----------------------------------
phiMn = phi*np.array(Mn)
phiPn = phi*np.array(Pn)

for i in range(len(Mn)-1):
    ax.plot([phiMn[i],phiMn[i+1]],[phiPn[i],phiPn[i+1]],'k')

ax.plot([phiMn[-1],0],[phiPn[-1],phi*A[1]],'k')
ax.plot([phiMn[ 0],0],[phiPn[ 0],phi*F[1]],'k')


pA = phi*np.array(A)
pB = phi*np.array(B)
pC = phi*np.array(C)
pD = phi*np.array(D)
pF = phi*np.array(F)

ax.plot([pA[0],pB[0],pC[0],pD[0],pF[0]],[pA[1],pB[1],pC[1],pD[1],pF[1]], 'ok')
# -----------------------------------------------------

# Líneas auxiliares ----------------------------------
ax.plot([0,C[0]],[phi*pA[1],phi*pA[1]] , '-r', linewidth=1)

ax.plot([0,C[0]],[0,0] , '-k', linewidth=0.8)

ax.grid(color='r', linestyle='-', linewidth=0.1)
ax.set_ylabel('P')
ax.set_xlabel('M')

#ax.plot([Mn[len(Mn)-1],0],[Pn[len(Pn)-1],P0])


plt.show()


#print(eps_s) 
#print(fs) 
#print(Fs)