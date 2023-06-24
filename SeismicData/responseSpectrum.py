# Cálculo del espectro de respuesta para un sismo.

import math
import matplotlib.pyplot as plt

import FunctionUtils as F

from ReadUtils import importAccelerogram

from Dynamic_Methods import u_exact_linealExitation as u
from Dynamic_Methods import v_exact_linealExitation as v
from Dynamic_Methods import a_exact_linealExitation as a

# Method: Piecewise exact solution (e.g. Clough 1999)

# INUPT DATA
#file = 'Accelerograms/Kobe.txt'
file = 'Accelerograms/Elcentro.dat'
#file = 'Accelerograms/Northridge.txt'
#file = 'Accelerograms/TestFree_ElCentro.dat'
#========================================================================================

# Spectrum computations ============================================================================

DATA = importAccelerogram(file,9.80665)

time = DATA[:][0] # Time discretization corresponding to Ag
ag   = DATA[:][1] # Ground Acceleration [m]/[s^2]

n = len(time)

xi = 0.05

Ti = 0 # [s]
Tf = 10 # [s]

step = 0.1
Tf += step
n_steps = int((Tf-Ti)/step)

T_array = [Ti+t*step for t in range(0,n_steps)]

if T_array[0]==0:T_array.pop(0)

dis_Spect = []
vel_Spect = []
acc_Spect = []

for T in T_array:
    # SOLUTION FOR SPECIFIC PERIOD ===================================

    w=(2*math.pi)/T #Natural frequency

    u_t = [] #To store the solution of each time interval
    v_t = []
    a_t = []

    u0 = 0
    v0 = 0

    global_max_u = 0
    global_max_v = 0
    global_max_a = 0

    N_points = 10

    for i in range(n-1):

        dt = time[i+1]-time[i]
        u_t.append(u(xi,w,u0,v0,ag[i],ag[i+1],dt)) #OK
        v_t.append(v(xi,w,u0,v0,ag[i],ag[i+1],dt))
        a_t.append(a(xi,w,u0,v0,ag[i],ag[i+1],dt)) #OK

        u0 = u_t[i](dt)
        local_max_v = max([v0,(u_t[i](dt)-u_t[i](0))/dt])

        #v0 = (u_t[i](dt)-u_t[i](0))/dt

        #print(v0,"---",v_t[i](dt),"---",v0-v_t[i](dt))
        v0 = v_t[i](dt)

        local_max_u = F.abs_maxF(u_t[i],0,dt,N_points)
        if global_max_u < local_max_u: global_max_u = local_max_u

        
        local_max_v = F.abs_maxF(v_t[i],0,dt,N_points)
        if global_max_v < local_max_v: global_max_v = local_max_v

        local_max_a = F.abs_maxF(a_t[i],0,dt,N_points)
        if global_max_a < local_max_a: global_max_a = local_max_a

        #print(local_max_u)

    print(T)

    dis_Spect.append(global_max_u)
    vel_Spect.append(global_max_v)
    acc_Spect.append(global_max_a)

# =====================================================================
# PLOT SPECTRA
fig, ax = plt.subplots(2,2)

ax[0,0].set_title('Acelerograma')
ax[0,1].set_title('Espectro de Desplazamiento')
ax[1,0].set_title('Espectro de Velocidad')
ax[1,1].set_title('Espectro de Aceleración')

ax[0,0].plot(time,ag) # Sismo. Acelerograma
ax[0,1].plot(T_array,dis_Spect) # Espectro de Desplazamientos.
ax[1,0].plot(T_array,vel_Spect) # Espectro de Velocidades.
ax[1,1].plot(T_array,acc_Spect) # Espectro de Aceleración.

ax[0,0].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax[0,1].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax[1,0].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax[1,1].grid(color="gray", which="both", linestyle=':', linewidth=0.5)

# =====================================================================
## PLOT LAST RESPONSE.

fig, ax2 = plt.subplots(2,2)

ax2[0,1].set_title('Desplazamiento')
ax2[0,0].set_title('Acelerograma')
ax2[1,0].set_title('Velocidad')
ax2[1,1].set_title('Aceleración')

ax2[0,0].plot(time,ag) # Sismo. Acelerograma
ax2[0,1] = F.plot_picewise_lambda_localX(ax2[0,1],u_t,time,100)
ax2[1,0] = F.plot_picewise_lambda_localX(ax2[1,0],v_t,time,100)
ax2[1,1] = F.plot_picewise_lambda_localX(ax2[1,1],a_t,time,100)

ax2[0,0].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax2[0,1].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax2[1,0].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax2[1,1].grid(color="gray", which="both", linestyle=':', linewidth=0.5)

# =====================================================================

plt.show()


    

    
























