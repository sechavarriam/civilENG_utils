# Cálculo del espectro de respuesta para un sismo.
#
# Method: Piecewise exact solution (e.g. Clough 1999)

import math
import matplotlib.pyplot as plt

from ReadUtils import importAccelerogram

# INUPT DATA
#file = 'Accelerograms/Kobe.txt'
file = 'Accelerograms/Elcentro.dat'
#file = 'Accelerograms/Northridge.txt'

DATA = importAccelerogram(file)


time = DATA[:][0] # Time discretization corresponding to Ag
ag   = DATA[:][1] # Ground Acceleration                     



n = len(time)
#========================================================================================




# Usefull function definitions =====================================================================

# Clough Eq 7.5.
def u(xi,w,u0,v0,P0,P1,h):
    
    alpha = (P1-P0)/h
    wD = w*math.sqrt(1-xi**2) #Damped Frequency

    A0 = P0/(w**2) - 2*xi*alpha/(w**3)
    A1 = alpha/(w**2)
    A2 = u0-A0
    A3 = (v0+(xi*w*A2)-A1)/wD
    
    return lambda t: A0+A1*t+(A2*math.cos(wD*t)+A3*math.sin(wD*t))*math.exp(-xi*w*t)

def v(xi,w,u0,v0,P0,P1,h):
    alpha = (P1-P0)/h
    wD = w*math.sqrt(1-xi**2) #Damped Frequency

    A0 = P0/(w**2) - 2*xi*alpha/(w**3)
    A1 = alpha/(w**2)
    A2 = u0-A0
    A3 = (v0+ (xi*w*A2) -A1)/wD

    B2 =  (wD*A3-xi*w*A2)
    B3 = -(wD*A2+xi*w*A3)

    return lambda t: A1*t+(B2*math.cos(wD*t)+B3*math.sin(wD*t))*math.exp(-xi*w*t)

def a(xi,w,u0,v0,P0,P1,h):
    
    alpha = (P1-P0)/h
    wD = w*math.sqrt(1-xi**2) #Damped Frequency

    A0 = P0/(w**2) - 2*xi*alpha/(w**3)
    A1 = alpha/(w**2)
    A2 = u0-A0
    A3 = (v0+ (xi*w*A2) -A1)/wD

    B2 =  (wD*A3-xi*w*A2)
    B3 = -(wD*A2+xi*w*A3)

    C2 =  (wD*B3-xi*w*B2)
    C3 = -(wD*B2+xi*w*B3)

    return lambda t: (C2*math.cos(wD*t)+C3*math.sin(wD*t))*math.exp(-xi*w*t)
    #C1 = -xi*w*(wD*A3-xi*w*A2)- wD*(wD*A2+xi*w*A3)
    #C2 = -wD*(wD*A3-xi*w*A2) + xi*w*(wD*A2+xi*w*A3)
    #return lambda t: C1*math.exp(-xi*w*t)*math.cos(wD*t)+ C2*math.exp(-xi*w*t)*math.sin(wD*t)

alpha = lambda P0,P1,h: (P1-P0)/h


def abs_maxF(f,x0,x1,n_points):
    # f is a lambda scalar valued function defined in [x0,x1]

    step = (x1-x0)/n_points
    x_eval = (xi*step for xi in range(0,n_points+1)) #Evaluation test points
    #x_eval = [xi*step for xi in range(0,n_points+1)] #Evaluation test points
    
    #for x in x_eval:
    #    print(x,f(x)) 

    #print("maximo------> ",max([abs(f(x)) for x in x_eval]) )
    #print("===============================")

    return max([abs(f(x)) for x in x_eval])
    

def fplot(ax,f,x0,x1,translation=0,n=5):
    # plots a lambda real function f in a defined interval [x0,x1] in n line segments, in a given axis object.
    
    step = (x1-x0)/n
    x_eval = [x0+xi*step for xi in range(0,n+1)] #Evaluation test points
    f_eval = [f(x) for x in x_eval]

    x_plot = [x+translation for x in x_eval]


    ax.plot(x_plot,f_eval,'b') 

    return ax

def plot_picewise_lambda(ax,F,X,N=5):
    n = len(F) # Extract interval number (lenX)=n+1)

    for t in range(n):
        ax = fplot(ax,F[t],X[t],X[t+1],0, N)

    return ax
     
def plot_picewise_lambda_localX(ax,F,X, N=5):
    n = len(F) # Extract interval number (lenX)=n+1)

    
    for t in range(n):
        ax = fplot(ax,F[t],0,X[t+1]-X[t],X[t],N)

    return ax

# Spectrum computations ============================================================================

xi = 0.05

Ti = 0 # [s]
Tf = 10 # [s]

step = 0.1
n_steps = int((Tf-Ti)/step)

T_array = [Ti+t*step for t in range(0,n_steps)]
T_array.pop(0)

dis_Spect = []
vel_Spect = []
acc_Spect = []

# Solution for Fixed T

for T in T_array:
    # SOLUTION FOR SPECIFIC PERIOD ===================================
    w=(2*math.pi)/T

    u_t = [] #To store the solution of each time interval
    v_t = []
    a_t = []

    u0 = 0
    v0 = 0

    global_max_u = 0
    global_max_v = 0
    global_max_a = 0

    N_points = 3

    for i in range(n-1):
        dt = time[i+1]-time[i]
        u_t.append(u(xi,w,u0,v0,ag[i],ag[i+1],dt))
        v_t.append(v(xi,w,u0,v0,ag[i],ag[i+1],dt))
        a_t.append(a(xi,w,u0,v0,ag[i],ag[i+1],dt))

        u0 = u_t[i](dt)
        v0 = v_t[i](dt)

        local_max_u = abs_maxF(u_t[i],0,dt,N_points)
        if global_max_u < local_max_u: global_max_u = local_max_u

        local_max_v = abs_maxF(v_t[i],0,dt,N_points)
        if global_max_v < local_max_v: global_max_v = local_max_v

        local_max_a = abs_maxF(a_t[i],0,dt,N_points)
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
# PLOT LAST RESPONSE.


fig, ax2 = plt.subplots(2,2)

ax2[0,1].set_title('Desplazamiento')
ax2[0,0].set_title('Acelerograma')
ax2[1,0].set_title('Velocidad')
ax2[1,1].set_title('Aceleración')

ax2[0,0].plot(time,ag) # Sismo. Acelerograma
ax2[0,1] = plot_picewise_lambda_localX(ax2[0,1],u_t,time,100)
ax2[1,0] = plot_picewise_lambda_localX(ax2[1,0],v_t,time,100)
ax2[1,1] = plot_picewise_lambda_localX(ax2[1,1],a_t,time,100)

ax2[0,0].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax2[0,1].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax2[1,0].grid(color="gray", which="both", linestyle=':', linewidth=0.5)
ax2[1,1].grid(color="gray", which="both", linestyle=':', linewidth=0.5)

#print(global_max_a)
#ax = plot_picewise_lambda_localX(ax,a_t,time,10)

#ax.plot(time,ag)


    # Max evaluat. Use ternary operator

#print(acc_Spect)

## TEST ==========================================
#f_test_1 =  lambda x: -(x-1)**2
#f_test_2 =  lambda t: t


#F = [f_test_1,f_test_2]
#t = [0,1.5,2]

#ax = plot_picewise_lambda(ax,F,t,100)
#ax = plot_picewise_lambda_localX(ax,F,t,100)

plt.show()


    

    
























