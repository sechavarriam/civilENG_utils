# This file contain a collection of methods and functions used in structural dynamics.

import math
import matplotlib.pyplot as plt


# Usefull function definitions ==========================================================

# Homogeneous solution for damped system (Free response)
def u_h(xi,w,u0,v0):
    wD = w*math.sqrt(1-xi**2)

    return lambda t: math.exp(-xi*w*t)*(u0*math.cos(wD*t) + ((v0-u0*xi*w)/wD)*math.sin(wD*t))

# =======================================================================================

def u_exact_linealExitation(xi,w,u0,v0,P0,P1,dt):
    
    alpha = (P1-P0)/dt
    wD = w*math.sqrt(1-xi**2) #Damped Frequency

    A0 = P0/(w**2) - (2*xi*alpha)/(w**3) #[m]
    A1 = alpha/(w**2)                    #[m]/[s]
    A2 = u0-A0                           #[m]
    A3 = (v0+(xi*w*A2)-A1)/wD            #[m]
    
    return lambda t: A0 + A1*t + (A2*math.cos(wD*t)+A3*math.sin(wD*t))*math.exp(-xi*w*t)

def v_exact_linealExitation(xi,w,u0,v0,P0,P1,h):

    alpha = (P1-P0)/h
    wD = w*math.sqrt(1-xi**2) #Damped Frequency

    A0 = P0/(w**2) - (2*xi*alpha)/(w**3) #[m]
    A1 = alpha/(w**2)                    #[m]/[s]
    A2 = u0-A0                           #[m]
    A3 = (v0+(xi*w*A2)-A1)/wD            #[m]

    B2 =  (wD*A3-xi*w*A2)                #[m]/[s]
    B3 = -(wD*A2+xi*w*A3)                #[m]/[s]

    return lambda t: A1 + (B2*math.cos(wD*t)+B3*math.sin(wD*t))*math.exp(-xi*w*t)

def a_exact_linealExitation(xi,w,u0,v0,P0,P1,h):
    
    alpha = (P1-P0)/h

    wD = w*math.sqrt(1 - xi**2) #Damped Frequency

    A0 = P0/(w**2) - 2*xi*alpha/(w**3) #[m]
    A1 = alpha/(w**2)                  #[m]/[s]
    A2 = u0-A0                         #[m]
    A3 = (v0+ (xi*w*A2) -A1)/wD        #[m]

    B2 =  (wD*A3-xi*w*A2)              #[m]/[s]
    B3 = -(wD*A2+xi*w*A3)              #[m]/[s]

    C2 =  (wD*B3-xi*w*B2)              #[m]/[s]^2
    C3 = -(wD*B2+xi*w*B3)              #[m]/[s]^2

    return lambda t: (C2*math.cos(wD*t)+C3*math.sin(wD*t))*math.exp(-xi*w*t)

# =======================================================================================

# Response to a Dirac impulse. Impulse response function h(t)
# Paultre (eq. 7.17)
def IRF(xi,w,m=1):
    # m: mass (=1 by default)
    wD = w*math.sqrt(1-xi**2) #Damped Frequency

    return lambda t: (math.exp(-xi*w*t)*math.sin(wD*t))/(wD*m)

# Duhamel Integral









