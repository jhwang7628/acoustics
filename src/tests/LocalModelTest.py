#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt
import math
import random
from math import pi

# def h(xi, omega, omega_d, m, t): 
#     return np.exp(-xi*omega*t) / (m*omega_d) * np.sin(omega_d*t)

def epsilon(xi, omega, h): 
    return np.exp(-xi*omega*h) 

def theta(omega_d, h): 
    return omega_d * h

def gamma(xi): 
    return math.asin(xi)

def Compute_coeff_qNew(epsilon, theta): 
    return 2.0*epsilon*np.cos(theta)

def Compute_coeff_qOld(epsilon): 
    return -epsilon*epsilon

def Compute_coeff_Q(epsilon, theta, gamma, omega, omega_d, m): 
    return 2.0/(3.0*omega*omega_d*m) * (epsilon*np.cos(theta+gamma) - epsilon**2*np.cos(2.0*theta+gamma))

def IIR(qOld, qNew, coeff_qOld, coeff_qNew, coeff_Q, Q): 
    return coeff_qNew*qNew + coeff_qOld*qOld + coeff_Q*Q

def EvaluateModel(xi, omega, omega_d, t0, t, q, qDot, flag): 
    if flag == 0: 
        return np.exp(-xi*omega*(t-t0)) * (q*np.cos(omega_d*(t-t0)) + (qDot+xi*omega*q)/omega_d * np.sin(omega_d*(t-t0)))
    else:  # polynomial approximation
        a = xi*omega*(t-t0)
        b = q
        c = (qDot+xi*omega*q)/omega_d
        x = omega_d * (t-t0)
        return (1.0 - a + a**2/2.0)*(b*(1.0 - x**2/2.0) + c*x)

if __name__ == "__main__": 
    alpha = 40.0
    beta = 1.0E-7
    m = 2300.0 
    f = 500.0
    omega = 2.0*pi*f
    xi = 0.5*(alpha/omega + beta*omega) 
    omega_d = omega * np.sqrt(1.0 - xi**2) 
    h = 1./176400.

    epsilon = epsilon(xi, omega, h)
    theta = theta(omega_d, h)
    gamma = gamma(xi) 

    coeff_qNew = Compute_coeff_qNew(epsilon, theta) 
    coeff_qOld = Compute_coeff_qOld(epsilon) 
    coeff_Q    = Compute_coeff_Q(epsilon, theta, gamma, omega, omega_d, m) 

    N = 1.0 / h + 1
    T = np.linspace(0, 1, N, endpoint=True)

    # build random model around t0
    modelStartTick = random.randint(0, 100)
    t0 = float(modelStartTick)*h
    Tmodel = np.linspace(t0, 1., len(T)-modelStartTick)
    qmodel = np.zeros(len(Tmodel))
    qmodel2= np.zeros(len(Tmodel))

    qOld = 0.
    qNew = 0.
    allq = np.zeros(len(T))
    count = 0
    countModel = 0
    for t in T: 
        if t < 1./40000.: 
            Q = 1.
        else: 
            Q = 0.
        q = IIR(qOld, qNew, coeff_qOld, coeff_qNew, coeff_Q, Q)
        qOld = qNew
        qNew = q
        allq[count] = q

        # build and evaluate model
        if (count >= modelStartTick): 
            q = allq[modelStartTick]
            # qDot = (allq[modelStartTick] - allq[modelStartTick-1]) / h  # first order
            qDot = (1.5*allq[modelStartTick] - 2.0*allq[modelStartTick-1] + 0.5*allq[modelStartTick-2])/h # second order
            qmodel[countModel]  = EvaluateModel(xi, omega, omega_d, t0, t, q, qDot, 0)
            qmodel2[countModel] = EvaluateModel(xi, omega, omega_d, t0, t, q, qDot, 1)
            countModel += 1 

        count += 1

    draw = range(1000)
    plt.figure() 
    plt.plot(T[draw], allq[draw], Tmodel[draw], qmodel[draw], 'rx', Tmodel[draw], qmodel2[draw], 'go') 
    plt.show()

    

