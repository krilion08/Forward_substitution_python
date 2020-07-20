# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy
#import matplotlib.pyplot as plt


#solve diagonal matrix LGS
def diagonal(A, b):
    for i in range(0,len(A)):
        x[i]=b[i]/(A[i][i])
    return (b)

#upper triangular
def tri_up(A, b):
    x[-1]=b[-1]/A[-1][-1]
    for i in range(1,len(A)):
        buf=0
        for j in range(0,i):
            buf=buf+A[-1-i][-1-j]*x[-1-j]
        x[-1-i]=(b[-1-i]-buf)/A[-1-i][-1-i]
    return(x)
    
#lower triangular 
def tri_low(A,b):
    x[0]=b[0]/A[0][0]
    for i in range(1, len(A)):
        buf=0
        for j in range(0,i):
            buf=buf+A[i][j]*x[j]
        x[i]=(b[i]-buf)/A[i][i]
    return(x)


x=np.zeros(3)
b=np.zeros(3)
A=np.zeros([3, 3])

b[0]=27*np.random.rand()
b[1]=3*np.random.rand()
b[2]=3378*np.random.rand()


A[2][2]=np.random.rand()
A[1][1]=np.random.rand()
A[1][2]=np.random.rand()
A[0][2]=np.random.rand()
A[0][1]=np.random.rand()
A[0][0]=np.random.rand()

k=np.linalg.solve(A,b)

x=tri_up(A, b)

y=x-k