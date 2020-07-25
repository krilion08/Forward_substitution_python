#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 22:59:25 2020

LUP decomposition of a matrix using the gaussian elimination method with partial pivoting

@author: MSeregin
"""

import numpy as np



#solve diagonal matrix LGS
def diagonal(A, b, x):
    for i in range(0,len(A)):
        x[i]=b[i]/(A[i][i])
    return (b)

#upper triangular
def tri_up(A, b, x):
    x[-1]=b[-1]/A[-1][-1]
    for i in range(1,len(A)):
        buf=0
        for j in range(0,i):
            buf=buf+A[-1-i][-1-j]*x[-1-j]
        x[-1-i]=(b[-1-i]-buf)/A[-1-i][-1-i]
    return(x)
    
#lower triangular 
def tri_low(A, b, x):
    x[0]=b[0]/A[0][0]
    for i in range(1, len(A)):
        buf=0
        for j in range(0,i):
            buf=buf+A[i][j]*x[j]
        x[i]=(b[i]-buf)/A[i][i]
    return(x)


#LUP decomposition
def LUP(U, b, x):
    N=len(U)
    P=np.identity(N,dtype='d')
    L=np.identity(N,dtype='d')
    #loop over all rows
    for i in range(0,N-1):
    #pivoting
        imax=i
        for j in range(i, N):
            if abs(U[imax][i])<abs(U[j][i]):
                imax=j
        #change rows if needed
        if imax!=i:
            #"elegant"
            U[[imax,i]]=U[[i,imax]]
            P[[imax,i]]=P[[i,imax]]
            #change the rows so, that the diagonal of L were not changed
            for k in range(0, i):
                buf=L[i][k]
                L[i][k]=L[imax][k]
                L[imax][k]=buf
    #gaussian elimination
        for k in range(i,N-1):
            #division by pivot -> element of L
            L[k+1][i]=U[k+1][i]/U[i][i]
            for j in range(i, N):
                U[k+1][j]=U[k+1][j]-L[k+1][i]*U[i][j]


    #solve lover triangular system
    x=tri_low(L, np.dot(P,b), x)
    #solve upper triangular system
    x=tri_up(U, x, x)
    return(x)