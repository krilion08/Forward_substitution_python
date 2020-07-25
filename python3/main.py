#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 23:52:30 2020

@author: MSeregin
"""

import sys
import numpy as np
import argparse

import LUP_Gaussian


parser=argparse.ArgumentParser(description='Solve system of linear equations Ax=b using the gaussian elimination method with partial pivoting')
parser.add_argument('matrix', metavar='A', type=str,nargs=1,
                    help='Matrix with coefficients')
parser.add_argument('vector', metavar='b', type=str,nargs=1,
                    help='Right-hand side vector')
parser.add_argument('-o','--output', metavar='x', type=str,nargs=1,
                    help='Soultion output file x')
parser.add_argument('-n','--nopiv', help='disable partial pivoting', default='False', action=('store_true'))
parser.add_argument('-d','--delimiter', metavar='del', type=str, nargs=1, help='set the delimiter', default=' ')

args=parser.parse_args()



#load matrix from stdin file
A = np.loadtxt(args.matrix[0], delimiter=args.delimiter[0])
b = np.loadtxt(args.vector[0], delimiter=args.delimiter[0])
print('A: \n',A)
print('b: \n',b)
x=np.zeros([len(A),1])
x=LUP_Gaussian.LUP(A, b, x)
print ('x: \n',x)