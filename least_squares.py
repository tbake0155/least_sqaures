#!/usr/bin/env python3
#
# Author:     Timothy Baker
# Date:       October 11, 2017
# Class:      CS-517
# Assignment: Machine Assignment 2

import sys
import sympy as sym
import numpy as num

def t(i, N) :
    return -1.0 + ((2.0 * i)/N)

def getT(N, T) :
    for i in range (0,int(N)+2) :
        T.append(float(t(i,N)))    
    return T

def betaKay(k, N) :
    if k == 0 :
        return 2
    else :
        return ((1.0+1.0/N)**2) * (1.0-((k/(N+1.0))**2)) * ((4.0-(1.0/(k*k)))**(-1))

def getBk(N, Bk) :
    for i in range (0,int(N)+2) :
        Bk.append(float(betaKay(i,N)))
    return Bk

def piKay(k, N) :
    t = sym.Symbol('t')
    if k == 0 :
        return 1.0
    elif k == -1.0 :
        return 0
    else :
        return (t)*piKay(k-1, N) - betaKay(k, N) * piKay(k-2, N)

def getPk(N, Pk) :
    for i in range(0, int(N)+2) :
        Pk.append(piKay(i, N))
    return Pk

def normL2 (Pii, Pij, T, N) :
    Fx = sym.Mul(Pii, Pij)
    t = sym.Symbol('t')
    summed = 0
    for k in range(0, len(T)) :
        summed = summed + Fx.subs(t, T[k])
    return summed*(2/N+1)

def normL (Pii, Pij, T, N) :
    Fx = Pii * Pij
    summed = 0
    for k in range(0, len(T)) :
        summed = summed + Fx
    return summed*(2/N+1)

def getNorms(N, T, Pk, Yk) :
    for i in range(0, int(N)+2) :
        Yk.append(normL2(Pk[i], Pk[i], T, N))
    return Yk

def displayPartOne(N, Bk, Yk, Mk) :
    print("\n\n")
    print("  ", '%-5s' % 'K', '%-11s' % "Bk", '%-11s' % "Yk", '%-11s' % "Mk+1", "\n")
    for i in range(0,int(N)+1) :
        print("  ", '%-5d' % i,'%-11f' % Bk[i],'%-11f' % Yk[i], '%-11f' % Mk[i+1])

def getC(i,j, Ck) :
    C = Ck[i][0][j]
    C = C.astype(float)
    return C

def displayPartTwo(N, FX, Ck, Error, infError) :
    for i in range(0,len(FX)):
        print("\n\n")
        print("  ", '%-5s' % FX[i])
        print("  ", '%-5s' % 'n', '%-11s' % "Cn", '%-11s' % "||en||2", '%-11s' % "||en||inf", "\n")
        for j in range(0,int(N)+1) :
            print("  ", '%-5d' % j,'%-11f' % getC(i,j, Ck),'%-11f' % Error[i][j], '%-11f' % 0.0)

def fOFt(N, Pk, T, Mk) :
    for i in range (0, int(N)+2) :
        tempList = []
        for j in range (0, int(N)+1) :
            if type(Pk[i]) == float :
                tempList.append(Pk[i])
            else :
                tempList.append(Pk[i].subs(sym.Symbol('t'), T[j]))
        Mk.append(max(tempList))
    return Mk

def getN() :
    if len(sys.argv) > 1 :
        return float(sys.argv[1])
    else :
        return 10.0

def getFX() :
    t = sym.Symbol('t')
    FX =[(sym.E**(-t)), (sym.ln(2+t)), (t+1)**sym.Rational(1,2), (t**2)**(sym.Rational(1,2))]
    return FX

def errorL2(N, T, Fx, B) :
    EL2 = []
    for j in range(0, len(Fx)) :
        normList = []
        for i in range(0, int(N)+1) :
            fx = Fx[j].subs(sym.Symbol('t'), i)
            Bi = B[j]
            normList.append(normL(Bi[i], fx, T, N))
        EL2.append(normList)
    return EL2

def linearize(N, A, B) :
    Ck = []
    for i in range(0, len(B)) :
        A = num.array(A)
        A = A.astype(float)
        Bi = B[i]
        Bi = num.array(Bi)
        Bi = Bi.astype(float)
        C = num.linalg.lstsq(A,Bi,rcond=-1)
        Ck.append(C)
    return Ck

def getB(N, T, FX, Pk) :
    Bs =[]
    for i in range(0, len(FX)) :
        tempList = []
        for j in range(0, int(N)+1) :
            tempList.append(normL2(FX[i], Pk[j], T, N))
        Bs.append(tempList)
    return Bs

def getA(Yk) :
    A = [[Yk[0],0,0,0,0,0,0,0,0,0,0],
        [0,Yk[1],0,0,0,0,0,0,0,0,0],
        [0,0,Yk[2],0,0,0,0,0,0,0,0],
        [0,0,0,Yk[3],0,0,0,0,0,0,0],
        [0,0,0,0,Yk[4],0,0,0,0,0,0],
        [0,0,0,0,0,Yk[5],0,0,0,0,0],
        [0,0,0,0,0,0,Yk[6],0,0,0,0],
        [0,0,0,0,0,0,0,Yk[7],0,0,0],
        [0,0,0,0,0,0,0,0,Yk[8],0,0],
        [0,0,0,0,0,0,0,0,0,Yk[9],0],
        [0,0,0,0,0,0,0,0,0,0,Yk[10]]]
    return A

def errorInf(N, T, FX, B) :
    return 0.0

def main() :
    N = getN()
    T  = getT(N, [])
    Bk = getBk(N, [])
    Pk = getPk(N, [])
    Yk = getNorms(N, T, Pk, [])
    Mk = fOFt(N, Pk, T, [])
    displayPartOne(N, Bk, Yk, Mk)

    FX = getFX()
    B = getB(N, T, FX, Pk)
    A = getA(Yk)
    Ck = linearize(N, A, B)
    Error = errorL2(N, T, FX, B)
    ErrorInf = errorInf(N, T, FX, B)
    displayPartTwo(N, FX, Ck, Error, ErrorInf)
    print("\n\n")

if __name__ == "__main__":
    main()

#end of file 

