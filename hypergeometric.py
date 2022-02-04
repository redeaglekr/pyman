import numpy as np
from mpmath import *
from scipy.special import *
def Ns(z,p,n):## to find N##
    a=(z-p)**2
    b=(z+p)**2
    k=(1-a)**((n+6)/4)
    bet=beta((n+8)/4,(1/2))
    ist=(((z**2)-(p**2))/a)*appellf1(0.5,1,0.5,(n+10)/4,(a-1)/a,(1-a)/(b-a))
    snd=hyp2f1(0.5,0.5,(n+10)/4,(1-a)/(b-a))
    N=((k*bet)*(ist-snd))/((b-a)**0.5)
    return N
def Ms(z,p,n):##to find M##
    a=(z-p)**2
    b=(z+p)**2
    k=(1-a)**((n+4)/4)
    ist1=(((z**2)-(p**2))/a)*appellf1(0.5,-(n+4)/4,1,1,(b-a)/(1-a),(a-b)/a)
    snd1=hyp2f1(-(n+4)/4,0.5,1,(b-a)/(1-a))
    M=k*(ist1-snd1)
    return M
def elik(k):
    m1=1.-k**2
    logm1 = log(m1)

    a1=0.44325141463
    a2=0.06260601220
    a3=0.04757383546
    a4=0.01736506451
    b1=0.24998368310
    b2=0.09200180037
    b3=0.04069697526
    b4=0.00526449639
    ee1=1.+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*(-logm1)
    ek = ee1+ee2

    a0=1.38629436112
    a1=0.09666344259
    a2=0.03590092383
    a3=0.03742563713
    a4=0.01451196212
    b0=0.5
    b1=0.12498593597
    b2=0.06880248576
    b3=0.03328355346
    b4=0.00441787012
    ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*logm1
    kk = ek1-ek2
    return [ek,kk]


##Note that  the Appells  and Gauss Hypergeometric functions are calculated using Sympy
    




   
    
