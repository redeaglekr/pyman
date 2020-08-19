##This is a program written for plotting analytical(Theoretical) Light curve using Mandel and Agol
##This can accomodate Linear limb darkening
###This uses all 11 conditions of table 1 in the paper
import sys
sys.path.append("/home/sreenivasan/cusat/project/codes")
from hypergeometric import Ns,Ms
from keplermodel import Zlist
from math import *
from mpmath import *
import numpy as np
def mlal(time,c1,c2,c3,c4,rat,inc,per,T0,p0):
    z=Zlist(time,per,inc,T0,rat)
    c=[c1,c2,c3,c4]
    c1=c[0]
    c2=c[1]
    c3=c[2]
    c4=c[3]
    c0=1-c1-c2-c3-c4
    k0=c0/4
    m=[]
    ##flux array which is going to be appended by the calculation given below
    m.append(k0)
    for i in range(1,5,1):
        k=c[i-1]/(i+4)
        m.append(k)
    om=np.sum(m)
    f=[]
    for y in range(0,len(z),1):
        if p0>0 and z[y]>(0.5+abs(p0-0.5)) and z[y]<(1+p0):
            t=[]
            n0=Ns(z[y],p0,0)*k0
            t.append(n0)
            for fj in range(1,5,1):
                j=Ns(z[y],p0,i)*(c[fj-1]/(fj+4))
                t.append(j)
            summ=np.sum(t)
            fg=(1/(2*np.pi*om))*(summ)
            hg=float(1-fg)
            f.append(hg)
            t.clear()
        elif  (p0>0 and p0<0.5 and z[y]>p0 and z[y]<(1-p0)):
            lm=[]
            for th in range(1,4,1):
                j=Ms(z[y],p0,th)*(c[i-1]/(i+4))
                lm.append(j)
            summ2=np.sum(lm)
            l=((p0)**2)*((1-(p0)**2)/(2-(z**2)))
            fn=(c0*((p0)**2))+(2*summ2)+(c4*l)
            fn2=fn/(4*om)
            fnm=(1-fn2)
            f.append(fnm)
            lm.clear()
        elif  (p0>0 and p0<0.5 and z[y]==1-p0):
            lm=[]
            for i in range(1,4,1):
                j=Ms(z[y],p0,i)*(c[i-1]/(i+4))
                lm.append(j)
            summ2=np.sum(lm)
            l=((p0)**2)*((1-(p0)**2)/(2-(z[y]**2)))
            fn=(c0*((p0)**2))+(2*summ2)+(c4*l)
            fn3=fn/(4*om)
            fnm=float(1-fn3)
            f.append(fnm)
            lm.clear()
        elif p0>0 and p0<0.5 and z[y]==p0:
            m5=[]
            so=k0*hyp2f1(0.5,-1,1,4*(p0)**2)
            m5.append(so)
            for i in range(1,5,1):
                j=(c[i-1]/(i+4))*(hyyp2f1(0.5,-(i+4)/4,1,4*(p0)**2))
                m5.append(j)
            c5=np.sum(m5)
            t2=c5/(2*om)
            fv=float(0.5+t2)
            f.append(fv)
            m5.clear()
        elif p0 == 0.5 and z [y]== 0.5:
            nm=[]
            mo=k0*(gamma(1.5))/(gamma(2))
            nm.append(mo)
            for i in range(1,5,1):
                g=((c[i-1])/4)*gamma(1.5+(i/2))/gamma(2+(i/4))
                nm.append(g)
            df=np.sum(nm)   
            sm=nm/(2*np.sqrt(np.pi)*om)
            f6=float(0.5+sm)
            f.append(f6)
            nm.clear()
        elif p0>0.5 and z[y]==p0:
            df=[]
            bo=ko*(beta(0.5,2)*hyp2f1(0.5,0.5,2.5,(1/4*p0**2)))
            df.append(bo)
            for i in range(1,5,1):
                qw=(c[i-1]/(i+4))*(beta(0.5,(i+8)/4)*hyp2f1(0.5,0.5,2.5+(i/4),1/(4*p0**2)))
                df.append(qw)
            dfg=np.sum(df)
            fi=dfg/(4*p0*np.pi*om)
            flo=float(0.5+fi)
            f.append(flo)
            df.clear()
        elif p0>0.5 and z[y]>=abs(1-p0) and z[y]<p0:#case 8
            jh=[]
            df=Ns(z[y],p0,0)*k0
            jh.append(df)
            for i in range(1,5,1):
                j=((c[i-1])*Ns(z,p0,i))/(i+4)
                jh.append(j)
            a=np.sum(jh)
            ds=float(-(a)/(4*om))
            f.append(ds)
            jh.clear()
        elif p0>0 and p0<1 and z[y]>0 and z[y]<(0.5-abs(p0-0.5)):
            ed=[]
            for i in range(1,4,1):
                mf=c[i-1]/(i+4)*(Ms(z[y],p0,i))
                ed.append(mf)
            bn=np.sum(ed)
            l=((p0)**2)*((1-(p0)**2)/(2-z**2))
            op=(c0*(1-(p0)**2))+(c4*(0.5-l))-(2*bn)
            fml=float(op/(4*om))
            f.append(fml)
            ed.clear()
        elif p0>0 and p0<1 and z[y]==0:
            tf=[]
            t0=k0*(1-(p0)**2)
            tf.append(t0)
            for i in range(1,5,1):
                ti=(c[i-1])*((1-(p0)**2)**((i+4)/4))/(i+4)
                tf.append(ti)
            sd=np.sum(tf)
            de=float(sd/(4*om))
            f.append(de)
            tf.clear()
        elif p0>1 and z[y]>=0 and z[y]<(p0-1):
            f.append(0)
        elif (p0>0 and z[y]>=(1+p0)):
            f.append(1.0002)
        elif (p0==0 and z[y]>=0 ):
            f.append(1.0002)
        else:
            f.append(1.0002)

    flux=np.array(f)
    return flux
   

        
            
        
    
        
        
        
        
        
        
    
    
