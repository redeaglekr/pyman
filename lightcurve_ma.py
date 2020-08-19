##This is a program written for plotting analytical(Theoretical) Light curve using Mandel and Agol
##This can accomodate Linear limb darkening
###This uses all 11 conditions of table 1 in the paper

from hypergeometric import Ns,Ms
from keplermodel import Zlist
import math as ma
import mpmath as mp
import numpy as np
def mlal2(time,rat,inc,per,T0,p0):
    #w=2*np.pi/(per)
    #time=time-T0
    #phase=time/(per);phase -=np.floor(phase)
    #trans=np.where(np.logical_or(phase>0.75,phase<0.25))[0]
    #z=rat*np.sqrt((np.sin(w*time))**2+(np.cos(inc*(np.pi/180))*np.cos(w*time))**2)
    z=Zlist(time,per,inc,T0,rat)
    c=[0.5807,-0.0336,0.4824,-0.2851]
    #c=[c1,c2,c3,c4]
    #c=[0,0,0,0]
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
            tk=[]
            n0=Ns(z[y],p0,0)*k0
            tk.append(n0)
            for q in range(1,5,1):
                s3=Ns(z[y],p0,q)*(c[q-1]/(q+4))
                tk.append(s3)
            summ=np.sum(tk)
            fg=(1/(2*np.pi*om))*(summ)
            hg=float(1-fg)
            f.append(hg)
            tk.clear()
           # print("case 2,z value=",z[y])
        elif  ((p0>0 and p0<0.5) and (z[y]>p0 and z[y]<(1-p0))):
            lm=[]
            for r in range(1,4,1):
                s4=Ms(z[y],p0,r)*(c[r-1]/(r+4))
                lm.append(s4)
            summ2=np.sum(lm)
            #print(summ2)
            l=(p0**2) *( 1-((p0**2)/2)-z[y]**2 )
            fn=(c0*((p0)**2))+(2*summ2)+(c4*l)
            fn2=fn/(4*om)
            fnm=float(1-fn2)
            f.append(fnm)
            lm.clear()
           # print("case 3,z value=",z[y])

        elif  ((p0>0 and p0<0.5) and (z[y]==1-p0)):
            lmc=[]
            for b in range(1,4,1):
                s5=Ms(z[y],p0,b)*(c[b-1]/(b+4))
                lmc.append(s5)
            summ3=np.sum(lmc)
            l5=((p0)**2)*((1-(p0)**2)/(2-(z[y]**2)))
            fn=(c0*((p0)**2))+(2*summ3)+(c4*l5)
            fn2=fn/(4*om)
            fni=float(1-fn2)
            f.append(fni)
            lmc.clear()
            #print("case 4,z value=",z[y])


        elif (p0>0 and p0<0.5) and (z[y]==p0):
            m5=[]
            so=k0*hyp2f1(0.5,-1,1,4*(p0)**2)
            m5.append(so)
            for v in range(1,5,1):
                s6=(c[v-1]/(v+4))*(hyyp2f1(0.5,-(v+4)/4,1,4*(p0)**2))
                m5.append(s6)
            c5=np.sum(m5)
            t2=c5/(2*om)
            fv=float(0.5+t2)
            f.append(fv)
            m5.clear()
          #  print("case 5,z value=",z[y])
           
        elif (p0 == 0.5 )and (z[y] == 0.5):
            nm3=[]
            mo=k0*(gamma(1.5))/(gamma(2))
            nm3.append(mo)
            for h in range(1,5,1):
                g=((c[h-1])/(h+4))*(gamma(1.5+(h/4))/gamma(2+(h/4)))
                nm3.append(g)
            df=np.sum(nm3)   
            sm=df/(2*np.sqrt(np.pi)*om)
            f6=float(0.5+sm)
            f.append(f6)
            nm3.clear()
           # print("case 6,z value=",z[y])
            
        elif (p0>0.5) and (z[y]==p0):
            gf=[]
            bo=ko*(beta(0.5,2)*hyp2f1(0.5,0.5,2.5,(1/(4*p0**2))))
            gf.append(bo)
            for x in range(1,5,1):
                qw=(c[x-1]/(x+4))*(beta(0.5,(x+8)/4)*hyp2f1(0.5,0.5,2.5+(x/4),1/(4*p0**2)))
                gf.append(qw)
            dfg=np.sum(gf)
            fi=dfg/(4*p0*np.pi*om)
            flo=float(0.5+fi)
            f.append(flo)
            gf.clear()
            #print("case 7,z value=",z[y])
            
        elif p0>0.5 and (z[y]>=abs(1-p0) and z[y]<p0):#case 8
            jhi=[]
            df3=Ns(z[y],p0,0)*k0
            jhi.append(df3)
            for u in range(1,5,1):
                s7=((c[u-1])*Ns(z[y],p0,u))/(u+4)
                jhi.append(s7)
            a=np.sum(jhi)
            ds=float(-(a)/(om))
            f.append(ds)
            jhi.clear()
            #print("case 8,z value=",z[y])
           
        elif (p0>0 and p0<1) and (z[y]>0 and z[y]<(0.5-abs(p0-0.5))):
            ed=[]
            for n in range(1,4,1):
                mf=(c[n-1]/(n+4))*(Ms(z[y],p0,n))
                ed.append(mf)
            bn=np.sum(ed)
            lm=((p0)**2)*((1-(p0)**2)/(2-z[y]**2))
            op=(c0*(1-(p0)**2))+(c4*(0.5-lm))-(2*bn)
            fml=float(op/(4*om))
            f.append(fml)
            ed.clear()
            #print("case 9,z value=",z[y])

        elif (p0>0 and p0<1) and (z[y]==0):
            ttf=[]
            t0=k0*(1-(p0)**2)
            ttf.append(t0)
            for s in range(1,5,1):
                ti=(c[s-1])*((1-(p0)**2)**((s+4)/4))/(s+4)
                ttf.append(ti)
            sd=np.sum(ttf)
            de=float(sd/(om))
            f.append(de)
            ttf.clear()
            #print("case 10,z value=",z[y])

        elif (p0>1) and  (z[y]>=0 and z[y]<(p0-1)):
            f.append(0)

        elif (p0>0 and z[y]>=(1+p0)):
            f.append(1.0002)

        elif (p0==0 and z[y]>=0 ):
            f.append(1.0002)


        else:
            f.append(1.0002)
    
    
    flux=np.array(f)     
    return flux



            
   

        
            
        
    
        
        
        
        
        
        
    
    
