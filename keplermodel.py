import numpy as np
def Zlist(t,per,inc,T0,a):
    z=np.empty(len(t))
    z[:]=(np.NAN)
    w=(2*np.pi)/(per)
    alpha=(90-inc)/(np.pi*180)
    time=t-T0
    phase=time/(per);phase -=np.floor(phase)
    trans=np.where(np.logical_or(phase>0.75,phase<0.25))[0]
    for i in trans:
        z[i]=a*np.sqrt((np.sin(w*(time[i])))**2+(np.cos(inc*(np.pi/180))*np.cos(w*(time[i])))**2)

    return z


