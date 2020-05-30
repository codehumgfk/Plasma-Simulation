import numpy as np

class Settings:
    t_start=0.0
    t_end=np.pi * 12837.0
    dt=2.0 * np.pi * 0.01
    q0=1.60217662e-19 # elementary chage in [C]
    m0=1.6726219e-27 #proton mass in [kg]
    re=6.3712e6 #Earth radius in [m]
    g0=-2.980592e4 # Earth's dipole g0 coeff. in [nT]

dt=Settings.dt
q0=Settings.q0
m0=Settings.m0
re=Settings.re
g0=Settings.g0

deg2rad=np.arctan(1.0)/45.0
t=Settings.t_start
t_end=Settings.t_end

pos_r=6.6 # initial particle position (R) in Re
pos_lat=0.0 # initial particle position (latitude) in degree
pos_long=0.0 # initial particle position (longitude) in degree
alpha=45.0 #initial particle pitch-angle in degree
bp=-g0*np.sqrt(1.00+3.00*np.sin(pos_lat*deg2rad)**2)/pos_r**3
bm=bp/np.sin(alpha*deg2rad)**2

qm=1.00 # charge / mass normalized by proton
b0=bm # normalized parameter for B in [nT]
v0=3.0e6 # normalized parameter for velocity in [m/s]
omega0=abs(q0*b0*1.0e-9/m0) # gyrofrequency at b0 [1/s]
t0=1.0e0/omega0 # normalized parameter for t in [sec]
r0=v0/omega0 # normalized parameter for x in [m]
e0=v0*b0 # normalized parameter for E in [V/m]

class Vector3:
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z

    def substitute(self,new):
        self.x=new.x
        self.y=new.y
        self.z=new.z

    def cross(self,b):
        return cross(self,b)

    def dot(self,b):
        return dot(self,b)

    @property
    def mag(self):
        return np.sqrt(self.x**2+self.y**2+self.z**2)

    def update():
        raise NotImplementedError()

    def __str__(self):
        return '{:15.6F} {:15.6F} {:15.6F}'.format(self.x,self.y,self.z)


def cross(v3,b):
    cx=v3.y*b.z-v3.z*b.y
    cy=v3.z*b.x-v3.x*b.z
    cz=v3.x*b.y-v3.y*b.x
    return Vector3(cx,cy,cz)

def dot(v3,b):
    return v3.x*b.x+v3.y*b.y+v3.z*b.z

def add(x,y):
    if not isinstance(y,(Vector3,Position,Velocity)):
        y=Vector3(y,y,y)
    zx=x.x+y.x
    zy=x.y+y.y
    zz=x.z+y.z
    return Vector3(zx,zy,zz)

def sub(x,y):
    if not isinstance(y,(Vector3,Position,Velocity)):
        y=Vector3(y,y,y)
    zx=x.x-y.x
    zy=x.y-y.y
    zz=x.z-y.z
    return Vector3(zx,zy,zz)

def rsub(x,y):
    return sub(y,x)

def mul(x,y):
    if not isinstance(y,(Vector3,Position,Velocity)):
        y=Vector3(y,y,y)
    zx=x.x*y.x
    zy=x.y*y.y
    zz=x.z*y.z
    return Vector3(zx,zy,zz)

def div(x,y):
    if not isinstance(y,(Vector3,Position,Velocity)):
        zx=x.x/y
        zy=x.y/y
        zz=x.z/y
        return Vector3(zx,zy,zz)

Vector3.__add__=add
Vector3.__radd__=add
Vector3.__sub__=sub
Vector3.__rsub__=rsub
Vector3.__mul__=mul
Vector3.__rmul__=mul
Vector3.__truediv__=div

class Position(Vector3):
    def update(self,v):
        new=self+v*dt
        self.substitute(new)

class Velocity(Vector3):
    def update(self,e,b):
        self.prev=Vector3(self.x,self.y,self.z)
        ta=b.mag**2*dt**2*0.25
        fac=2.0*dt*0.5/(ta+1.0)
        v1=self+e*qm*dt*0.5
        v3=v1+v1.cross(b)*qm*dt*0.5
        v2=v1+v3.cross(b)*qm*fac
        new=v2+e*qm*dt*0.5
        self.substitute(new)

class Er(Vector3):
    def __init__(self,theta,phi):
        self.x=np.sin(theta)*np.cos(phi)
        self.y=np.sin(theta)*np.sin(phi)
        self.z=np.cos(theta)

class Etheta(Vector3):
    def __init__(self,theta,phi):
        self.x=np.cos(theta)*np.cos(phi)
        self.y=np.cos(theta)*np.sin(phi)
        self.z=-np.sin(theta)

class Ephi(Vector3):
    def __init__(self,phi):
        self.x=-np.sin(phi)
        self.y=np.cos(phi)
        self.z=0

def getdb(x):
    xx=x*r0/re
    rxy=np.sqrt(xx.x**2+xx.y**2)
    if xx.mag > 1.0:
        theta=np.arccos(xx.z/xx.mag)
        if rxy > 0.0:
            phi=np.arctan2(xx.y/rxy,xx.x/rxy)
        else:
            phi=0.0
    else:
        theta=np.arcsin(rxy)
        b=Vector3(0.0,0.0,np.sqrt(np.cos(theta)**2*3.0+1.0)*g0/b0)
    er=Er(theta,phi)
    etheta=Etheta(theta,phi)
    ephi=Ephi(phi)
    b=(2.0*np.cos(theta)*er+np.sin(theta)*etheta)*g0/xx.mag**3/b0
    return b

theta=pos_lat*deg2rad
phi=pos_long*deg2rad
er=Er(theta,phi)

ca=np.cos(alpha*deg2rad)
sa=np.sin(alpha*deg2rad)

x=pos_r*er*re/r0
x=Position(x.x,x.y,x.z)
b=getdb(x)
eb=b/b.mag
ephi=Ephi(phi)
e=Vector3(0.0,0.0,0.0)
v=eb*ca+ephi*sa
v=v-(e+v.cross(b))*qm*dt*0.50
v=Velocity(v.x,v.y,v.z)

file=open('dipole_py.dat','w')

while t <= t_end:
    b=getdb(x)
    v.update(e,b)
    vn=(v+v.prev)*0.5
    pene=vn.mag**2*0.5
    #txt='{:15.6F} '.format(t)+str(x)+' '+str(vn)+' {:15.6F} \n'.format(pene) #normalized
    txt='{:15.6F} '.format(t*t0)+str(x*r0/re)+' '+str(vn*v0*1.0e-3)+' {:15.6F} \n'.format(pene)
    file.write(txt)
    x.update(vn)
    t+=dt
file.close()
