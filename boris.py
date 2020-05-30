import numpy as np

class Settings:
    t_start=0.0
    t_end=np.pi * 100.0
    dt=0.01

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

def add(x,y):
    if not isinstance(y,(Vector3,Position,Velocity)):
        y=Vector3(y,y,y)
    zx=x.x+y.x
    zy=x.y+y.y
    zz=x.z+y.z
    return Vector3(zx,zy,zz)

def mul(x,y):
    if not isinstance(y,(Vector3,Position,Velocity)):
        y=Vector3(y,y,y)
    zx=x.x*y.x
    zy=x.y*y.y
    zz=x.z*y.z
    return Vector3(zx,zy,zz)

Vector3.__add__=add
Vector3.__radd__=add
Vector3.__mul__=mul
Vector3.__rmul__=mul

class Position(Vector3):
    def update(self,v):
        new=self+v*Settings.dt
        self.substitute(new)

class Velocity(Vector3):
    def update(self,e,b):
        self.prev=Vector3(self.x,self.y,self.z)
        ta=b.mag**2*Settings.dt**2*0.25
        fac=2.0*Settings.dt*0.5/(ta+1.0)
        v1=self+e*Settings.dt*0.5
        v3=v1+v1.cross(b)*Settings.dt*0.5
        v2=v1+v3.cross(b)*fac
        new=v2+e*Settings.dt*0.5
        self.substitute(new)
    
def geteb(x):
    e=Vector3(0.0,0.0,0.0)
    b=Vector3(0.0,0.0,1.0)
    return e,b

def getvn(v):
    vnx=(v.x+v.prev.x)*0.5
    vny=(v.y+v.prev.y)*0.5
    vnz=(v.z+v.prev.z)*0.5
    return Velocity(vnx,vny,vnz)


x=Position(0.0,0.0,0.0)
v=Velocity(0.0,1.0,0.0)
e,b=geteb(x)
v.update(e,b)
file=open('boris_py.dat','w')
t=Settings.t_start

while t<=Settings.t_end:
    e,b=geteb(x)
    v.update(e,b)
    vn=(v+v.prev)*0.5
    pene=vn.mag**2*0.5
    txt='{:15.6F} '.format(t)+str(x)+' '+str(vn)+' {:15.6F} \n'.format(pene)
    file.write(txt)
    x.update(vn)
    t+=Settings.dt
file.close()
