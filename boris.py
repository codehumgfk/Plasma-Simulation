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

def dot(v1,v2):
    vx=v1.x*v2.x
    vy=v1.y*v2.y
    vz=v1.z*v2.z
    return Vector3(vx,vy,vz)

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
        #self.getv1(e,b)
        #self.getv3(e,b)
        #self.getv2(e,b)
        v1=self+e*Settings.dt*0.5
        v3=v1+v1.cross(b)*Settings.dt*0.5
        v2=v1+v3.cross(b)*fac
        #self.x=v2.x+e.x*dt*0.5
        #self.y=v2.y+e.y*dt*0.5
        #self.z=v2.z+e.z*dt*0.5
        new=v2+e*Settings.dt*0.5
        self.substitute(new)
    '''
    def getv1(self,e,b):
        v1x=self.x+e.x*dt*0.5
        v1y=self.y+e.y*dt*0.5
        v1z=self.z+e.z*dt*0.5
        self.v1=Vector3(v1x,v1y,v1z)

    def getv2(self,e,b):
        v1=self.v1
        v3=self.v3
        v2x=v1.x+(v3.y*b.z-v3.z*b.y)*self.fac
        v2y=v1.y+(v3.z*b.x-v3.x*b.z)*self.fac
        v2z=v1.z+(v3.x*b.y-v3.y*b.x)*self.fac
        self.v2=Vector3(v2x,v2y,v2z)

    def getv3(self,e,b):
        v1=self.v1
        v3x=v1.x+(v1.y*b.z-v1.z*b.y)*dt*0.5
        v3y=v1.y+(v1.z*b.x-v1.x*b.z)*dt*0.5
        v3z=v1.z+(v1.x*b.y-v1.y*b.x)*dt*0.5
        self.v3=Vector3(v3x,v3y,v3z)
    '''

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
    #bt2=b.x**2+b.y**2+b.z**2
    #ta=b.mag**2*dt**2*0.25
    #v.fac=2.0*dt*0.5/(ta+1.0)
    v.update(e,b)
    #vn=getvn(v)
    vn=(v+v.prev)*0.5
    #pene=(vn.x**2+vn.y**2+vn.z**2)*0.5
    pene=vn.mag**2*0.5
    txt='{:15.6F} '.format(t)+str(x)+' '+str(vn)+' {:15.6F} \n'.format(pene)
    #file.write(t,x.x,x.y,x.z,vn.x,vn.y,vn.z,pene,'\n')
    file.write(txt)
    x.update(vn)
    t+=Settings.dt
file.close()
