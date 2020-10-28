import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

startWithRotation = True
initialSpin = 1.8

SphereConst = 4/3*np.pi
GRAVITATIONAL_CONSTANT = 0.8

class Particle():
    def __init__(self,pos,v,m,rho):
        self.corr_m = 1
        self.corr_rad = 1
        self.CollTh = .001
        self.position = pos
        self.vel = v #np.array(3)
        self.mass = m*self.corr_m
        self.density = rho
        self.radius = self.corr_m*SphereConst*rho*m**(1/3) # Corr
        self.angular_mom = (2/5)*self.mass*self.radius*self.vel
        self.SelfGpotentialE = -(3/5)*GRAVITATIONAL_CONSTANT*self.mass**2/self.radius
        self.KineticE = 0.5*self.mass*np.sum((self.vel**2))
        self.temp = self.KineticE/(1.5*1.380649*10e-23)
        self.combineWith=False
        self.acc = []

    def unit_vec(array):
        array = array/np.sqrt(np.sum(array[:]**2))

    def __del__(self):
        #print("Body removed")
        None

    def force(self,other):
        d = np.sqrt(np.sum((self.position - other.position)**2))
        force = GRAVITATIONAL_CONSTANT*self.mass*other.mass/d**2
        self.acc.append(force / self.mass)
        
    def move(self, time):
        self.acc = np.sum(np.array(self.acc))
        self.position += self.vel + self.acc*time
        self.acc = []
    
    def attract(self, other, change_CT = False):
        d = np.sqrt(np.sum((self.position - other.position)**2))
        if d > 1e25:
            return 2

        #if change_CT:                      #CHANGE CT
        #    self.CollTh = self.CollTh + 1                      #CHANGE CT
        elif d<=((self.radius + other.radius)*self.CollTh):
            self.merge(other)
            #print("Mergre success")
            return 1
        
        else:
            force = GRAVITATIONAL_CONSTANT*self.mass*other.mass/d**2
            accel1 = force / self.mass
            accel2 = force / other.mass

            ds_hat = (self.position - other.position)/d
            self.vel -= accel1 * ds_hat
            other.vel -= accel2 * ds_hat
            #print("Attraction success")
            return 0

    def merge(self, other):
        #print("Merging bodies....")
        self.mass = (self.mass + other.mass)*self.corr_m
        self.position = (self.mass*self.position+other.mass*other.position)/(self.mass+other.mass)
        self.vel = self.mass*self.vel + other.mass*other.vel/(self.mass+other.mass)
        self.angular_mom = (self.angular_mom + other.angular_mom)
        self.radius = (self.radius + other.radius)*self.corr_rad#*(1e-6*self.angular_mom)
        self.density = self.mass/(SphereConst*self.radius**(1/3))

class space():
    def __init__(self, n, iter):
        self.n = n
        time = np.linspace(0,100000,10)
        v, pos = np.zeros((self.n,3))*1e1,np.random.random((self.n,3))*1.5*10e2,
        m, rho = np.random.random(self.n)*1e-2, np.random.random(self.n)*1e6
        self.particle=[]
        self.temp = 0
        for i in range((self.n)):
            self.particle.append(Particle(pos[i],v[i],m[i],rho[i]))
            if startWithRotation:
                self.particle[i].vel[1] += initialSpin*np.random.random(1)
                self.particle[i].vel[2] += initialSpin*np.random.random(1)
        self.varm, self.varrad, self.vardensity = [],[],[]
        self.varSelfGpotentialE, self.varKineticE, self.vartemp = [],[],[]

        print("Particles initialized: ",str(self.n)," particles",str(iter)," iterations")
        self.start(time)

    def start(self,time):
        CT = False
        for l in time:
            #print("\nIter: ",l+1)
            self.update(l,CT)
            for mov in range(self.n):
                for acc in range(self.n):
                    if acc!=mov:
                        self.particle[mov].force(self.particle[acc])
                self.particle[mov].move(l)
#   Make CT a function of energy and not loop                


    def update(self,l,CT):
        i,j = 0,0
        while i<self.n :
            while j<self.n :
                if j!=i :
                    if self.particle[i].attract(self.particle[j],CT)==1 :
                        for k in range(j+1, self.n):
                            self.particle[k-1] = self.particle[k]
                            self.n -= 1
                        del self.particle[-1]
                        self.n = len(self.particle)
                    
                    elif (self.particle[i].attract(self.particle[j],CT)==2) :
                        for k in range(j+1, self.n):
                            self.particle[k-1] = self.particle[k]
                            self.n -= 1
                        del self.particle[-1]
                        self.n = len(self.particle)
                j+=1
            i+=1
        print("Particles left: ", self.n)
        #self.plot()
        self.variation()

    def variation(self):
        varm,varrad,vardensity,varSelfGpotentialE,varKineticE,vartemp =0,0,0,0,0,0
        for i in range(self.n):
            varm += self.particle[i].mass
            varrad += self.particle[i].radius
            vardensity += self.particle[i].density
            varSelfGpotentialE += self.particle[i].SelfGpotentialE
            varKineticE += self.particle[i].KineticE
            vartemp += self.particle[i].temp
        self.vartemp.append(vartemp/self.n)
        self.varm.append(varm/self.n)
        self.varrad.append(varrad/self.n)
        self.vardensity.append(vardensity/self.n)
        self.varSelfGpotentialE.append(varSelfGpotentialE/self.n)
        self.varKineticE.append(varKineticE/self.n)

    def plot(self):
        m,rad,density,angular_mom,vel= [],[],[],[],[]
        #x1, y1, z1, x2, y2, z2 = [],[],[],[],[],[]
        for i in range(self.n):
            m.append(np.log10(self.particle[i].mass))
            #HEEEEEEEEREEE

       #plt.title('Mass')
       #plt.hist(np.array(m))
       #plt.show()
        plt.figure(figsize=(18,3))
        plt.semilogy(self.varm,label='mass(gm)')
        plt.semilogy(self.varrad,label='radius(m)')
        plt.semilogy(self.vardensity,label='density')
        #plt.semilogy(self.varSelfGpotentialE,label='Grav Pot')
        #plt.semilogy(self.varKineticE,label='KE')
        #plt.semilogy(self.vartemp,label='temp')
        plt.legend(loc=1)
        plt.show()
        #THERE
