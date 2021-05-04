import math
import numpy as np
import scipy.constants as sp
import copy
import time

L = 0 # Lower
U = 1 # Upper

class Fields: 
    def __init__(self, e, h, j):
        self.e = e
        self.h = h
        self.j = j

    def get(self):
        return (self.e, self.h, self.j)

class SolverPartial:

    _timeStepPrint = 100

    def __init__(self, mesh, options, probes, sources, data):
        self.options = options
        
        self._mesh = copy.deepcopy(mesh)
       
        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIds(box)
            Nx = abs(ids)
            p["mesh"] = {"origin": box[L], "steps": abs(box[U]-box[L]) / Nx}
            p["indices"] = ids
            p["time"]   = [0.0]
            p["values"] = [np.zeros((1,Nx[1]))]

        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIds(box)
            source["index"] = ids

        self.data = data
        #self.k = 1e5
        self.k = sp.e / sp.hbar
        self.l = len(self.data['ap'])

        e = np.zeros( mesh.pos.size )
        x = np.linspace(0,1,len(e))

        # Condicion inicial 1

        e[:] = SolverPartial._gaussian(x[:], data["InitialPulse"]["Amp"], data["InitialPulse"]["Mean"], data["InitialPulse"]["Std"])

        self.old = Fields(e,
                          h = np.zeros( mesh.pos.size-1 ),
                          j = np.zeros( ( self.l , mesh.pos.size ) , dtype=complex ) )


    def solve(self, finalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        numberOfTimeSteps = int(finalTime / dt)
        for n in range(numberOfTimeSteps):
            self._updateE(t, dt)
            t += dt/2.0
            self._updateH(t, dt)
            t += dt/2.0
            self._updateProbes(t)
    
            if n % self._timeStepPrint == 0 or n+1 == numberOfTimeSteps:
                remaining = (time.time() - tic) * \
                    (numberOfTimeSteps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("    Step: %6d of %6d. Remaining: %2.0f:%02.0f"% (n, \
                    numberOfTimeSteps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))

    def _dt(self):
        return self.options["cfl"] * self._mesh.steps() / sp.speed_of_light  

    def timeStep(self):
        return self._dt()

    def getProbes(self):
        res = self._probes
        return res

    def _kp(self, coord):
        dt=self._dt()
        a=np.multiply(self.k,self.data['ap'])
        return (1 + complex(a[coord][0],a[coord][1]) * dt / 2) / (1 - complex(a[coord][0],a[coord][1]) *dt / 2)

    def _betap(self, coord):
        dt=self._dt()
        a=np.multiply(self.k,self.data['ap'])
        c=np.multiply(0.0001*self.k,self.data['cp'])
        return ( sp.epsilon_0 * complex(c[coord][0],c[coord][1]) * dt ) / ( 1 - complex(a[coord][0],a[coord][1]) * dt / 2 )

    def _arr(self):

        vkp = np.zeros ( self.l , dtype=complex )
        vbetap = np.zeros ( self.l , dtype=complex )

        for i in range(self.l):
            vkp[i] = self._kp(i)
            vbetap[i] = self._betap(i)
        
        self.vkp = vkp
        self.vbetap = vbetap 

        return (self.vkp, self.vbetap)

    def _updateE(self, t, dt):

        (e,h,j) = self.old.get()

        eNew = np.zeros( self.old.e.shape )
        jNew = np.zeros( self.old.j.shape , dtype=complex)
        SumPol = np.zeros( self.old.j.shape , dtype=complex)

        self._arr()
        for i in range(self.l):
            SumPol[i]=(1+self.vkp[i])*j[i,:]
        Suma=np.sum(SumPol,axis=0)

        Rvbetap = np.real(np.sum(self.vbetap))
        Rvkp = np.real(np.sum(self.vkp))
        
        dim = int((len(e)-1)/2)

        cE = dt / sp.epsilon_0 / self._mesh.steps()

        sigma = sp.mu_0 / sp.epsilon_0 * self.data['medio']['condc']

        c1=(2 * sp.epsilon_0 * self.data['medio']['epsinf'] + 2 * Rvbetap - sigma * dt) / \
            (2 * sp.epsilon_0 * self.data['medio']['epsinf'] + 2 * Rvbetap + sigma * dt)  


        # Se cambian los límites que se deseen para indicar el medio dispersivo
        for i in range(1,len(e)-1):
            if  0 < i < 0.6*dim:

                eNew[i] = e[i] + cE * (h[i] - h[i-1])

            else :
                eNew[i] = c1 * e[i] + \
                2 * dt * ( h[i] - h[i-1] - np.real( Suma[i] ) ) / \
                        (2 * sp.epsilon_0 * self.data['medio']['epsinf'] + 2 * Rvbetap + sigma * dt)

                for W in range(self.l):
                    jNew[W,i] = self._kp(W) * j[W,i] + self._betap(W) * ( eNew[i] - e[i] ) / dt
        j=jNew 

        # Boundary conditions
        for b in range(len(self._mesh.bounds)):
            bound = self._mesh.bounds[b]
            if b == 0:
                ind = 0
            else:
                ind = -1
            if bound == "pec":
                eNew[ind] = 0.0

        e[:] = eNew[:]
    
        return (c1,Rvbetap,Rvkp,Suma)
        
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.h.shape )
        (e, h, j) = self.old.get()
        cH = dt / sp.mu_0 / self._mesh.steps()
        hNew[:] = h[:] + cH * (e[1:] - e[:-1])
        h[:] = hNew[:]

    def _updateProbes(self, t):
        for p in self._probes:
            if "samplingPeriod" not in p or \
               "samplingPeriod" in p and \
               (t/p["samplingPeriod"] >= len(p["time"])):
                p["time"].append(t)
                ids = p["indices"]
                values = np.zeros(ids[U]-ids[L])
                values[:] = self.old.e[ ids[0]:ids[1] ]
                p["values"].append(values)

    @staticmethod
    def _gaussian(x, amp, delay, spread):
        return amp*np.exp( - ((x-delay)**2 / (2*spread**2)) )