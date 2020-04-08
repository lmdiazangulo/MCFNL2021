import math
import numpy as np
import scipy.constants as sp
import copy
import time

L = 0 # Lower
U = 1 # Upper

class Fields: 
    def __init__(self, e, h):
        self.e = e
        self.h = h
    
    def get(self):
        return (self.e, self.h)

class Solver:
    
    __timeStepPrint = 100

    def __init__(self, mesh, options, probes, sources):
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

        self.old = Fields(e = np.zeros( mesh.pos.size ),
                          h = np.zeros( mesh.pos.size-1 ) )

    def solve(self, dimensionalFinalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        numberOfTimeSteps = \
            int(dimensionalFinalTime * sp.speed_of_light / dt)
        for n in range(numberOfTimeSteps):
        
            self._updateE(t, dt)
            t += dt/2.0

            self._updateH(t, dt)
            t += dt/2.0

            self._updateProbes(t)
    
            if n % self.__timeStepPrint == 0 or n+1 == numberOfTimeSteps:
                remaining = (time.time() - tic) * \
                    (numberOfTimeSteps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("    Step: %6d of %6d. Remaining: %2.0f:%02.0f"% (n, \
                    numberOfTimeSteps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))

    def _dt(self):
        return self.options["cfl"] * self._mesh.steps() / math.sqrt(2.0)  

    def timeStep(self):
        return self._dt() / sp.speed_of_light

    def getProbes(self):
        res = self._probes
        return res

    def _updateE(self, t, dt):
        (e, h) = self.old.get()
        eNew = np.zeros( self.old.e.shape )
        eNew[1:-1] = e[1:-1] + dt * (h[1:] - h[:-1])
        
        # Boundary conditions
        for bound in self._mesh.bounds:
            if bound == "pec":
                eNew[ 0] = 0.0
                eNew[-1] = 0.0
            else:
                raise ValueError("Unrecognized boundary type")

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    delay  = sp.speed_of_light * magnitude["gaussianDelay"]
                    spread = sp.speed_of_light * magnitude["gaussianSpread"]
                    id = source["index"]
                    eNew[id] += Solver._gaussian(t, delay, spread) * dt
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
            else:
                raise ValueError("Invalid source type: " + source["type"])

        e[:] = eNew[:]
        
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.h.shape )
        (e, h) = self.old.get()
        hNew[:] = h[:] + dt * (e[1:] - e[:-1])
        h[:] = hNew[:]
            
    def _updateProbes(self, t):
        for p in self._probes:
            dimensionalTime = t/sp.speed_of_light
            writeStep = "samplingPeriod" in p \
                and (dimensionalTime/p["samplingPeriod"] >= len(p["time"]))
            writeStep = writeStep or "samplingPeriod" not in p
            if writeStep:
                p["time"].append(dimensionalTime)
                ids = p["indices"]
                values = np.zeros(ids[U]-ids[L])
                values[:] = self.old.e[ ids[0]:ids[1] ]
                p["values"].append(values)

    @staticmethod
    def _gaussian(x, delay, spread):
        return np.exp( - ((x-delay)**2 / (2*spread**2)) )
