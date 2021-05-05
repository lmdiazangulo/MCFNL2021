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
    
    _timeStepPrint = 100

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
            p["valuesh"] = [np.zeros((1,Nx[1]))]

        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIds(box)
            source["index"] = ids

        self.old = Fields(e = np.zeros( mesh.pos.size ),
                          h = np.zeros( mesh.pos.size-1 ) )

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

    def _updateE(self, t, dt):
        (e, h) = self.old.get()
        eNew = np.zeros( self.old.e.shape )
        cE = dt / sp.epsilon_0 / self._mesh.steps()
        eNew[1:-1] = e[1:-1] + cE * (h[1:] - h[:-1])
        
        # Boundary conditions
        for b in range(len(self._mesh.bounds)):
            bound = self._mesh.bounds[b]
            if b == 0:
                ind = 0
            else:
                ind = -1
            if bound == "pec":
                eNew[ind] = 0.0
            elif bound == "pmc":
                if ind == 0:
                    eNew[ind] = e[ind] + cE * 2.0 * h[ind]
                elif ind == -1:
                    eNew[ind] = e[ind] - cE * 2.0 * h[ind]
            elif bound == "mur":
                c0 = sp.speed_of_light
                dx = self._mesh.steps()
                if ind == 0:
                    eNew[0] = e[1] + \
                       (c0 * dt - dx)/(c0 * dt + dx)  * (eNew[1]-e[0]) 
                elif ind == -1:
                    eNew[ind] = e[ind-1] + \
                       (c0 * dt - dx)/(c0 * dt + dx)  * (eNew[ind-1]-e[ind])
            else:
                raise ValueError("Unrecognized boundary type")

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    eNew[source["index"]] += Solver._gaussian(t, \
                        magnitude["gaussianDelay"], \
                        magnitude["gaussianSpread"] ) * \
                            dt / self._mesh.steps() * sp.speed_of_light
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
            else:
                raise ValueError("Invalid source type: " + source["type"])

        e[:] = eNew[:]
        
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.h.shape )
        (e, h) = self.old.get()
        cH = dt / sp.mu_0 / self._mesh.steps()
        hNew[:] = h[:] + cH * (e[1:] - e[:-1])

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    hNew[source["index"]-1 ] -= Solver._gaussian((t-dt/2), \
                        magnitude["gaussianDelay"], \
                        magnitude["gaussianSpread"] ) * \
                            dt / self._mesh.steps() * sp.speed_of_light / \
                                np.sqrt(sp.mu_0 / sp.epsilon_0)
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
            else:
                raise ValueError("Invalid source type: " + source["type"])

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
                valuesh = np.zeros(ids[U]-ids[L])
                valuesh[:] = self.old.h[ ids[0]:ids[1] ]
                p["values"].append(values)
                p["valuesh"].append(valuesh)

    @staticmethod
    def _gaussian(x, delay, spread):
        return np.exp( - ((x-delay)**2 / (2*spread**2)) )
