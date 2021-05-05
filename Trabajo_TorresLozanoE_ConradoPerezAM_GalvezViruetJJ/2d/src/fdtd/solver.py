import math
import numpy as np
import scipy.constants as sp
import copy
import time

from fdtd.common import X, Y, L, U
from fdtd.sources import TMn_Ex, TMn_Ey, TMn_Hz, gaussian, step, fc_f 

def subsId(id):
    if id is None:
        return -1
    else:
        return id-1

class Solver:
    
    class Fields: 
        def __init__(self, ex, ey, hz):
            self.ex = ex
            self.ey = ey
            self.hz = hz
        
        def get(self):
            return (self.ex, self.ey, self.hz)

    __timeStepPrint = 5000

    def __init__(self, mesh, options, probes, sources, materials):
        self.options = options
        
        self._mesh = copy.deepcopy(mesh)
        (dX, dY) = self._mesh.steps()
        self.posEx = [np.delete(mesh.pos[X] - dX/2, 0), mesh.pos[Y]]
        self.posEy = [mesh.pos[X], np.delete(mesh.pos[Y] - dX/2, 0)]
        self.posHz = [np.delete(mesh.pos[X] - dX/2, 0), np.delete(mesh.pos[Y] - dX/2, 0)] 

        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIdx(box)
            Nxy = abs(ids[Y] - ids[X])
            p["mesh"] = {"origin": box[L], "steps": (dX, dY), "posEx": self.posEx, \
                            "posEy": self.posEy, "posHz": self.posHz}
            p["indices"] = ids
            p["time"]   = [0.0]
            p["values"] = [np.zeros((Nxy[X], Nxy[Y]))]
            p["valuese_mod"] = [np.zeros((Nxy[X], Nxy[Y]))]
            p["valuese_x"] = [np.zeros((Nxy[X], Nxy[Y]))]
            p["valuese_y"] = [np.zeros((Nxy[X], Nxy[Y]))]

        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIdx(box)
            source["index"] = ids
        
        self.old = self.Fields(
            ex = np.zeros( (mesh.pos[X].size-1, mesh.pos[Y].size  ) ),
            ey = np.zeros( (mesh.pos[X].size,   mesh.pos[Y].size-1) ),
            hz = np.zeros( (mesh.pos[X].size-1, mesh.pos[Y].size-1) ) )

        self._materials = copy.deepcopy(materials)
        self.mu = np.zeros((mesh.pos[X].size, mesh.pos[Y].size))
        self.epsilon = np.zeros((mesh.pos[X].size, mesh.pos[Y].size))
        for material in self._materials:
            box = self._mesh.elemIdToBox(material["elemId"])
            id = mesh.toIdx(box)
            self.mu[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] = material["mu_rel"]
            self.epsilon[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] = material["epsilon_rel"]

    def _dt(self):
        '''
        Gives the time step for the fdtd algorithm taking into account the 
        cfl condition for stability.
        '''
        return self.options["cfl"] * min(self._mesh.steps()) / math.sqrt(2.0)  

    def timeStep(self):
        '''
        Gives the time step for the fdtd algorithm In SI units.
        '''
        return self._dt() / sp.speed_of_light

    def getProbes(self):
        res = self._probes
        return res

    # ======================= UPDATE E =============================
    def _updateE(self, t, dt, overFields = None):
        eNew = (np.zeros( self.old.ex.shape ),
                np.zeros( self.old.ey.shape ) )
        (ex, ey, h) = self.old.get()
        e = (ex, ey)

        (dX, dY) = self._mesh.steps()
        A = dX * dY * np.array([self.epsilon[:-1,1:-1], self.epsilon[1:-1,:-1]], dtype=object) 
        eNew[X][:,1:-1] = e[X][:,1:-1] + dt/A[X]*dX * (h[:,1:] - h[:,:-1])
        eNew[Y][1:-1,:] = e[Y][1:-1,:] - dt/A[Y]*dY * (h[1:,:] - h[:-1,:])

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                (initEx, initEy, initHz) = source["fields"]
                if magnitude["type"] == "gaussian":
                    break

                elif magnitude["type"] == "TMgauss":
                    id = source["index"]

                    epsilon = self.epsilon[id[L][X], id[L][Y]]
                    mu = self.mu[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]

                    c0 = sp.speed_of_light # To change from SI units to algortihm units

                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    
                    n = source["mode"]
                    freq = source["frequency"] * 2 * np.pi / c0
                    intens = magnitude["sinIntensity"]
                    lon_y = (id[U][Y] - id[L][Y]) * dY
                    
                    fc_f(n, lon_y, self.mu[:-1,:-1], self.epsilon[:-1,:-1], source["frequency"])

                    if initEx == 1:
                        (xEx, yEx) = self._mesh.IdxToPos(id,self.posEx)
                        NxEx = len(xEx)
                        NyEx = len(yEx)
                        xEx = np.tile(xEx[...,None], (1, NyEx))
                        yEx = np.tile(yEx, (NxEx, 1))
                        tEx = t + dt/2
                        eNew[X][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                          TMn_Ex(tEx, xEx, yEx, n, intens, freq, mu, epsilon, lon_y) \
                              * gaussian(tEx, delay, spread)

                    if initEy == 1:
                        (xEy, yEy) = self._mesh.IdxToPos(id,self.posEy)
                        NxEy = len(xEy)
                        NyEy = len(yEy)
                        xEy = np.tile(xEy[...,None], (1, NyEy))
                        yEy = np.tile(yEy, (NxEy, 1))
                        tEy = t + dt/2
                        eNew[Y][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                         TMn_Ey(t, xEy, yEy, n, intens, freq, mu, epsilon, lon_y) \
                         * gaussian(tEy, delay, spread) 

                elif magnitude["type"] == "TMstep":
                    id = source["index"]

                    epsilon = self.epsilon[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]
                    mu = self.mu[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]

                    c0 = sp.speed_of_light # To change from SI units to algortihm units

                    n = source["mode"]
                    freq = source["frequency"] * 2 * np.pi / c0
                    intens = magnitude["sinIntensity"]
                    tlim = magnitude["stepTimeLimit"]
                    lon_y = (id[U][Y] - id[L][Y]) * dY

                    fc_f(n, lon_y, self.mu[:-1,:-1], self.epsilon[:-1,:-1], source["frequency"])

                    if initEx == 1:
                        (xEx, yEx) = self._mesh.IdxToPos(id,self.posEx)
                        NxEx = len(xEx)
                        NyEx = len(yEx)
                        xEx = np.tile(xEx[...,None], (1, NyEx))
                        yEx = np.tile(yEx, (NxEx, 1))
                        tEx = t + dt/2
                        eNew[X][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                         TMn_Ex(tEx, xEx, yEx, n, intens, freq, mu, epsilon, lon_y) \
                         * step(tEx, tlim * dt)

                    if initEy == 1:
                        (xEy, yEy) = self._mesh.IdxToPos(id,self.posEy)
                        NxEy = len(xEy)
                        NyEy = len(yEy)
                        xEy = np.tile(xEy[...,None], (1, NyEy))
                        yEy = np.tile(yEy, (NxEy, 1))
                        tEy = t + dt/2
                        eNew[Y][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                         TMn_Ey(tEy, xEy, yEy, n, intens, freq, mu, epsilon, lon_y) \
                         * step(tEy, tlim * dt)
                     
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
            else:
                raise ValueError("Invalid source type: " + source["type"])

        # Boundary conditions
        for bound in self._mesh.bounds:
            xy = bound.orientation()
            (lx, ux) = (bound.arrayIdx(L,X), \
                        bound.arrayIdx(U,X))
            (ly, uy) = (bound.arrayIdx(L,Y), \
                        bound.arrayIdx(U,Y))

            # In the json file: [ left, right]
            #                   [ down,    up]

            if isinstance(bound, self._mesh.BoundPEC):
                 eNew[xy][lx:ux,ly:uy] = 0.0
            
            elif isinstance(bound, self._mesh.BoundMUR):
                c0 = 1 # Light speed in algorithm units
                if xy == Y: # Left and right: we change Ey   
                    if lx == 0: # Left
                        eNew[Y][0,:] = e[Y][1,:] + \
                         (c0 * dt - dY)/(c0 * dt + dY)  * (eNew[Y][1,:]-e[Y][0,:])
                    elif lx == -1: # Right
                        eNew[Y][-1,:] = e[Y][-2,:] + \
                         (c0 * dt - dY)/(c0 * dt + dY)  * (eNew[Y][-2,:]-e[Y][-1,:])
                else: # Up and down: we change Ex   
                    if lx == 0: # Down
                        eNew[X][:,0] = e[X][:,1] + \
                         (c0 * dt - dX)/(c0 * dt + dX)  * (eNew[X][:,1]-e[X][:,0])
                    elif lx == -1: # Up
                        eNew[X][:,-1] = e[X][:,-2] + \
                         (c0 * dt - dX)/(c0 * dt + dX)  * (eNew[X][:,-2]-e[X][:,-1])  
            
            else:
                raise ValueError("Unrecognized boundary type")


        # Subgridding and updating
        e[X][:] = eNew[X][:]
        e[Y][:] = eNew[Y][:]  

    # ======================= UPDATE H =============================
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.hz.shape )
        (ex, ey, h) = self.old.get()
        
        (dX, dY) = self._mesh.steps()
        A = dX * dY * self.mu[:-1,:-1]
              
        hNew[:,:] = h[:,:] \
                     - dt/A * dY * ey[1:,  :] \
                     + dt/A * dX * ex[ :, 1:] \
                     + dt/A * dY * ey[:-1,   :] \
                     - dt/A * dX * ex[  :, :-1]
        
        # Source terms
        for source in self._sources:
            (initEx, initEy, initHz) = source["fields"]
            if initHz == 1:
                if source["type"] == "dipole":
                    magnitude = source["magnitude"]

                    if magnitude["type"] == "gaussian":
                        c0 = sp.speed_of_light # To change from SI units to algorithm units
                        delay  = c0 * magnitude["gaussianDelay"]
                        spread = c0 * magnitude["gaussianSpread"]
                        id = source["index"]
                        hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                         gaussian(t, delay, spread)*dt

                    elif magnitude["type"] == "TMgauss":
                        id = source["index"]

                        epsilon = self.epsilon[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]
                        mu = self.mu[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]

                        c0 = sp.speed_of_light # To change from SI units to algorithm units

                        delay  = c0 * magnitude["gaussianDelay"]
                        spread = c0 * magnitude["gaussianSpread"]
                    
                        n = source["mode"]
                        freq = source["frequency"] * 2 * np.pi / c0
                        intens = magnitude["sinIntensity"]
                        lon_y = (id[U][Y] - id[L][Y]) * dY

                        fc_f(n, lon_y, self.mu[:-1,:-1], self.epsilon[:-1,:-1], source["frequency"])

                        (xHz, yHz) = self._mesh.IdxToPos(id,self.posHz)
                        NxHz = len(xHz)
                        NyHz = len(yHz)
                        xHz = np.tile(xHz[...,None], (1, NyHz))
                        yHz = np.tile(yHz, (NxHz, 1))
                        hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                         TMn_Hz(t, xHz, yHz, n, intens, freq, mu, epsilon, lon_y) \
                         * gaussian(t, delay, spread)

                    elif magnitude["type"] == "TMstep":
                        id = source["index"]

                        epsilon = self.epsilon[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]
                        mu = self.mu[id[L][X]:id[U][X], id[L][Y]:id[U][Y]]

                        c0 = sp.speed_of_light # To change from SI units to algorithm units

                        n = source["mode"]
                        freq = source["frequency"] * 2 * np.pi / c0
                        intens = magnitude["sinIntensity"]
                        tlim = magnitude["stepTimeLimit"]
                        lon_y = (id[U][Y] - id[L][Y]) * dY

                        fc_f(n, lon_y, self.mu[:-1,:-1], self.epsilon[:-1,:-1], source["frequency"])

                        (xHz, yHz) = self._mesh.IdxToPos(id,self.posHz)
                        NxHz = len(xHz)
                        NyHz = len(yHz)
                        xHz = np.tile(xHz[...,None], (1, NyHz))
                        yHz = np.tile(yHz, (NxHz, 1))
                        hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                         TMn_Hz(t, xHz, yHz, n, intens, freq, mu, epsilon, lon_y) \
                         * step(t, tlim * dt)
                     
                    else:
                        raise ValueError(\
                        "Invalid source magnitude type: " + magnitude["type"])
                else:
                    raise ValueError("Invalid source type: " + source["type"])
        
        # Updating
        h[:] = hNew[:]

    # Results in SI units
    def _updateProbes(self, t):
        for p in self._probes:
            dimensionalTime = t/sp.speed_of_light
            writeStep = "samplingPeriod" in p \
                and (dimensionalTime/p["samplingPeriod"] >= len(p["time"]))
            writeStep = writeStep or "samplingPeriod" not in p
            if writeStep:
                p["time"].append(dimensionalTime)
                idx = p["indices"]

                # Calculating magnetic fields values
                values = np.zeros(tuple(idx[U]-idx[L]))
                values[:,:] = \
                    self.old.hz[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ] \
                    / np.sqrt(sp.mu_0/sp.epsilon_0)

                # Calculating electric field values at magnetic field positions
                # Ex values without first raw    
                valuesexup = \
                    np.array(self.old.ex[ idx[L][X]:idx[U][X], (idx[L][Y]+1):(idx[U][Y]+1) ])
                # Ex values without last raw
                valuesexdown = \
                    np.array(self.old.ex[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ])
                # Mean values, Ex same position as Hz
                valuesex = (valuesexup + valuesexdown)/2               
                # Ey values without first column
                valueseyright = \
                    self.old.ey[ (idx[L][X]+1):(idx[U][X]+1), idx[L][Y]:idx[U][Y] ]
                # Ey values without last comlumn
                valueseyleft = \
                    self.old.ey[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ]
                # Mean values, Ey same position as Hz
                valuesey = (valueseyright + valueseyleft)/2

                # Adding values to the results
                p["values"].append(values)
                p["valuese_mod"].append(np.array([list(map(lambda x,y: np.sqrt(x**2 +y**2), valuesex[i], valuesey[i])) for i in range(0,len(valuesex))]))
                p["valuese_x"].append(valuesex)
                p["valuese_y"].append(valuesey)

    # Solving time loop
    def solve(self, dimensionalFinalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        numberOfTimeSteps = \
            int(dimensionalFinalTime * sp.speed_of_light / dt)
        for n in range(numberOfTimeSteps):
        
            self._updateE(t, dt, self.old)
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