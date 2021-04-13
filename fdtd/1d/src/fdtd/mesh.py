import numpy as np
import copy
import math

L = 0 # Lower
U = 1 # Upper

class Mesh:

    def __init__(self, coordinates, elements, grid):
        self.coordinates = coordinates
        self.elements = elements
        
        if "elemId" in grid:
            box = self.elemIdToBox(grid["elemId"])
        else:
            raise ValueError("Grid data must contain \"elemId\" or \"box\".")

        Lx = abs(box[U] - box[L])
        dx = grid["steps"]
        self.pos =  np.linspace(box[L], box[U], num=int(Lx/dx)+1, endpoint=True)
        self.bounds = grid["bounds"]

    def steps(self):
        return self.pos[1]-self.pos[0]  # Assumes everything is equispaced.

    def elemIdToBox(self, id):
        return ( self.coordinates[ self.elements[id][0] ], \
                 self.coordinates[ self.elements[id][1] ] )

    def toIds(self, coords):
        if type(coords) != tuple and type(coords) != list:
            coords = np.array(coords)
        ids = np.empty((0,1), int)
        for coord in np.unique(coords):
            nearest = (np.abs(self.pos - coord)).argmin()
            ids = np.append(ids, nearest)
        return ids

    def snap(self, coords):
        res = []
        for coord in coords:
            id = np.array(self.toIds(coord)[0], dtype=int)
            res.append( np.array([ self.pos[id] ]) )
        return res