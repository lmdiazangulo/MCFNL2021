from __future__ import division

import numpy as np
import copy
import math

from fdtd.common import X, Y, L, U

class Mesh:

    def __init__(self, coordinates, elements, grid):
        self.coordinates = coordinates
        self.elements = elements
        
        if "elemId" in grid:
            box = self.elemIdToBox(grid["elemId"])
        elif "box" in grid:
            box = grid["box"]
        else:
            raise ValueError("Grid data must contain \"elemId\" or \"box\".")

        (Lx, Ly) = abs(box[U] - box[L]) # Gives the longitude of the grid subtracting the corners of the rectangle
        (dx, dy) = grid["steps"]
        (self.dx, self.dy) = (dx, dy)
        
        self.pos =  \
            (np.linspace(box[L][X], box[U][X], num=int(Lx/dx)+1, endpoint=True),
             np.linspace(box[L][Y], box[U][Y], num=int(Ly/dy)+1, endpoint=True) )

        self.bounds = []
        if "bounds" in grid:
            for xy in range(len(grid["bounds"])):
                for lu in range(len(grid["bounds"][xy])):
                    if grid["bounds"][xy][lu] == "pec":
                        self.bounds.append(Mesh.BoundPEC().idsAs(lu, xy))
                    elif grid["bounds"][xy][lu] == "mur":
                        self.bounds.append(Mesh.BoundMUR().idsAs(lu, xy))
    

    def steps(self):
        '''
        Gives a tuple with the real steps of the mesh in the x and y axis, 
        that can vary from the .json gird steps if longitude/step isn't an integer.
        '''
        return (self.pos[X][1]-self.pos[X][0], self.pos[Y][1]-self.pos[Y][0])


    def origin(self):
        '''
        Gives a tuple with the coordinates of the origin of the mesh.
        '''
        return (self.pos[X][0], self.pos[Y][0])


    def elemIdToBox(self, id):
        '''
        Generates a domain for the mesh depending of the type of element given.
        Given a diagonal, generates a rectangle: a tuple with the left bottom 
        corner as first element and the right top corner as second element.
        '''
        return ( np.array(self.coordinates[ self.elements[id][0] ]), \
                 np.array(self.coordinates[ self.elements[id][1] ]) )


    def toIdx(self, coords):
        '''
        Gives the first and last index of both x-axis and y-axis, of the domain 
        given in coords.
        '''
        if type(coords) != tuple and type(coords) != list:
            coords = [coords]
        idx = np.empty((0,2), int)
        for coord in coords:
            nearest = ((np.abs(self.pos[X] - coord[X])).argmin(),
                       (np.abs(self.pos[Y] - coord[Y])).argmin())
            idx = np.vstack((idx, np.asarray(nearest).astype(int)))
        return idx


    def snap(self, coords):
        res = []
        for coord in coords:
            id = self.toIdx(coord)[0]
            res.append( np.array([ self.pos[X][id[X]], self.pos[Y][id[Y]] ]) )
        return res


    def IdxToPos(self, id, pos):
        '''
        Gives the positions on the x-axis and y-axis of the mesh
        that correspond to the array index given in id.
        '''
        pos = np.array(pos, dtype=object)
        return [pos[X][id[L][X]:id[U][X]], \
            pos[Y][id[L][Y]:id[U][Y]]] 

    
    class Bound:
        def __init__(self, ids=None):
            self.ids = ids
            
        def orientation(self):
            for xy in range(2):
                x = xy
                y = (xy + 1) % 2
                if self.ids[L][x] == self.ids[U][x] and \
                   self.ids[L][y] != self.ids[U][y]:
                    return y
            raise ValueError("Error getting orientation")

        def lu(self):
            if (self.ids[L][X], self.ids[U][X]) == (0, 0) or \
               (self.ids[L][Y], self.ids[U][Y]) == (0, 0) :
                return L
            elif (self.ids[L][X], self.ids[U][X]) == (-1, -1) or \
                 (self.ids[L][Y], self.ids[U][Y]) == (-1, -1):
                return U
            else:
                raise ValueError("Error getting lower/upper bound")

        def arrayIdx(self, lu, xy):
            if lu == L:
                return self.ids[lu][xy]
            else:
                if self.ids[U][xy] == -1:
                    return None
                else:
                    return self.ids[U][xy] + 1

        def idsAs(self, lu, xy):
            if xy is X and lu is L:
                self.ids = (np.array([ 0,  0], int), 
                            np.array([ 0, -1], int))
            if xy is X and lu is U:
                self.ids = (np.array([-1,  0], int), 
                            np.array([-1, -1], int))
            if xy is Y and lu is L:
                self.ids = (np.array([ 0,  0], int), 
                            np.array([-1,  0], int))
            if xy is Y and lu is U:
                self.ids = (np.array([ 0, -1], int), 
                            np.array([-1, -1], int))
            
            return self

    class BoundPEC(Bound):
        def __init__(self, ids=None):
            Mesh.Bound.__init__(self, ids)

    class BoundMUR(Bound):
        def __init__(self, ids=None):
            Mesh.Bound.__init__(self, ids)