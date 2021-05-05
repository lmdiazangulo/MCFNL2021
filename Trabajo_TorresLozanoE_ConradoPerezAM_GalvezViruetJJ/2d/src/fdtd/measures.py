import numpy as np
import copy

from fdtd.common import X, Y, L, U

class Measures:
    def __init__(self, mesh, data, measures):
        self._mesh = copy.deepcopy(mesh)
        self._data = copy.deepcopy(data[0])
        self._measures = copy.deepcopy(measures)
        self.port_inc = self._measures['port_inc']
        self.box_inc = self._mesh.elemIdToBox(self.port_inc["elemId"])
        self.ids_inc = mesh.toIdx(self.box_inc)
        self.port_trans = self._measures['port_trans']
        self.box_trans = self._mesh.elemIdToBox(self.port_trans["elemId"])
        self.ids_trans = mesh.toIdx(self.box_trans)
        self.port_refl = self._measures['port_refl']
        self.box_refl = self._mesh.elemIdToBox(self.port_refl["elemId"])
        self.ids_refl = mesh.toIdx(self.box_refl)

    def Ports(self,n):
        """ Function that calculates the Power per unit of length in the z-axis for the specified port.
        | Input:
        | - n: Integer that specifies the port whose power is measured.
        | Output: 
        | - Power: Power in port.
        """
        if n == 0:
            ids =  self.ids_inc
            x_width = float(ids[U][X]-ids[L][X])
        elif n == 1:
            ids = self.ids_refl
            x_width = float(ids[U][X]-ids[L][X])
        elif n == 2:
            ids = self.ids_trans  
            x_width = float(ids[U][X]-ids[L][X])
        else: raise Exception("0, 1 and 2 ports from left to right")
        self.fields = np.array([{"Hz": np.array([np.array([k for k in j[ids[L][Y]:ids[U][Y]]]) for j in self._data["values"][i][ids[L][X]: ids[U][X]]]) ,\
            "Ex": np.array([np.array([k for k in j[ids[L][Y]:ids[U][Y]]]) for j in self._data["valuese_x"][i][ids[L][X]: ids[U][X]]]) ,\
            "Ey": np.array([np.array([k for k in j[ids[L][Y]:ids[U][Y]]]) for j in self._data["valuese_y"][i][ids[L][X]: ids[U][X]]])}\
            for i in range(0,len(self._data["values"]))])
        self.Power = np.array([sum(sum((-1)*i["Ey"]*i["Hz"]))*self._mesh.dy*(1/x_width) for i in self.fields])    
        return self.Power
        