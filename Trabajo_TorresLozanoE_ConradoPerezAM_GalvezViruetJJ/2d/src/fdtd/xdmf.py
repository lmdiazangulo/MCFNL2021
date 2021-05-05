try:
    import lxml.etree as et
except AttributeError:
    pass

import numpy as np

X = 0 # Cartesian indices
Y = 1


class Xdmf:
    def __init__(self, basename, format):
        self.root = et.Element("Xdmf", Version="3.0")
        self.domain = et.SubElement(self.root, "Domain")
        self.basename = basename
        self.format = format

    def tostring(self):
        return et.tostring(self.root, \
                           pretty_print=True, \
                           xml_declaration=True)

    def _createDataItem(self, elem, text, dimensions, \
                        dataType="Float", format="XML"):
        dataItem = et.SubElement(elem, "DataItem", \
                      DataType=dataType, Dimensions=dimensions, Format=format)
        dataItem.text = text
        return dataItem

    def _createAttribute(self, elem, valuesArray, name, \
                         attributeType="Scalar",
                         center="Cell"):
        attr = et.SubElement(elem, "Attribute", \
                             AttributeType=attributeType, \
                             Center=center, \
                             Name=name)
        valueStr =  ' '.join(map(str,  valuesArray.flatten(order='F')))

        if valuesArray.dtype == "float64":
            valueType = "Float"
        elif valuesArray.dtype == "int32" or valuesArray.dtype == "int64":
            valueType = "Int"
        else:
            raise ValueError("Unrecognized dtype when writing DataItem")

        self._createDataItem(attr, valueStr, 
                             dataType=valueType, 
                             dimensions=' '.join(map(str, valuesArray.shape)))
        return attr

    def _createGridWithAttributes(self, frame, probe, timeStep, level):
        mesh   = probe["mesh"]
        values = probe["values"][timeStep]
        Nxy = values.shape
        grid = et.SubElement(frame, "Grid", 
                           Name="L" + str(level), 
                           GridType="Uniform")
        et.SubElement(grid, "Topology", TopologyType="3DCoRectMesh", \
                        Dimensions="2 %d %d"%(Nxy[Y]+1, Nxy[X]+1))
        geom = et.SubElement(grid, "Geometry", GeometryType="ORIGIN_DXDYDZ")
        origin = mesh["origin"]
        self._createDataItem(geom, \
                        "%f %f %f"%(0.0, origin[Y], origin[X]), \
                        dimensions="3")
        steps = mesh["steps"]
        self._createDataItem(geom, \
                        "0.1 " + ' '.join(map(str, mesh["steps"])), \
                        dimensions="3")

        # Values
        self._createAttribute(grid, values, name="Hz")
        
    def add(self, p):
        timeSeries = et.SubElement(self.domain, \
                                   "Grid", \
                                   GridType="Collection", \
                                   CollectionType="Temporal")

        for n in range(len(p["time"])):
            frame = et.SubElement(timeSeries, "Grid", \
                                  Name=p["name"], \
                                  GridType="Collection", \
                                  CollectionType="Spatial")
            et.SubElement(frame, "Time", Value=str(p["time"][n]))
            self._createGridWithAttributes(frame, p, n, 0)
            
