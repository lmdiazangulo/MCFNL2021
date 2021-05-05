import unittest
import numpy as np
from mesh import Mesh
from solver import FDTD, Utilities, Source

class TestSolver(unittest.TestCase):

    def test_FFT(self):
        
        pulso=Source('gauss',40,12,20)

        #Mesh with a specific permitivity, no conductivity
        et1k1_test1=FDTD(Mesh(200,0.001,4,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)[0]
        et1k2_test1= FDTD(Mesh(200,0.001,4,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)[1]
        et2k1= FDTD(Mesh(200,0.001,1,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)[0]
        et2k2= FDTD(Mesh(200,0.001,1,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)[1]
        #Mesh with no permitivity, no conductivity
        et1k1_test2= FDTD(Mesh(200,0.001,1,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)[0]   
        et1k2_test2= FDTD(Mesh(200,0.001,1,0,110,140), \
                      pulso,5e-9).FDTDLoop(40,160)[1]    


        result1,result2= Utilities().FFT(et1k1_test1,et2k1, et1k2_test1, et2k2)  
        result3,result4= Utilities().FFT(et1k1_test2,et2k1, et1k2_test2, et2k2)  
        

        for i in range(100):
            #Mesh without conductivity should return R*R+T*T=1
            self.assertLess(result1[i], 1)       
            self.assertAlmostEqual(result1[i]*result1[i]+result2[i]*result2[i],1,delta=0.1)
            #Mesh without material should return R=0, T=1
            self.assertEqual(result3[i],0)
            self.assertEqual(result4[i],1)
            






if __name__ == '__main__':
    unittest.main()