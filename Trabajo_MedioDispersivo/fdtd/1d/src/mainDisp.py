import json
import argparse
import os.path
import sys

from fdtd.mesh import Mesh
from fdtd.DispMedia import SolverComplete
from fdtd.PartialDisp import SolverPartial
from fdtd.solver import Solver 
from fdtd.viewerDisp import Animator

import matplotlib.pyplot as plt

print("=== Python FDTD 1D")

parser = argparse.ArgumentParser(description='Python FDTD 1D')
parser.add_argument('-i', '--input', nargs=1, type=str)
args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

inputFilename = ''.join(args.input).strip()
print("--- Reading file: %s"%(inputFilename))
data = json.load(open(inputFilename))

print('--- Initializing mesh')
mesh = Mesh(data["coordinates"], data["elements"], data["grid"])

print('--- Initializing solver')
solver = Solver(mesh, data["options"], data["probes"], data["sources"], data)
solverComp = SolverComplete(mesh, data["options"], data["probes"], data["sources"], data)
solverPart = SolverPartial(mesh, data["options"], data["probes"], data["sources"], data)

print('--- Solving')
#solver.solve(data["options"]["finalTime"])
#solverComp.solve(data["options"]["finalTime"])
solverPart.solve(data["options"]["finalTime"])

print('--- Visualizing')
#Animator(mesh, solver.getProbes()[0])
#Animator(mesh, solverComp.getProbes()[0])
Animator(mesh, solverPart.getProbes()[0])

print('=== Program finished')

