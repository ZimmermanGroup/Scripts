#TO USE
#need to set up pygsm
#do: module load zimmerman; module load pygsm
#this uses python 3.5, and anaconda

from pyGSM.utilities import manage_xyz as mxyz
from pyGSM.utilities import elements 
from pyGSM.coordinate_systems import topology as top
from sys import argv
from collections import deque

infile=argv[1]

def countAtom(element, connection):
  counter=0
  for item in connection[1:]:
    if item[0] == element:
      counter=counter+1
  return counter

coords = mxyz.read_xyz(infile) # read in xyz
geom = mxyz.xyz_to_np(coords) # convert cartesians to np array
atom_symbols = mxyz.get_atoms(coords)
element_table = elements.ElementData()
atoms = [element_table.from_symbol(atom) for atom in atom_symbols]

topology=top.Topology.build_topology(geom, atoms)
bonds = topology.edges() #gives all bonding relationships
elements = topology.e() #gives elemental symbol for each atom (elements[5] gives the element of 5)


print("building connections from "+infile)

connections=[]
#connections should look like:
#[[atom_name,[element, index], [element, index], ...] ...]

for i,atom in enumerate(atom_symbols):
  atom_connect = [atom]
  for bond in bonds:
    if i in bond:
      for a in bond:
        if a != i:
          atom_connect.append([elements[a],a])
  connections.append(atom_connect)

#print("built connections")

molecGraph = {}

for line in enumerate(connections):
  centrAtom = line[0]+1 #getting rid of 0-indexing
  adjacentAtoms = (line[1][1:]) #getting near atoms
  indexes = []
  for i in range(len(adjacentAtoms)):
    atomIndex = adjacentAtoms[i][1] + 1
    indexes.append(atomIndex)
  molecGraph[centrAtom] = indexes

#print(molecGraph)

def find_atomsforQMMM(inGraph, startAtoms, borderAtoms):
  visited = []
  queue = deque([])
  boundaries = []
  for atom in startAtoms:
    queue.append(int(atom))
  for atom in borderAtoms:
    boundaries.append(int(atom))
  
  while queue:
    node = queue.popleft()
    if node not in visited:
      visited.append(node)
      neighbours = inGraph[node]
      for neighbour in neighbours:
        if neighbour not in boundaries:
          queue.append(neighbour)
  return visited

#you need specify at least one QM atom in molecule which fully in QM to get QM atoms for this mol also
print("Specify one QM atom for molecule which is fully in QM")
qmAtomStart = (input("Enter QM border atom indexes: ")).split()
mmAtomStart = (input("Enter MM border atom indexes: ")).split()

qmAtoms = find_atomsforQMMM(molecGraph, qmAtomStart, mmAtomStart)
mmAtoms = find_atomsforQMMM(molecGraph, mmAtomStart, qmAtomStart)

qmAtoms.sort()
mmAtoms.sort()

print("QM atoms: ")
outString = ''
for atom in qmAtoms:
  outString += str(atom) + " "
print(outString)
