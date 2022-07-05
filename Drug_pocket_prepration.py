
from pathlib import Path
import numpy as np 
import pandas as pd
import os


import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


if os.name == 'nt':
	HOME = Path(r'C:\Users\Mohamed\jupyter_notebooks\researches\biosolveit_project\PocketSimilarityApproach')
else:
	HOME = Path('/content/PocketSimilarityApproach')

path = HOME / Path('DrugPocket')
orginal_path = HOME / Path('OrginalPDBs')
het = HOME / Path('HET_of_approved_drugs.txt')

pdb_pockets = HOME / Path('results')

with open(het,'r') as f:
	drugs_het = f.read().split(',')

###################################################################
cov_rads = {  'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Kr': 1.03, 'X' : 0.00}
# relative atomic masses of elements (in atomic mass units [amu]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_masses = {    'H' : 1.00794, 'D': 1.00794, 'C' : 12.0107, 'O' : 15.9994, 'N' : 14.0067,
  'F' : 18.9984, 'P' : 30.9738, 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
  'I' : 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
  'Be': 9.01218, 'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
  'Si': 28.0855, 'K' : 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
  'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
  'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
  'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X' : 0.00000}
###########################functions###############################
def get_com(at_types, coords):
    n_atoms = len(at_types)
    com = [0.0 for p in range(3)]
    mass = 0.0
    for i in range(n_atoms):
        at_mass = at_masses[at_types[i]]
        mass += at_mass
        for j in range(3):
            com[j] += at_mass * coords[i][j]
    for p in range(3):
        com[p] /= mass
    return com

def aa_pdb(pocket, protein_name):
	print('started')
	os.chdir(path / protein_name / 'pockets')
	with open(pocket) as pk:
		aas = []	
		for line in pk.readlines():
			if line.startswith('ATOM'):
				aas.append(line.split()[5])
		aas = list(set(aas))
		os.chdir(path / protein_name)
		parser = PDBParser()
		protein_name1 = protein_name + '.pdb'
		structure = parser.get_structure(protein_name, protein_name1)
		class Pocket(Select):
			def accept_residue(self, residue):
				if str(residue.get_id()[1]) in aas:
					return 1
				else:
					return 0
		io = PDBIO()
		io.set_structure(structure)
		pocket_name = protein_name.replace('out','') + pocket.replace('pocket','_').replace('_atm', '') + '.pdb'
		os.chdir(pdb_pockets)
		io.save(pocket_name, Pocket())
###################################################################

os.chdir(path)
pdblst = os.listdir()
columns = ['Name','Drug', 'pocket NUM','COM', 'dist']
PocketDF = pd.DataFrame(columns = columns)
for pdb in pdblst:
	os.chdir(path)
	with open('log.txt','a+') as a:
		a.write(pdb + '\n')
	os.chdir(path / pdb / 'pockets')
	pockets = list(a for a in os.listdir() if a.endswith('.pdb'))
	
	aas = []
	for pocket in pockets:
		os.chdir(path / pdb / 'pockets')
		at_types = []
		coords = []
		with open(pocket) as pk:
			aas = []	
			for line in pk.readlines():
				if line.startswith('ATOM'):
					aas.append(line.split()[5])
			aas = list(set(aas))
			os.chdir(path / pdb)
			protein_name = pdb + '.pdb'
			parser = PDBParser()
			name = pdb.replace('_out', '')
			structure = parser.get_structure(name, protein_name)
			for rez in structure.get_residues():
				if str(rez.get_id()[1]) in aas and rez.get_full_id()[3][0] == ' ' :
					for atom in rez.get_atoms():
						atome = atom.full_id[4][0][0]
						coord = atom.get_coord()
						at_types.append(atome) 
						coords.append(coord)
		if 'R' in at_types:
			com = 'not calculated'
		else:
			try:
				com = get_com(at_types, coords)
			except ZeroDivisionError:
				com = 'ZeroDivisionError'
		PocketDF = PocketDF.append({'Name':pdb[:4], 'Drug':pdb[5:8],'pocket NUM': pocket, 'COM':com}, ignore_index=True)
				
os.chdir(path)
PocketDF.to_csv('PocketDF.csv')
##############################getting coms of Drugs###########################################
os.chdir(orginal_path)
pdblst = os.listdir()
columns = ['Name','Drug','COM','cls pocket']
DrugDF = pd.DataFrame(columns = columns)
for pdb in pdblst:
  name = pdb.replace('.pdb','')
  parser = PDBParser()
  structure = parser.get_structure(name, pdb)
  inc = []
  for rez in structure.get_residues():
    at_types = []
    coords = []
    if rez.resname in drugs_het:
      for atom in rez.get_atoms():
        atome = atom.full_id[4][0][0]
        coord = atom.get_coord()
        at_types.append(atome) 
        coords.append(coord)
      if 'R' in at_types:
        com = 'not calculated'
      else:
        com = get_com(at_types, coords)
        DrugDF = DrugDF.append({'Name':pdb[:4], 'Drug':rez.resname, 'COM':com}, ignore_index=True)  
    
DrugDF.to_csv('DrugDF.csv')

##########################################matching distaces##############################################

os.chdir(path)
pockets = pd.read_csv('PocketDf.csv')
os.chdir(orginal_path)
drugs = pd.read_csv('DrugDF.csv')
for indexi, drug in drugs.iterrows():
	dist = 1000
	for index, row in pockets.loc[pockets['Name'] == drug.Name].iterrows():
		pkNum = row['pocket NUM']
		if row.COM != 'ZeroDivisionError':
			if drug.COM != 'not calculated':
				pkCOM = np.array([float(a) for a in row.COM.replace('[','').replace(']','').replace(',','').split() ])
				drCOM = np.array([float(a) for a in drug.COM.replace('[','').replace(']','').replace(',','').split() ])
				diff_vector = drCOM - pkCOM
				new = np.sqrt(np.sum(diff_vector * diff_vector))
				pockets['dist'][index] = new
				if new < dist:
					dist = new
					print(new)
					print(row.Name)
					drugs['cls pocket'][indexi] = pkNum
os.chdir(path)
pockets.to_csv('pockets_dist.csv')
drugs.to_csv('Drugs.csv')


##############################getting pdb form closest pocket#################################
os.chdir(path)
df = pd.read_csv('Drugs.csv')
df.dropna(subset = ['cls pocket'], inplace=True)

for index, row in df.iterrows():
	## Old file nameing
	# protein_name = row.Name +'_' + row.Drug + '_out' 
	# pocket_name = row['cls pocket']
	## New file nameing
	protein_name = row.Name + '_out' 
	pocket_name = row['cls pocket']
	print(protein_name)
	try:
		os.chdir(path / protein_name)
		aa_pdb(pocket_name, protein_name)
	except FileNotFoundError :
		with open('FileNotFoundError.txt', 'a+') as a:
			a.write(protein_name)



