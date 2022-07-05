from Bio.PDB import Select, PDBParser, PDBIO
from pathlib import Path
import os

path = Path(input('input the dirctory of Fpocket results: '))
out = Path(input('enter the emity dirctory of the out put: '))

############################function####################################
def aa_pdb(pocket):
	with open(pocket) as pk:
		aas = []	
		for line in pk.readlines():
			if line.startswith('ATOM'):
				aas.append(line.split()[5])
		aas = list(set(aas))
		os.chdir(path / protein)
		protein_name = name + '_out.pdb'
		parser = PDBParser()

		structure = parser.get_structure(name, protein_name)
		class Pocket(Select):
			def accept_residue(self, residue):
				if str(residue.get_id()[1]) in aas:
					return 1
				else:
					return 0
		io = PDBIO()
		io.set_structure(structure)
		pocket_name = name + pocket.replace('pocket','_').replace('_atm', '')
		os.chdir(out)
		io.save(pocket_name, Pocket())
######################################################################
os.chdir(path)
mainlist = os.listdir()
for protein in mainlist:
	name = protein[:6]
	os.chdir(path / protein / 'pockets')
	pockets = list(a for a in os.listdir() if a.endswith('.pdb')) 
	for pocket in pockets:
		os.chdir(path / protein / 'pockets')
		aas = []
		aa_pdb(pocket)
		


