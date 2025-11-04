import os
import parser
from io import StringIO
import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter
import numpy as np
import sklearn.metrics
import oddt
from oddt import fingerprints
from oddt.toolkits.rdk import Molecule
from Bio.PDB import PDBParser, PDBIO, Select
import time

## ADDED: replacing MDAnalysis PDB output by BioPython (to avoid errors with RDKit) #
standard_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY','HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER','THR', 'TRP', 'TYR', 'VAL'}

# added classes for selections with BioPython
class ProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in standard_residues
				
class LigandOnly(Select):
    def accept_residue(self, residue):
        return residue.get_resname() == "LIG"
#####################################################################################

class score:
    def __init__(self, vars):
        self.__dict__ = vars


    def reference(self):
        if not os.path.exists('reference_protein.pdb') or not os.path.exists('reference_ligand.pdb'):

            # moved to BioPython to avoid errors with ODDT/RDKit #
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("complex", self.solvpdb)
        
            io = PDBIO()
            io.set_structure(structure)
            io.save("reference_protein.pdb", select=ProteinOnly())
            io.save("reference_ligand.pdb", select=LigandOnly())
            ######################################################

        ref = self.ref_fingerprint('reference_protein.pdb', 'reference_ligand.pdb')

        update = {'ref': ref}
        return update



    def ref_fingerprint(self, protein_file, ligand_file):
        protein = next(oddt.toolkit.readfile('pdb', protein_file))
        protein.protein = True

        ligand = next(oddt.toolkit.readfile('pdb', ligand_file))
        
        fp = fingerprints.InteractionFingerprint(ligand, protein, strict=self.strict)

        return fp

        

    def score(self, topology, trajectory, temperature):
        if not os.path.exists('frame_pdbs'):
            os.mkdir('frame_pdbs')
            
        u = mda.Universe(topology, trajectory)
        
        mp_score = []

        for i,ts in enumerate(u.trajectory):
            mp_score.append([u, i])

        outscore = self.parallelizer.run(mp_score, self.calc_ifp, 'Calculating IFPs')

        return outscore

        # 135 rdkit fuori e oddt dentro
        # 85 normal
        # 481 rdkit e oddt dentro

    
    def calc_ifp(self, u, i):
        u.trajectory[i]

        # moved to BioPython to avoid errors with ODDT/RDKit #
        # v0.1 now avoids intermediary files (all on memory)
        protein_file = f'frame_pdbs/protein_{i}.pdb'       
        
        pdb_buffer = StringIO()
        writer = PDBWriter(pdb_buffer)
        writer.write(u.atoms)
        pdb_buffer.seek(0)
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(f"protein_{i}", pdb_buffer)
        writer.close()
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(protein_file, select=ProteinOnly())
        #####################################################

        # original code from here onwards
        u_ligand = u.select_atoms('resname LIG')
        ligand_file = f'frame_pdbs/ligand_{i}.pdb'

        with mda.Writer(ligand_file, u_ligand.n_atoms) as w:
            w.write(u_ligand)

        p = next(oddt.toolkit.readfile('pdb', protein_file))
        p.protein = True

        l = next(oddt.toolkit.readfile('pdb', ligand_file))

        fp = fingerprints.InteractionFingerprint(l, p, strict=self.strict)

        l_plif_temp=[]

        l_plif_temp.append(self.ref)
        l_plif_temp.append(fp)
        matrix = np.stack(l_plif_temp, axis=0)
        idx = np.argwhere(np.all(matrix[..., :] == 0, axis=0))
        matrix_dense = np.delete(matrix, idx, axis=1)
        x = matrix_dense[0].reshape(1,-1)
        y = matrix_dense[1].reshape(1,-1)
        sim_giovanni = float(sklearn.metrics.pairwise.cosine_similarity(x, y))
        sim = round(sim_giovanni * -1,2)

        os.system(f'rm -r {protein_file} {ligand_file}')

        return sim



