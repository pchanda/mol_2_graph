import re
from rdkit import Chem
import rdkit


def get_mol_from_graph(line) :

        x = line.split('\t')
        mol_smiles = x[1]
        mol_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(mol_smiles),True) #make sure this is cannonical
        atom_info = x[2]
        bond_info = x[3]
        regex = r"\{(.*?)\}"
        matches = re.finditer(regex, atom_info , re.MULTILINE | re.DOTALL)
        atom_info_list = []
        for matchNum, match in enumerate(matches):
            for groupNum in range(0, len(match.groups())):
                atom_info_list.append(match.group(1))

        
        #create mol from the graph
        mol = Chem.RWMol()
        natoms = int(atom_info_list[0])
        node_to_idx = {}
        idx_to_info = {}
        for atom_index in range(natoms):
            atom_info = atom_info_list[atom_index + 1]
            z = atom_info.split(':')
            ID = z[0]
            typ = z[1]
            smiles = z[2]
            order = z[3]
            sym = z[4]
            at_no = z[5]
            f_c = z[6]
            hyb = z[7]
            i_hc = z[8]
            e_hc = z[9]
            is_aro = z[10]

            a=Chem.Atom(int(at_no))
            #a.SetChiralTag(chiral_tags[node])
            a.SetFormalCharge(int(f_c))
            if is_aro=='true':
                a.SetIsAromatic(True)
            else:   
                a.SetIsAromatic(False)
            if hyb in rdkit.Chem.rdchem.HybridizationType.names:
                a.SetHybridization(rdkit.Chem.rdchem.HybridizationType.names[hyb])
            else:
                a.SetHybridization(rdkit.Chem.rdchem.HybridizationType.names['OTHER'])
            #a.SetHybridization(hyb_map[hyb])
            a.SetNumExplicitHs(int(e_hc))
            idx = mol.AddAtom(a)
            node_to_idx[ID] = idx
            idx_to_info[idx] = [ID,typ,smiles,order]
            print(idx,ID,typ,smiles,order)
        
        bond_info_arr = bond_info.split(',')
        nbonds = int(bond_info_arr[0])
        for bond_index in range(nbonds):
            bond = bond_info_arr[bond_index + 1]
            y = bond.split(':')
            edge = y[0]
            btype = y[1]
            y = edge.split('-')
            atom1_id = y[0]
            atom2_id = y[1]
            atom1_idx = node_to_idx[atom1_id]
            atom2_idx = node_to_idx[atom2_id]
            if(atom1_idx < atom2_idx):
                #print('Adding',atom1_idx,atom2_idx,btype)
                mol.AddBond(atom1_idx,atom2_idx,rdkit.Chem.rdchem.BondType.names[btype])

        Chem.SanitizeMol(mol)   

        new_smiles = Chem.MolToSmiles(mol)

        print('smiles=',mol_smiles,'\n')
        print('new_smiles=',new_smiles,'\n')
        assert new_smiles == mol_smiles
        print('graph and molecule matched')
        return mol,mol_smiles,idx_to_info


#with open('cluster_1_graphs.txt','r') as ff:
#    for line in ff:
#        mol = get_mol_from_graph(line)
