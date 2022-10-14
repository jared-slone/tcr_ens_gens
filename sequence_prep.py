import pandas as pd
from anarci import anarci
import argparse
import os
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.SeqUtils import seq1

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

imgt_loops = [(27,38),(56,66),(105,117)]

def jared_renumber(df,chain,start,stop):
    temp = df.iloc[df.index[df['idx']==start].tolist()[0]-2:df.index[df['idx']==stop].tolist()[0]+3]
    nums = (200+df.index[df['idx']==start].tolist()[0],200+df.index[df['idx']==stop].tolist()[0])
    for index, row in temp.iterrows():
        chain[(' ', row['idx'], ' ')].id = (' ', 200+index, ' ')
    return nums

def imgt_renumber(df,chain):
    #disgusting workaround ... see https://github.com/biopython/biopython/issues/1551
    for residue in chain:
        residue.id = (residue.id[0],500+residue.id[1],residue.id[2])
    #actual freaking code
    i = 0
    for residue in chain:
        if (residue.get_resname() in d.keys()) and (i<len(df.index)) and (d[residue.get_resname()] == df['acid'].iloc[i]):
            residue.id = (residue.id[0],df['idx'].iloc[i],residue.id[2])
            i += 1

def final_renumber(chain,loop_points):
    final_loop_points = []
    i = 0
    for residue in chain:
        if residue.id[1] in loop_points.keys():
            final_loop_points.append(i)
        residue.id = (residue.id[0],i,residue.id[2])
        i+=1
    return final_loop_points

def prep(pdb,alpha_seq=None,beta_seq=None,out_dir=''):

# def prep():

    # parser = argparse.ArgumentParser()
    # parser.add_argument('pdb', help = 'the target pdb file')
    # parser.add_argument('--alpha_seq', dest ='alpha_seq', default = None, help = 'Sequence of TCR alpha chain')
    # parser.add_argument('--beta_seq', dest ='beta_seq', default = None, help = 'Sequence of TCR beta chain')
    # parser.add_argument('--out_dir', dest ='out_dir', default = '', help = 'Directory to save intermediate files')
    # args = parser.parse_args()
    #
    # pdb = args.pdb
    # alpha_seq = args.alpha_seq
    # beta_seq = args.beta_seq
    # out_dir = args.out_dir

    if alpha_seq is None or beta_seq is None:
        pdbparser = PDBParser()
        structure = pdbparser.get_structure('original_structure', pdb)
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
        alpha_seq = chains['D']
        beta_seq = chains['E']

    # #call ANARCI and define loops
    TCR_seq = [("alpha",alpha_seq),("beta",beta_seq)]
    anarci(TCR_seq,scheme='imgt',output=True,outfile=out_dir+'TCR_anarci.txt')
    os.system("tr -s '[:blank:]' ',' < "+out_dir+"TCR_anarci.txt > "+out_dir+"TCR_anarci_stripped.txt")

    renumbered_TCR = pd.read_csv(out_dir+"TCR_anarci_stripped.txt",names=['chain','idx','acid'],comment='#')
    renumbered_TCR = renumbered_TCR[renumbered_TCR['acid'] != '-'].dropna()
    renumbered_TCR_A = renumbered_TCR[renumbered_TCR['chain']=='A'].reset_index()
    renumbered_TCR_B = renumbered_TCR[renumbered_TCR['chain']=='B'].reset_index()

    io = PDBIO()
    parser = PDBParser()
    tcr_pdb = parser.get_structure('original_structure', pdb)

    numsA = []
    imgt_renumber(renumbered_TCR_A,tcr_pdb[0]['D'])
    for loop in imgt_loops:
        numsA.append(jared_renumber(renumbered_TCR_A,tcr_pdb[0]['D'],loop[0],loop[1]))

    numsB = []
    imgt_renumber(renumbered_TCR_B,tcr_pdb[0]['E'])
    for loop in imgt_loops:
        numsB.append(jared_renumber(renumbered_TCR_B,tcr_pdb[0]['E'],loop[0],loop[1]))

    #hard-coded for now
    loop_points_A = {numsA[0][0]:'cdr1',numsA[0][1]:'cdr1_',numsA[1][0]:'cdr2',numsA[1][1]:'cdr2_',numsA[2][0]:'cdr3',numsA[2][1]:'cdr3_'}
    loop_points_B = {numsB[0][0]:'cdr1',numsB[0][1]:'cdr1_',numsB[1][0]:'cdr2',numsB[1][1]:'cdr2_',numsB[2][0]:'cdr3',numsB[2][1]:'cdr3_'}
    final_loop_points_A = final_renumber(tcr_pdb[0]['D'],loop_points_A)
    final_loop_points_B = final_renumber(tcr_pdb[0]['E'],loop_points_B)

    io.set_structure(tcr_pdb)
    io.save('renumbered_'+pdb,preserve_atom_numbering=True)

    return {'a1':(final_loop_points_A[0],final_loop_points_A[1],'D'), 'a2':(final_loop_points_A[2],final_loop_points_A[3],'D'), 'a3':(final_loop_points_A[4],final_loop_points_A[5],'D'),
            'b1':(final_loop_points_B[0],final_loop_points_B[1],'E'), 'b2':(final_loop_points_B[2],final_loop_points_B[3],'E'), 'b3':(final_loop_points_B[4],final_loop_points_B[5],'E')}


# if __name__ == '__main__':
#     prep()
