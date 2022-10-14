import sys
from biopandas.pdb import PandasPdb
import pandas as pd
import os

sys.path.insert(0,'..')
from generate_TCR_ensemble import *
from sequence_prep import *

from openmm.app import *
from openmm import *
from openmm.unit import *

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from biopandas.pdb import PandasPdb

from matplotlib import pyplot as plt
import numpy as np

from time import time
import mdtraj as md
import plotly.express as px


def calculate_rmsd(template_file_name, result_file_name):

	# Load template file
	template_peptide = PandasPdb()
	template_peptide.read_pdb(template_file_name)

	peptide_result = PandasPdb()
	peptide_result.read_pdb(result_file_name)
	# print(len(peptide_result.df['ATOM']),len(template_peptide.df['ATOM']))

	peptide_result_ = PandasPdb()
	template_peptide_ = PandasPdb()

	peptide_result_.df['ATOM'] = peptide_result.df['ATOM'].merge(template_peptide.df['ATOM'][['atom_name','residue_number']], on=['atom_name','residue_number'])
	template_peptide_.df['ATOM'] = template_peptide.df['ATOM'].merge(peptide_result.df['ATOM'][['atom_name','residue_number']], on=['atom_name','residue_number'])
	# print(len(peptide_result.df['ATOM']),len(template_peptide.df['ATOM']))

	ca_rmsd = PandasPdb.rmsd(template_peptide_.df['ATOM'], peptide_result_.df['ATOM'], s='c-alpha') # all atoms, including hydrogens
	backbone_rmsd = PandasPdb.rmsd(template_peptide_.df['ATOM'], peptide_result_.df['ATOM'], s='main chain') # all atoms, including hydrogens
	all_atom_rmsd = PandasPdb.rmsd(template_peptide_.df['ATOM'], peptide_result_.df['ATOM'], s=None) # all atoms, including hydrogens

	return ca_rmsd, backbone_rmsd, all_atom_rmsd

###
# CREATE ENSEMBLES
###
def create_ensembles_for_benchmarking():
	out_dir = 'temp/'
	pdb = 'temp.pdb'
	cdr_loops = prep(pdb,out_dir=out_dir)
	new_pdb = 'renumbered_'+pdb
	prepare_for_RCD(new_pdb,cdr_loops,target_dir=out_dir)
	run_RCD(n_samples=1000,n_keep=10,target_dir=out_dir)
	n_models = None
	loop_pdbs = [out_dir+'rcd_files/renumbered_tempa1_closed.pdb',out_dir+'rcd_files/renumbered_tempa2_closed.pdb',out_dir+'rcd_files/renumbered_tempa3_closed.pdb',
                 out_dir+'rcd_files/renumbered_tempb1_closed.pdb',out_dir+'rcd_files/renumbered_tempb2_closed.pdb',out_dir+'rcd_files/renumbered_tempb3_closed.pdb']
	combine_best_loops(n_models,new_pdb,loop_pdbs,target_dir=out_dir)
	for file in os.listdir(out_dir):
		if file[-4:] == '.pdb':
			renumber_pdb(out_dir+file)
	for file in os.listdir(out_dir):
		if file[-9:] == 'clean.pdb':
			add_sidechains(out_dir+file)
			minimizeEnergy(out_dir+file[:-4]+'sidechains.pdb')

def create_cyclic_ensembles_for_benchmarking():
	out_dir = 'temp/'
	pdb = 'input_equi_1.pdb'
	cdr_loops = prep(pdb,out_dir=out_dir)
	new_pdb = 'renumbered_'+pdb
	# cyclic_RCD(new_pdb,cdr_loops,out_dir=out_dir)
	stocastic_RCD(new_pdb,cdr_loops,out_dir=out_dir)
	print(['*']*100)
	# for file in os.listdir(out_dir):
	# 	if file[-4:] == '.pdb':
	# 		renumber_pdb(out_dir+file)
	# os.chdir('..')
	for file in os.listdir(out_dir):
		if file[-9:] == 'hains.pdb':
			# add_sidechains(out_dir+file)
			minimizeEnergy(out_dir+file)

def calculate_stats():
	out_dir = 'temp/'
	out_f = out_dir+'energy.txt'
	template_file_name = out_dir+'renumbered_tempcleansidechains_fin.pdb' #TRY with RENUMBERD pdb
	with open(out_f,'w+') as of:
		for file in os.listdir(out_dir):
			if file[-8:] == '_fin.pdb':

				energy = get_comfirmation_energy(out_dir+file)
				ca_rmsd, backbone_rmsd, all_atom_rmsd = calculate_rmsd(template_file_name,out_dir+file)

				of.write(file+','+str(energy)+','+str(ca_rmsd)+','+str(backbone_rmsd)+','+str(all_atom_rmsd)+'\n')
	res = np.genfromtxt(out_f,delimiter=',')
	plt.scatter(res[:,1],res[:,2])
	plt.show()

def get_md_rmsd(pdb):
	tcr_traj = md.load('../../../5brz_md/rep0/full_traj_aligned.dcd',top='../../../5brz_md/input_equi_1.pdb')
	ref_comformation = md.load(pdb)
	ind_ref = ref_comformation.topology.select_atom_indices('heavy')
	ind_traj = tcr_traj.topology.select_atom_indices('heavy')
	tcr_traj = tcr_traj.atom_slice(ind_traj)
	ref_comformation = ref_comformation.atom_slice(ind_ref)
	rmsds = md.rmsd(tcr_traj,ref_comformation,0)
	fig = px.line(x=range(tcr_traj.n_frames), y=rmsds, title='RMSDs to the reference')
	fig.show()
