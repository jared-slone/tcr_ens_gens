import os
import shutil
import numpy as np
import pandas as pd
from sys import stdout
import copy
from time import sleep

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from openmm.app import *
from openmm import *
from openmm.unit import *

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from biopandas.pdb import PandasPdb

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def prepare_for_RCD(pdb,cdr_loops,target_dir=''):
    #assumes pdb file is already acceptably numbered
    #pdb : file path to tcr pdb
    #cdr_loops : dictionary of the form {loop name:indicies} indicies=[start,stop,chain]
    io = PDBIO()
    parser = PDBParser()
    tcr_structure = parser.get_structure('initial_structure',pdb)

    # print(pdb)

    for loop in cdr_loops.keys():
        shutil.copy(pdb,pdb[:-4]+loop+'.pdb')

        with open(target_dir+'rcd_targets.txt','a+') as out_file:
            out_file.write(pdb[:-4]+loop+'.pdb '+str(cdr_loops[loop][0])+' '+str(cdr_loops[loop][1])+' '+cdr_loops[loop][2]+' ')
            for i in range(cdr_loops[loop][0],cdr_loops[loop][1]+1):
                out_file.write(d[tcr_structure[0][cdr_loops[loop][2]][i].get_resname()])
            out_file.write('\n')

def run_RCD(n_cores=4,n_samples=1000,n_keep=100,dd=0.5,t=0.90,target_dir='',rcd_file_location='../RCD_required_files',target_file_name='rcd_targets.txt'):
    os.system('mpirun.openmpi -np '+str(n_cores)+' '+rcd_file_location+'/bin/rcd_mpi_gnu '+target_dir+target_file_name+' -n '+str(n_samples)+' -r -t '+str(t)+' -d '+str(dd)+' --linear -x '+rcd_file_location+'/dunbrack.bin --energy_file '+rcd_file_location+'/korp6Dv1.bin --loco_best '+str(n_keep)+' --bench -o '+target_dir+'rcd_files')

def cyclic_RCD(pdb,cdr_loops,iterations=10,n_cores=4,step_samples=200,n_keep=1,dd=0.5,t=0.90,rcd_file_location='../../RCD_required_files',target_dir='',target_file_name='rcd_targets.txt',out_dir=''):
    shutil.copyfile(pdb, out_dir+pdb)
    os.chdir(out_dir)

    io = PDBIO()
    parser = PDBParser()
    tcr_structure = parser.get_structure('initial_structure',pdb)

    loop_location_strings = []
    for loop in cdr_loops.keys():
        temp_str = ' '+str(cdr_loops[loop][0])+' '+str(cdr_loops[loop][1])+' '+cdr_loops[loop][2]+' '
        for i in range(cdr_loops[loop][0],cdr_loops[loop][1]+1):
            temp_str = temp_str+d[tcr_structure[0][cdr_loops[loop][2]][i].get_resname()]
        loop_location_strings.append(temp_str)

    # print(loop_location_strings)
    # sleep(15)

    pdb_ = copy.deepcopy(pdb)
    for i in range(iterations):
        n = 0
        for loop_string in loop_location_strings:
            with open(target_dir+target_file_name,'w') as of:
                of.write(pdb_+loop_string)
            run_RCD(n_cores=n_cores,n_samples=step_samples,n_keep=n_keep,dd=dd,t=t,rcd_file_location=rcd_file_location,target_dir=target_dir,target_file_name=target_file_name) #TODO: increse number of samples for cdr3 loops

            loop_model = parser.get_structure('loop',target_dir+'rcd_files/'+pdb_[:-4]+'_closed.pdb')[0]
            for residue in loop_model[list(loop_model.get_chains())[0].id]:
                del tcr_structure[0][list(loop_model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(loop_model.get_chains())[0].id].add(residue)
            io.set_structure(tcr_structure)
            io.save(pdb[:-4]+str(i)+str(n)+'.pdb',preserve_atom_numbering=False)
            pdb_ = pdb[:-4]+str(i)+str(n)+'.pdb'
            renumber_pdb(pdb_)
            pdb_ = pdb_[:-4]+'clean.pdb'
            add_sidechains(pdb_)
            pdb_ = pdb_[:-4]+'sidechains.pdb'
            n += 1

            # print([' '+str(n)+' ']*100)
            # sleep(5)

def stocastic_RCD(pdb,cdr_loops,iterations=100,n_cores=4,step_samples=200,n_keep=1,dd=0.5,t=0.90,rcd_file_location='../../RCD_required_files',target_dir='',target_file_name='rcd_targets.txt',out_dir=''):
    shutil.copyfile(pdb, out_dir+pdb)
    os.chdir(out_dir)

    io = PDBIO()
    parser = PDBParser()
    tcr_structure = parser.get_structure('initial_structure',pdb)

    loop_location_strings = []
    for loop in cdr_loops.keys():
        temp_str = ' '+str(cdr_loops[loop][0])+' '+str(cdr_loops[loop][1])+' '+cdr_loops[loop][2]+' '
        for i in range(cdr_loops[loop][0],cdr_loops[loop][1]+1):
            temp_str = temp_str+d[tcr_structure[0][cdr_loops[loop][2]][i].get_resname()]
        loop_location_strings.append(temp_str)

    # Assume that cdr3 loops are at indecies 2 & 5 of cdr_loops
    p = [1/12,1/12,1/3,1/12,1/12,1/3] #weighting the cdr3 loops by 4x becasue that's what TCR_dist does.. so why not?

    pdb_ = copy.deepcopy(pdb)
    for i in range(iterations):
        loop_string = loop_location_strings[np.random.choice(6,p=p)]
        with open(target_dir+target_file_name,'w') as of:
            of.write(pdb_+loop_string)
        run_RCD(n_cores=n_cores,n_samples=step_samples,n_keep=n_keep,dd=dd,t=t,rcd_file_location=rcd_file_location,target_dir=target_dir,target_file_name=target_file_name) #TODO: increse number of samples for cdr3 loops

        loop_model = parser.get_structure('loop',target_dir+'rcd_files/'+pdb_[:-4]+'_closed.pdb')[0]
        for residue in loop_model[list(loop_model.get_chains())[0].id]:
            del tcr_structure[0][list(loop_model.get_chains())[0].id][residue.id]
            tcr_structure[0][list(loop_model.get_chains())[0].id].add(residue)
        io.set_structure(tcr_structure)
        io.save(pdb[:-4]+str(i)+'.pdb',preserve_atom_numbering=False)
        pdb_ = pdb[:-4]+str(i)+'.pdb'
        renumber_pdb(pdb_)
        pdb_ = pdb_[:-4]+'clean.pdb'
        add_sidechains(pdb_)
        pdb_ = pdb_[:-4]+'sidechains.pdb'

    os.chdir('..')

def combine_loops(n_models,initial_pdb,loop_pdbs,target_dir=''):
    io = PDBIO()
    parser = PDBParser()
    tcr_structure = parser.get_structure('initial_structure',initial_pdb)

    loop_structures = [parser.get_structure(pdb[:-4],pdb) for pdb in loop_pdbs]
    #TODO, try different ways of comnining
    for i in range(n_models):
        for structure in loop_structures:
            model = structure[i%len(structure)]
            # for residue in model[repr(list(model.get_chains())[0])]:
            for residue in model[list(model.get_chains())[0].id]:
                # print(residue)
                # print(tcr_structure[0][list(model.get_chains())[0].id][residue.id])
                # tcr_structure[0][list(model.get_chains())[0].id][residue.id] = residue
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)
        #TODO reorder chains so residues are next to eachother
        io.set_structure(tcr_structure)
        io.save(target_dir+initial_pdb[:-4]+str(i)+'.pdb',preserve_atom_numbering=False)

def combine_best_loops(n_models,initial_pdb,loop_pdbs,target_dir=''):
    io = PDBIO()
    parser = PDBParser()
    tcr_structure = parser.get_structure('initial_structure',initial_pdb)

    loop_structures = [parser.get_structure(pdb[:-4],pdb) for pdb in loop_pdbs]
    for a3 in range(10):
        for b3 in range(10):
            a1 = np.random.randint(9)
            a2 = np.random.randint(9)
            b1 = np.random.randint(9)
            b2 = np.random.randint(9)

            model = loop_structures[0][a1]
            for residue in model[list(model.get_chains())[0].id]:
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)

            model = loop_structures[1][a2]
            for residue in model[list(model.get_chains())[0].id]:
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)

            model = loop_structures[2][a3]
            for residue in model[list(model.get_chains())[0].id]:
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)

            model = loop_structures[3][b1]
            for residue in model[list(model.get_chains())[0].id]:
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)

            model = loop_structures[4][b2]
            for residue in model[list(model.get_chains())[0].id]:
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)

            model = loop_structures[5][b3]
            for residue in model[list(model.get_chains())[0].id]:
                del tcr_structure[0][list(model.get_chains())[0].id][residue.id]
                tcr_structure[0][list(model.get_chains())[0].id].add(residue)

            io.set_structure(tcr_structure)
            io.save(target_dir+initial_pdb[:-4]+str(a3)+str(b3)+'.pdb',preserve_atom_numbering=False)

def get_comfirmation_energy(pdb):
    # Read PDB
    pdb = PDBFile(pdb)
    top = pdb.getTopology()
    positions = np.array(pdb.positions)
    numAtoms = len(positions)
    positions = np.reshape(positions, (3*numAtoms,1))

    # Create the ForceField
    # forcefield = ForceField('amber14-all.xml', 'amber/tip3pfb.xml')
    forcefield = ForceField('amber/ff14SB.xml', 'amber/phosaa14SB.xml')
    # forcefield = ForceField('charmm/charmm36_nowaters.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=None)

    # Rest (Integrator + Platform)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    platform = Platform.getPlatformByName('CPU')

    # Create Simulation
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kilojoule_per_mole
    return energy

def renumber_pdb(pdb_file):
    with open(pdb_file, 'r') as in_file:
        lines = in_file.readlines()
    D_lines = [line for line in lines if((line[:4]=="ATOM")and(len(line)>=26)and(line[21]=='D'))]
    E_lines = [line for line in lines if((line[:4]=="ATOM")and(len(line)>=26)and(line[21]=='E'))]
    D_lines.sort(key = lambda x: int(x[22:26].strip()))
    E_lines.sort(key = lambda x: int(x[22:26].strip()))
    with open(pdb_file[:-4]+'clean.pdb', 'w') as out_file:
        for line in D_lines:
            out_file.write(line)
        for line in E_lines:
            out_file.write(line)
        out_file.write("END\n")

def add_sidechains(pdb_filename, add_hydrogens=True, remove_heterogens=False, add_solvent=False, keep_IDs=True):
    fixer = PDBFixer(filename=pdb_filename)
    fixer.findMissingResidues()
    if remove_heterogens:
        fixer.removeHeterogens(True)#  True keeps water molecules while removing all other heterogens, REVISIT!
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(7.0) # Ask Mauricio about those
    if add_solvent:
        fixer.addSolvent(fixer.topology.getUnitCellDimensions()) # Ask Mauricio about those

    with open(pdb_filename[:-4]+'sidechains.pdb', 'w') as out_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_file, keepIds=keep_IDs)

def minimizeConf(pdb_filename, best_energy, device='CPU'):
    # Read PDB
    pdb = PDBFile(pdb_filename)
    top = pdb.getTopology()
    positions = np.array(pdb.positions)
    numAtoms = len(positions)
    positions = np.reshape(positions, (3*numAtoms,1))

    # Create the ForceField
    forcefield = ForceField('amber14-all.xml', 'amber/tip3pfb.xml')
    # forcefield = ForceField('amber/ff14SB.xml', 'amber/phosaa14SB.xml')
    # forcefield = ForceField('charmm/charmm36_nowaters.xml')
    # forcefield = ForceField('charmm36.xml')#,'implicit/obc2.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=None,hydrogenMass=4)

    # Adding Forces?
    force_constant = 1000
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    # protein_particles = md.load(filename).top.select("backbone")

    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)
    pdb_df = ppdb.df['ATOM']
    protein_particles = [atomind - 1 for atomind in pdb_df[pdb_df['atom_name'].isin(["N", "O", "C", "CA"])]['atom_number'].tolist()]
    particle_indices = []
    for protein_particle in protein_particles:
        particle_indices.append(force.addParticle(int(protein_particle), modeller.positions[protein_particle]) )
    system.addForce(force)

    # Enumerate forces?
    forces = system.getForces()
    for i, f in enumerate(forces):
        f.setForceGroup(i)

    # Rest (Integrator + Platform)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    platform = Platform.getPlatformByName(device)

    # Create Simulation
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # Minimize energy
    simulation.minimizeEnergy()
    simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True,
                                temperature=True, progress=False, remainingTime=True, speed=True,
                                totalSteps=250000, separator='\t'))

    # Write results to a new file if energy is small enough
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kilojoule_per_mole
    if energy < best_energy:
        best_energy = energy
        # path, file = os.path.split(pdb_filename)
        r = PDBReporter(pdb_filename[:-4]+'_fin.pdb', 1)
        r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))

    return best_energy

def minimizeEnergy(pdb_filename,num_tries=1,out_file='energias.txt'):
    best_energy = float("inf")
    for minimization_effort in range(num_tries):
        best_energy = minimizeConf(pdb_filename, best_energy)
    with open(out_file,'a+') as of:
        of.write(pdb_filename+','+str(best_energy)+'\n')

#TODO: filter out comformations with steric clashes
#only check loop clashes, distance of <0.35 nm
def detect_clash(model1,model2,threshold=0.1225):
    coords1 = []
    for chain in model1:
        for residue in chain:
            for atom in residue:
                coords1.append(atom.get_vector().get_array())
    coords1 = np.array(coords1)
    coords2 = []
    for chain in model2:
        for residue in chain:
            for atom in residue:
                coords2.append(atom.get_vector().get_array())
    coords2 = np.array(coords2)
    min_d = lambda x,Y: np.min(np.sum((Y-x)**2,1))
    for coor in coords1:
        if min_d(coor,coords2) < threshold:
            return True
    return False

# def test_temp():
#     io = PDBIO()
#     parser = PDBParser()
#     tcr_structure1 = parser.get_structure('initial_structure','renumbered_tempb2_closed.pdb')
#     tcr_structure2 = parser.get_structure('initial_structure2','renumbered_tempb3_closed.pdb')
#     for i in range(100):
#         for j in range(100):
#             print(detect_clash(tcr_structure2[i],tcr_structure1[j]))
#             print(detect_clash(tcr_structure2[i],tcr_structure2[j]))

def test_stuff1():
    pdb = 'renumbered_TCR.pdb'
    cdr_loops = {'a1':(26,31,'D'), 'a2':(49,55,'D'), 'a3':(90,102,'D'),
                 'b1':(25,29,'E'), 'b2':(47,52,'E'), 'b3':(90,100,'E')}
    prepare_for_RCD(pdb,cdr_loops)

def test_stuff2():
    n_models = 10
    initial_pdb = 'renumbered_TCR.pdb'
    loop_pdbs = ['test_run/renumbered_TCRa1_closed.pdb','test_run/renumbered_TCRa2_closed.pdb','test_run/renumbered_TCRa3_closed.pdb',
                 'test_run/renumbered_TCRb1_closed.pdb','test_run/renumbered_TCRb2_closed.pdb','test_run/renumbered_TCRb3_closed.pdb']
    combine_loops(n_models,initial_pdb,loop_pdbs)

def test_stuff3():
    filename = 'test/renumbered_TCR0.pdb'
    renumber_pdb(filename)

def test_stuff4():
    filename = 'test/renumbered_TCR0clean.pdb'
    add_sidechains(filename)

def test_stuff5():
    filename = 'test/renumbered_TCR0cleansidechains.pdb'
    minimizeEnergy(filename)
