import sys, os, glob, math, argparse, re, pickle, json
import Bio.PDB as bp
from Bio.PDB.PDBParser import PDBParser
import openbabel as ob
from pybel import *
import renumber_pdb as repdb
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict
import pandas as pd



def Args():
	sysargparser = argparse.ArgumentParser()
	sysargparser.add_argument("-d","--directory",type=str,help="Specify folder containing the docking results", required=True)
	sysargparser.add_argument("-t","--target",type=str,help="Specify your protein target", required=True)
	sysargparser.add_argument("-i","--input",type=str,help="Specify (in lower case) the input file format of the docked molecules", required=True)
	sysargs = sysargparser.parse_args()
	return sysargs

def KinaseInfo(structure):
	global targets, use_angle
	if sysargs.target == "EGFR":
		targets = [structure["A"][793]["N"]]
		use_angle = True
	elif sysargs.target == "BRAF":
		targets = [structure["A"][532]["N"],structure["A"][532]["O"],structure["A"][530]["O"]]
		use_angle = False

def GetReceptorAtoms(receptor):
	structure_id = os.path.basename(receptor).replace(".pdb","")
	structures = parser.get_structure(structure_id, receptor)
	structure = structures[0]
	structure_atoms = bp.Selection.unfold_entities(structure, "A")
	KinaseInfo(structure)
	return structure_atoms

def HBondCheck(Hbond_LigAtoms,pose_atoms,posenum,affinity,target):
	if use_angle == True:
		if Hbond_LigAtoms:
			for donor in Hbond_LigAtoms:
				if 140.0 <= AngleFinder(donor,targethydrogen,target) <= 180.0:
					Pass = affinity[posenum]
				else:
					Pass = 0.0
		else:
			Pass = 0.0
	else:
		Pass = 0.0
		for hb_atom in Hbond_LigAtoms:
			for atom in pose_atoms:
				if ("H" in atom.get_id()) and hb_atom - atom <= 1.1:
					Pass = affinity[posenum]
	return Pass

def HydrogenBondFilter(receptor_atoms, molecule,affinity):
	ligand_id = os.path.basename(molecule).replace(".pdb","")
	ligand = parser.get_structure(ligand_id, molecule)
	checklist = []
	global targethydrogen, use_angle
	### Loop through list of all atoms we know contribute to hydrogen bonds from the protein. At least one is essential for function so poses only pass if they form at least one with the defined targets
	for target in targets:
		target_res = target.get_parent()
		protein_hydrogen = [atom for atom in bp.Selection.unfold_entities(target.get_parent(),"A") if ("H" in atom.get_id()) and target - atom <= 1.1]
		if protein_hydrogen:

			use_angle = True
			targethydrogen = protein_hydrogen[0]
		else:
			use_angle = False
		PosePassInfo = []
		### Grab the hydrogen bond info of each pose of each drug. Unfold the atoms of the drug and the combine it with the atoms of the the receptor. Then Search is performed on close atoms
		for posenum, pose in enumerate(ligand):
			pose_atoms = bp.Selection.unfold_entities(pose, "A")
			all_atoms = receptor_atoms + pose_atoms
			ns = bp.NeighborSearch(all_atoms)
			### Search for atoms that are within X angstrom of defined target. Then check whether found atoms are ligand atoms and if they are capable of hydrogen bond formation
			close_atoms = ns.search(target.coord, 3.6) 
			Hbond_LigAtoms = [atom for atom in close_atoms if (atom in pose_atoms) and (atom.get_name().translate(None,"0123456789") in Filter_atoms)]
			### Check if H-bond is within defined parameters per pose and put that in a list per drug. Poses that pass have their ligand efficiency written, poses that don't have 0.0

			Pass = HBondCheck(Hbond_LigAtoms,pose_atoms,posenum,affinity,target)
			PosePassInfo.append(Pass)
		### In case of multiple possible hydrogen bonds a list is kept of poses that have passed at least one of them
		if checklist:
			for posenum, entry in enumerate(PosePassInfo):
				if (entry != 0.0) and (checklist[posenum] == 0.0):
					checklist[posenum] = entry
		else:
			checklist = PosePassInfo

		
	return checklist

def AngleFinder(Atom1,Atom2,Atom3):
	vector1 = Atom1.get_vector()
	vector2 = Atom2.get_vector()
	vector3 = Atom3.get_vector()
	angle = bp.calc_angle(vector1,vector2,vector3)

	return math.degrees(angle)

def FixReceptor(receptor,refseq):
	out_PDB_renum = re.sub('.pdb','',receptor)+"_renumbered.pdb"
	repdb.renumber_noInputAlign(receptor,refseq,selection="chain A",outfile=out_PDB_renum,newAA=None,first=1)
	babelreceptor = readfile("pdb",out_PDB_renum).next()
	babelreceptor.OBMol.AddPolarHydrogens()
	h_receptor = out_PDB_renum.replace(".pdb","_addH.pdb")
	receptorstring = babelreceptor.write("pdb")
	# babelreceptor.write("pdb",h_receptor,overwrite = True)
	sumvar = "999999"
	Hcount = 0
	with open(h_receptor,"w") as pdbfile:
		for line in receptorstring.splitlines():
			if "H " in line:
				resnum = line[22:26]
				if resnum == sumvar:
					Hcount += 1
				else:
					sumvar = resnum
					Hcount = 1

				unique_atom = "H" + str(Hcount).ljust(2)
				newline = line[0:13] + unique_atom + line[16:] + "\n"
				pdbfile.write(newline)
			else:
			 pdbfile.write(line+"\n")
	os.remove(out_PDB_renum)

	return h_receptor

def MakePDB(molecules):
	if ".sdf" in molecules:
		poses = list(readfile("sdf",molecules))
		pdbfilename = molecules.replace(".sdf",".pdb")
	elif ".pdb" in molecules:
		poses = list(readfile("pdb",molecules))
		pdbfilename = molecules.replace(".pdb","_fix.pdb")

	### Normalize affinity data against Heavy Atom acount of corresponding drug. This gives a more comparable value to rank on. (Ligand efficiency). Write these values to a list, index corresponding to pose number
	HatomCount = poses[0].calcdesc()['atoms']
		
	affinitydata = [(float(pose.data["minimizedAffinity"]))/(HatomCount) if pose.data.has_key("minimizedAffinity") else 0 for pose in poses] ### grab affinity info for all conformers from SDF file and save it to list
	atomcountdict = {}
	modelcounter = 1
	### Rewrite PDB file to remove connect records and give each atom a unique name. This is necessary for Biopython to differentiate between them
	with open(pdbfilename,"w") as pdbfile:
		for pose in poses:
			pose.OBMol.AddPolarHydrogens()
			posestring = pose.write("pdb")
			pdbfile.write("MODEL " + str(modelcounter).rjust(8) + "\n")

			for line in posestring.splitlines():
				if "HETATM" in line:
					atomtype = line[13]

					if atomtype in atomcountdict:
						atomcountdict[atomtype] += 1

					else:
						atomcountdict[atomtype] = 1

					unique_atom = str(atomcountdict[atomtype]).ljust(2)
					conformer_name = "M" + str(modelcounter).ljust(2)
					newline = line[0:14] + unique_atom + line[16:17] + conformer_name + line[20:] + "\n"
					pdbfile.write(newline)

			pdbfile.write("ENDMDL\n")
			atomcountdict = {}
			modelcounter += 1

		pdbfile.write("END")
	return affinitydata,pdbfilename

def NewDataParser():
	### Loop through supplied directory and grab all required files as seen below
	for root, dirs, files in os.walk(path):
		ligands = [os.path.join(root, name) for name in files if ("." + sysargs.input in name) and ("osimertinib" not in root) and ("receptor" not in name)]
		for name in files:
			if name == "receptor.pdb":
				receptor = os.path.join(root, name)
				fixed_receptor = FixReceptor(receptor,refseq)
				receptor_atoms = GetReceptorAtoms(fixed_receptor)
				model = os.path.basename(root).upper()

				if ligands:
					for ligand in ligands:
						affin,fixed_ligand = MakePDB(ligand)
						PosePassInfo = HydrogenBondFilter(receptor_atoms, fixed_ligand,affin)
						os.remove(fixed_ligand)
						ligname = os.path.basename(ligand).replace("." + sysargs.input,"")
						modeldictionary[model][ligname] = PosePassInfo
	return modeldictionary

def OldDataParser():
	### Loop through supplied directory and grab all required files as seen below
	for root, dirs, files in os.walk(path):
		ligands = [os.path.join(root, name) for name in files if ("." + sysargs.input in name) and ("receptor" not in name)]
		for name in files:
			if "receptor.pdb" in name:
				receptor = os.path.join(root, name)
				fullname = os.path.basename(root).upper()
				mutation = fullname.split("_EGFR")[0]
				model = fullname.split("-")[1]
				fixed_receptor = FixReceptor(receptor,refseq)
				receptor_atoms = GetReceptorAtoms(fixed_receptor)
				if ligands:
					for ligand in ligands:
						ligname = os.path.basename(ligand).replace("." + sysargs.input,"")
						affin, fixed_ligand = MakePDB(ligand)
						PosePassInfo = HydrogenBondFilter(receptor_atoms,fixed_ligand,affin)
						os.remove(fixed_ligand)
						modeldictionary[mutation][model][ligname] = PosePassInfo
	return modeldictionary


def recursive_dd():
	return defaultdict(recursive_dd)

def PlotHeatmap(dictionary,mutation):
	fig, ax = plt.subplots()
	df = pd.DataFrame.from_dict(dictionary)
	df_swap = df.transpose()
	print df_swap.sum()
	ax.set_title(mutation)
	sns.heatmap(df,vmin=0,vmax=25,cmap="RdYlGn",center=12,annot=True)
	plt.show()

def Rank(modeldictionary,simple):
	for model in modeldictionary.keys():
		maxposenum = 25

		for drug in modeldictionary[model].keys():
			data = modeldictionary[model][drug]
			if maxposenum > len(data):
				for x in range(0, (maxposenum - len(data))):
					data.append(0.0)
				modeldictionary[model][drug] = data
			else:
				maxposenum = len(data)

			passed = np.count_nonzero(data)
			if simple == True:
				modeldictionary[model][drug] = passed
			else:
				data.append(passed)
				modeldictionary[model][drug] = data
			ranked = pd.Series(data).sort_values(ascending=True)[:5]
			topfive = ranked.index.tolist()
			


if __name__ == "__main__":
	pickled = False
	simple = True
	if pickled == True:
		modeldictionary = pickle.load(open("validation_save.p","rb"))
	else:

		parser = PDBParser(PERMISSIVE=1)
		Filter_atoms= ["S","O","F","N"]
		sysargs = Args()
		path = sysargs.directory
		### Grab the correct canonical sequence for the supplied kinase
		if sysargs.target == "EGFR":
			refseq = "../templates/EGFR/EGFR_canonical.fasta"
		elif sysargs.target == "BRAF":
			refseq = "../templates/BRAF/BRAF_canonical.fasta"
		else:
			print "Supplied target not matching template proteins, did you make a typo?"
			exit()
		modeldictionary = recursive_dd()
		modeldictionary = NewDataParser()
		print modeldictionary
		pickle.dump(modeldictionary,open("validation_save.p","wb"))
	

	sumdict = {}
	Rank(modeldictionary,simple)
	PlotHeatmap(modeldictionary, "V600E")

	# for mutation in modeldictionary:
	# 	df = pd.DataFrame.from_dict(modeldictionary[mutation])
	# 	df_swap = df.transpose()
	# 	print df_swap
	# 	print df_swap.axes

	# 	df_swap = df.transpose()
	# 	sumdict[mutation] = df_swap.sum()

	# print sumdict

	# PlotHeatmap(sumdict,"all")


		
	
	# PlotHeatmap(modeldictionary,"G779F")
