import sys, os, math, argparse, re, pickle
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB as bp
import openbabel as ob
from pybel import *
import renumber_pdb as repdb
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,TextBox
import numpy as np
import seaborn as sns
from collections import defaultdict
import pandas as pd


def Args():
	sysargparser = argparse.ArgumentParser()
	sysargparser.add_argument("-d","--directory",type=str,help="Specify folder containing the docking results", required=True)
	sysargparser.add_argument("-r","--refseq_folder",type=str, help="Specify folder containing the canonical sequences of the Kinases",required=True)
	sysargparser.add_argument("-t","--target",type=str,help="Specify your protein target e.g. EGFR", required=True)
	sysargparser.add_argument("-i","--input",type=str,help="Specify (in lower case) the 3-letter input file format of the docked molecules e.g. sdf", required=True)
	sysargparser.add_argument("-s","--sort",help="use flag if you want to sort data by most passed poses",action="store_true")
	sysargparser.add_argument("-a","--appendix",type=str,help="Specify suffix used to differentiate docked molecules from non-docked e.g. _dock", required=True)
	sysargparser.add_argument("-p","--pickled",help="Run this flag after analysing the data at least once, this flag will load in the previously analysed data without rerunning the analysis", action="store_true")
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
	if "sdf" in sysargs.input:
		affinitydata = [(float(pose.data["minimizedAffinity"]))/(HatomCount) if pose.data.has_key("minimizedAffinity") else 0 for pose in poses] ### grab affinity info for all conformers from SDF file and save it to list
	else:
		affinitydata = []
	atomcountdict = {}
	modelcounter = 1
	### Rewrite PDB file to remove connect records and give each atom a unique name. This is necessary for Biopython to differentiate between them
	with open(pdbfilename,"w") as pdbfile:
		for pose in poses:
			pose.OBMol.AddPolarHydrogens()
			posestring = pose.write("pdb")
			pdbfile.write("MODEL " + str(modelcounter).rjust(8) + "\n")

			for line in posestring.splitlines():
				if ("pdb" in sysargs.input) and ("minimizedAffinity" in line):
					affinitydata.append((float(line[25:]))/(HatomCount))
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
		ligands = [os.path.join(root, name) for name in files if (sysargs.appendix + "." + sysargs.input in name) and ("osimertinib" not in root) and ("receptor" not in name)]
		for name in files:
			if name == "receptor.pdb":
				receptor = os.path.join(root, name)
				fixed_receptor = FixReceptor(receptor,refseq)
				receptor_atoms = GetReceptorAtoms(fixed_receptor)
				model = os.path.basename(root).upper()

				if ligands:
					for ligand in ligands:
						print ligand
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

def ConstructDataframe(dictionary):
	df = pd.DataFrame.from_dict(dictionary)
	if sysargs.sort:
		total = df.transpose().sum()/len(df.columns)
		df['Average'] = total
		df_sorted = df.sort_values(by=["Average"])
		return df_sorted
	return df


def PlotHeatmap(dictionary,mutation):
	global spos, ax, ax_basic,molnum,fig, winval
	winval = "1"
	fig, ax = plt.subplots()
	plotdata = ConstructDataframe(dictionary)
	yaxis = list(plotdata.index)
	ax = sns.heatmap(plotdata,cmap="RdYlGn",yticklabels=yaxis,vmax = 1.0,linewidths=0.1,center=0.5,annot=True)
	ax.set_title(mutation)
	for label in ax.get_yticklabels():
		label.set_rotation(0)
	ax_basic = ax.axis()
	sliderbox = [0.2, 0.02, 0.12, 0.02]
	axpos_slider = plt.axes(sliderbox, facecolor='yellow')
	winvalbox = [0.5, 0.02, 0.12, 0.02]
	axpos_winval = plt.axes(winvalbox, facecolor='yellow')
	swinval = TextBox(axpos_winval, 'Window',initial=winval)
	ax.axis([ax_basic[0], ax_basic[1],0, int(winval)])
	fig.canvas.draw_idle()
	spos = Slider(axpos_slider, 'Pos', 0, len(plotdata.index),valinit=0)
	swinval.on_submit(windowupdate)
	spos.on_changed(update)
	plt.show()
### Update events to make Heatmap scrollable
def update(val):
	pos = int(spos.val)
	ax.axis([ax_basic[0], ax_basic[1],pos, pos+int(winval)]) # [xmin xmax ymin ymax]
	fig.canvas.draw_idle()
	return()
def windowupdate(text):
	global winval
	winval = int(text)
	pos = int(spos.val)
	ax.axis([ax_basic[0], ax_basic[1],pos, pos+int(winval)]) # [xmin xmax ymin ymax]
	fig.canvas.draw_idle()
	return()

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
				modeldictionary[model][drug] = (float(passed)/float(len(data)))
			else:
				data.append(passed)
				modeldictionary[model][drug] = data
			ranked = pd.Series(data).sort_values(ascending=True)[:5]
			topfive = ranked.index.tolist()
			


if __name__ == "__main__":
	sysargs = Args()
	simple = True
	if sysargs.pickled == True:
		modeldictionary = pickle.load(open("validation_save.p","rb"))
		sysargs = Args()
	else:

		parser = PDBParser(PERMISSIVE=1)
		Filter_atoms= ["S","O","F","N"]
		path = sysargs.directory
		### Grab the correct canonical sequence for the supplied kinase
		canon_files = [f for f in os.listdir(sysargs.refseq_folder) if (os.path.isfile(os.path.join(sysargs.refseq_folder, f))) and (sysargs.target in os.path.basename(f))]
		print canon_files
		if canon_files:
				refseq = os.path.join(sysargs.refseq_folder,canon_files[0])
		else:		
			print "Cannot find canonical sequence, is your specified folder correct?"
			exit()
		modeldictionary = recursive_dd()
		modeldictionary = NewDataParser()
		print modeldictionary
		pickle.dump(modeldictionary,open("validation_save.p","wb"))
	

	sumdict = {}
	Rank(modeldictionary,simple)


	PlotHeatmap(modeldictionary,"All Dock")

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
