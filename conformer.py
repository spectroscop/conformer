#import time
#import pickle
import os
import datetime
import argparse

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize
rdDepictor.SetPreferCoordGen(True)
#import rdkit
#print(rdkit.__version__)
#2020.09.4

class Settings:
    InputFiles   = ''    # Structure input files - *.mol files
    OutputFolder = ''    # folder to print dp4 output to - default is cwd
    OutputFiles  = ''
    RMSDcutoff   = 0.3   # Initial RMSD threshold for pruning or 0.1 or 0.74
    ConfLimit    = 99    # Max numbers of conformers allowed per structure for DFT stages
    add_ref      = True
    tautomer     = False # generate tautomer
    macrocycle   = False # molecule into macrocycle
    verbose      = False # verbose for check generator conformer
    ff           = 'MMFF94s' # verbose for check generator conformer
    numThreads   = 0

class Isomer:
    def __init__(self):
        self.Molecules = []             # Molecules
        self.Conformers = []            # from conformational search, list of atom coordinate lists
        self.MMEnergies = []            # Corresponding MM energies in kj/mol

#=========================

settings = Settings()

#=========================

def reorderTautomers(m):
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(m)
    tauts = enumerator.Enumerate(m)
    smi_canon = Chem.MolToSmiles(canon)
    smi_tauts = [Chem.MolToSmiles(x) for x in tauts]
    stpl = sorted((x,y) for x,y in zip(smi_tauts,tauts) if x!=smi_canon)
    res = [canon]
    res += [y for x,y in stpl]
    
    return res

#==============

def print_failure_causes(counts):
    for i,k in enumerate(rdDistGeom.EmbedFailureCauses.names):
        print(k,counts[i])
    # in v2022.03.1 two names are missing from `rdDistGeom.EmbedFailureCauses`:
    print('LINEAR_DOUBLE_BOND',counts[i+1])
    print('BAD_DOUBLE_BOND_STEREO',counts[i+2])    
    '''
INITIAL_COORDS:            generation of the initial coordinates from the random distance matrix (default) or from a set of random coordinates (when using random coordinate embedding) failed.
FIRST_MINIMIZATION:        the initial optimization of the atom positions using the distance-geometry force field failed to produce a low-enough energy conformer. The check here has thresholds for both average energy per atom and the individual atom energies. I’m not providing the threshold values here since the energies from the distance-geometry force field are not physically meaningful - the threshold values are not interpretable.
CHECK_TETRAHEDRAL_CENTERS: at least one tetrahedral C and N centers either has a volume around it which is too small or is outside the volume defined by its neighbors
MINIMIZE_FOURTH_DIMENSION: the minmization to force the values of the fourth-dimensional component of each atom position failed
ETK_MINIMIZATION:          after the minimization with the ET and/or K terms, at least one atom which should have been planar was not
FINAL_CHIRAL_BOUNDS:       the neighborhood of an atom with specified chirality was too distorted (it violated distance constraints)
FINAL_CENTER_IN_VOLUME:    an atom with specified chirality was outside of the volume defined by its neighbors
LINEAR_DOUBLE_BOND:        one of the end atoms of a double bond had a linear geometry
BAD_DOUBLE_BOND_STEREO:    the stereochemistry of a double bond with specified stereochemistry was wrong in the generated conformer
    '''
    return None

#================

def MMFFOptimizeMolecule (settings, mol):

    if settings.macrocycle:
        param_etkdg = rdDistGeom.ETKDGv3 ()
    else:
        param_etkdg = rdDistGeom.srETKDGv3 ()
    param_etkdg.randomSeed           = 0xf00d
    param_etkdg.pruneRmsThresh       = settings.RMSDcutoff
    param_etkdg.onlyHeavyAtomsForRMS = True
    param_etkdg.numThreads           = settings.numThreads
    param_etkdg.trackFailures        = True
    param_etkdg.enforceChirality     = True

#    metric_etkd = []
#    t1 = time.time()
    mol = Chem.AddHs (mol, addCoords=True)
    conf = rdDistGeom.EmbedMultipleConfs (mol, settings.ConfLimit, param_etkdg)
#    t2 = time.time()
    if settings.verbose:
        counts = param_etkdg.GetFailureCounts ()
        print_failure_causes (counts)
    if conf==-1:
        print ("Structure molecule error")
        return None
#    metric_etkdg_res.append((t2-t1,mol.GetNumConformers(),mol,mol.GetIntProp('RTB')))
#    with open('./results/random_coords_expt.pkl','wb+') as outf:
#       pickle.dump((metric_etkd,),outf)

    rdMolTransforms.CanonicalizeMol(mol, normalizeCovar=True )
    #AlignMolConformers used to calculate the RMSD
    rmslist = []
    rms = []
    rms = rdMolAlign.AlignMolConformers (mol, RMSlist=rmslist)
 
    if settings.ff == 'MMFF94s':
        AllChem.MMFFSanitizeMolecule (mol)
        mmff_props = AllChem.MMFFGetMoleculeProperties (mol, mmffVariant='MMFF94s')
        #AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=settings.numThreads, maxIters=1000, mmffVariant='MMFF94s')
    elif  settings.ff == 'MMFF94':
        AllChem.MMFFSanitizeMolecule (mol)
        mmff_props = AllChem.MMFFGetMoleculeProperties (mol, mmffVariant='MMFF94')
        #AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=settings.numThreads, maxIters=1000, mmffVariant='MMFF94')
    '''
   The UFF is an all atom force field, so it considers all the atoms in the 
periodic table. Having transition metal, for example, means to use UFF.
   MMFF94 considers the atoms: C, H, N, O, F, Si, P, S, Cl, Br, and ions: 
Fe+2, Fe+3, F-, Cl-, Br-, Li+, Na+, K+, Zn+2, Ca+2, Cu+1, Cu+2, and Mg+2.
   Furthemore MMFF94 has two versions, the MMFF94 and MMFF94(s). If you 
consider crystal structures, you should use MMFF94s (s=static) while for 
“docking” simulations, use the MMFF94 variant. The difference between the 
two is due to a high pyramidal puckering on the N atoms in the MMFF94 when 
optimizing crystal derived structures. Using the MMFF94s avoids such a problem.
    '''

    if settings.add_ref:
        res = []
        for cid in conf:
#            if settings.ff == 'MMFF94s':
#                ff = AllChem.MMFFGetMoleculeForceField (mol, mmff_props, confId=cid)
#            elif settings.ff == 'UFF':
#                ff = AllChem.UFFGetMoleculeForceField (mol, confId=cid)
            if settings.ff == 'UFF':
                ff = AllChem.UFFGetMoleculeForceField (mol, confId=cid)
            else:
                ff = AllChem.MMFFGetMoleculeForceField (mol, mmff_props, confId=cid)
            ff.Initialize ()
            ff.Minimize (maxIts=1000)
            e = ff.CalcEnergy ()
            res.append ((cid, e))
        
        sorted_res = sorted (res, key=lambda x:x[1])
        #AlignMolConformers used to calculate the RMSD
        rmslist = []
        rms = rdMolAlign.AlignMolConformers(mol, RMSlist=rmslist)
        print ('Generate ' + str (mol.GetNumConformers()) + ' conformers.' )

    return mol, conf, sorted_res

#==============

def writesdf (settings, Isomers):
    w = Chem.SDWriter (settings.OutputFiles)
    for i_ in range(0, len(Isomers.Molecules)):
        if settings.tautomer:
            Isomers.Molecules[i_].SetProp ('tautomers_id', str (i_))
        if settings.add_ref:
            for cid, e in Isomers.MMEnergies[i_]:
                Isomers.Molecules[i_].SetProp ('CID',    str (cid))
                Isomers.Molecules[i_].SetProp ('Energy', str (e))
                w.write (Isomers.Molecules[i_], confId=cid)
        else:
            Isomers.Molecules[i_].SetProp ('CID', '-1')
            Isomers.Molecules[i_].SetProp ('Energy', '')
            w.write (Isomers.Molecules[i_])
    w.close ()

#==============

def confgen (settings):

    mol   = Chem.MolFromMolFile(settings.InputFiles)
    if mol is None:
        print (settings.InputFiles + "- not file MOL")
        return None
    refmol = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    if refmol.GetNumHeavyAtoms()>50:
        print ("\n Molecule heavy more 50 atoms")
        return None

    # calculates the explicit and implicit valences on all atoms. This generates exceptions for atoms in higher-than-allowed valence states. This step is always performed, but if it is “skipped” the test for non-standard valences will not be carried out.
    refmol.UpdatePropertyCache(strict=False) # re-calculates the explicit and implicit valences on
    # next, you probably want to at least do a partial sanitization so that the molecule is actually useful:
    try:
         Chem.SanitizeMol(refmol,\
                            Chem.SanitizeFlags.SANITIZE_CLEANUP|\
                            Chem.SanitizeFlags.SANITIZE_SYMMRINGS|\
                            Chem.SanitizeFlags.SANITIZE_KEKULIZE|\
                            Chem.SanitizeFlags.SANITIZE_FINDRADICALS|\
                            Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|\
                            Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|\
                            Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|\
                            Chem.SanitizeFlags.SANITIZE_ADJUSTHS,\
                            catchErrors=True)
    except:
        print ("\n not clean molecule")
            
    '''
cleanUp: standardizes a small number of non-standard valence states. The clean up operations are:
    Neutral 5 valent Ns with double bonds to Os are converted to the zwitterionic form. Example: N(=O)=O -> [N+](=O)O-]
    Neutral 5 valent Ns with triple bonds to another N are converted to the zwitterionic form. Example: C-N=N#N -> C-N=[N+]=[N-]
    Neutral 5 valent phosphorus with one double bond to an O and another to either a C or a P are converted to the zwitterionic form. Example: C=P(=O)O -> C=[P+]([O-])O
    Neutral Cl, Br, or I with exclusively O neighbors, and a valence of 3, 5, or 7, are converted to the zwitterionic form. This covers things like chlorous acid, chloric acid, and perchloric acid. Example: O=Cl(=O)O -> [O-][Cl+2][O-]O
SANITIZE_SYMMRINGS: ???
Kekulize: converts aromatic rings to their Kekule form. Will raise an exception if a ring cannot be kekulized or if aromatic bonds are found outside of rings.
assignRadicals: determines the number of radical electrons (if any) on each atom.
setAromaticity: identifies the aromatic rings and ring systems (see above), sets the aromatic flag on atoms and bonds, sets bond orders to aromatic.
setConjugation: identifies which bonds are conjugated
setHybridization: calculates the hybridization state of each atom
adjustHs: adds explicit Hs where necessary to preserve the chemistry. This is typically needed for heteroatoms in aromatic rings. The classic example is the nitrogen atom in pyrrole.
    '''

    # if many fragments, get the "parent" (the actual mol we are interested in) 
    parent_mol = rdMolStandardize.FragmentParent(refmol)

    # try to neutralize molecule
    uncharger = rdMolStandardize.Uncharger() # annoying, but necessary as no convenience method exists
    uncharged_mol = uncharger.uncharge(parent_mol)

    # Create isomer data structures
    Isomers = Isomer()
    
    if settings.tautomer:
        # note that no attempt is made at reionization at this step
        # nor at ionization at some pH (rdkit has no pKa caculator)
        # the main aim to to represent all molecules from different sources
        # in a (single) standard way, for use in ML, catalogue, etc.
         
        tautomer = rdMolStandardize.TautomerEnumerator()    #init generator tautomer
        try:
            tautomers_cinf = tautomer.Canonicalize(uncharged_mol)               #canonical tautomer
        except:
            print('\n Tautomers not canonization')
            return None
        #tauts = tautomer.Enumerate(mol)                   #generate tautomer
        Isomers.Molecules = reorderTautomers(tautomers_cinf)                  #generate tautomer into sort canonic order
        print ('Generate ' + str (len(Isomers.Molecules)) + ' tautomers.' )
    else:
        Isomers.Molecules = [uncharged_mol]
    
    for i_ in range(0, len(Isomers.Molecules)):
        mol, conformers, MMEnergies = MMFFOptimizeMolecule (settings, Isomers.Molecules[i_])
        Isomers.Molecules[i_]=mol
        Isomers.Conformers.append(conformers)
        Isomers.MMEnergies.append(MMEnergies)
    
    writesdf (settings, Isomers)
    return None

#=======================

def run():

    # These are then overridden by any explicit parameters given through the command line
    parser = argparse.ArgumentParser(description='generate conformer')
    parser.add_argument('-i', '--InputFile', help="inputfile MOL")
    parser.add_argument('-o', '--OutputFolder', help="Directory for output conforme name-cinformer.sdf", default=settings.OutputFolder)
    parser.add_argument('-r', "--RMSDcutoff", help='Retain conformations, default - 0.3', type=float, default=settings.RMSDcutoff)
    parser.add_argument('-c', "--ConfLimit", help="Specify maximum number of conformers per structure, default - 99", type=int)
    parser.add_argument('-a', "--add_ref", help='Write calck energy force field, default - true', type=bool, default=True)
    parser.add_argument('-m', "--macrocycle", help='Use the macrocycle torsions from ETKDGv3, default - false', type=bool, default=False)
    parser.add_argument("-t", "--tautomer", help='Generate tautomer', type=bool, default=False)
    parser.add_argument("-f", "--forcefield", help='force field MMFF94s (default), MMFF94 or UFF for metaloorganic', default=settings.ff)
    parser.add_argument("-v", "--verbose", type=bool, default=False)

    args = parser.parse_args()
 
    if args.InputFile is not None:
        settings.InputFiles = os.path.abspath(os.path.expanduser(args.InputFile))
        if not os.path.isfile(settings.InputFiles):
            print("\n No input files were found please use -h for help with input options quitting...")
            quit()
        suffix = os.path.splitext(os.path.basename(settings.InputFiles))[1]
        if suffix.lower() != '.mol':
            print("\n Input files not MOL extended")
            quit()
        print(' Input files: ', settings.InputFiles)
    else:
        print ("Not input file")
        return None

    if args.OutputFolder is not None:
    #    settings.OutputFolder = Path(args.OutputFolder)
        settings.OutputFolder = os.path.abspath(os.path.expanduser(args.OutputFolder))
        if not os.path.isdir(settings.OutputFolder):
            print("\n No directory for output were found please use -h for help with input options quitting...")
            quit()
    else:
        settings.OutputFolder = os.path.dirname(settings.InputFiles)

    print(' Output folder: ', settings.OutputFolder)
    settings.OutputFiles = settings.OutputFolder + '/' + os.path.splitext(os.path.basename(settings.InputFiles))[0] + '.sdf'
    if os.path.isfile(settings.OutputFiles):
        print("\n Output file already exists")
        quit()

    settings.RMSDcutoff = args.RMSDcutoff
    if args.ConfLimit:
        if args.ConfLimit > 99:
            settings.ConfLimit = 99
        elif args.ConfLimit < 0:
            settings.ConfLimit = 99
        else:
            settings.ConfLimit = args.ConfLimit
    settings.add_ref    = args.add_ref
    settings.macrocycle = args.macrocycle
    settings.tautomer   = args.tautomer
    if args.forcefield:
        if args.forcefield == 'MMFF94s':
            settings.ff = 'MMFF94s'
        elif args.forcefield == 'MMFF94':
            settings.ff = 'MMFF94s'
        elif args.forcefield == 'UFF':
            settings.ff = 'UFF'
        else:
            settings.ff = 'MMFF94s'
    settings.verbose    = args.verbose

    #-----------------
    jobdir = os.getcwd()
    print(" Current working directory: " + os.getcwd())
    print (datetime.datetime.now())

    os.chdir(settings.OutputFolder)

    confgen(settings)

    print (datetime.datetime.now())
    print(' Process completed successfully.')

    os.chdir(jobdir)
    
    return None

#=========================

if __name__=='__main__':
    run()
