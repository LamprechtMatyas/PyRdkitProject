import rdkit
import rdkit.Chem
from rdkit.Chem import FragmentCatalog
from rdkit import RDConfig
import os


# vypsani fragmentu z FunctionalGroups.txt
fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
fparams = FragmentCatalog.FragCatParams(1, 3, fName)
print()
x = fparams.GetNumFuncGroups()
print(x)
print("fragmenty y FunctionalGroups.txt:")
for i in range(0, x):
    print(fparams.GetFuncGroup(i).GetProp('_Name'))

print()


# nacteni a vypsani fragmentu jedne molekuly
m = rdkit.Chem.MolFromSmiles('n1ccccc1')
print('Nyni fragmenty pro molekulu: ' + rdkit.Chem.MolToSmiles(m))
fcat = FragmentCatalog.FragCatalog(fparams)
fcgen = FragmentCatalog.FragCatGenerator()
pocetFragmentu = fcgen.AddFragsFromMol(m, fcat)
for i in range(0, pocetFragmentu):
    print(fcat.GetEntryDescription(i))


m1 = rdkit.Chem.MolFromSmiles('CCC')
fcat = FragmentCatalog.FragCatalog(fparams)
fcgen = FragmentCatalog.FragCatGenerator()
pocetFragmentu = fcgen.AddFragsFromMol(m1, fcat)
print()
print('Pro molekulu: ' + rdkit.Chem.MolToSmiles(m1))
for i in range(0, pocetFragmentu):
    print(fcat.GetEntryDescription(i))

print()


# nacteni a vypsani fragmentu pro vice molekul
print('Vice molekul:')
suppl = rdkit.Chem.SDMolSupplier('data/5ht3ligs.sdf')
for mol in suppl:
    print(rdkit.Chem.MolToSmiles(mol))

print()
for mol in suppl:
    fcat = FragmentCatalog.FragCatalog(fparams)
    fcgen = FragmentCatalog.FragCatGenerator()
    print(rdkit.Chem.MolToSmiles(mol))
    pocetFragmentu = fcgen.AddFragsFromMol(mol, fcat)
    for i in range(0, pocetFragmentu):
        print(fcat.GetEntryDescription(i))
    print()





