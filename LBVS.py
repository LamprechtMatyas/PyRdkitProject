import rdkit
import rdkit.Chem
import gzip
import random
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Scoring import Scoring
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import AllChem


# from molecules_arr it picks elements with indices in numbers_arr
def set_choosing(molecules_arr, numbers_arr):
    arr = []
    for i in range(len(numbers_arr)):
        arr.append(molecules_arr[numbers_arr[i]])
    return arr


# selects active molecules from train set, that means this molecules that have in array active value [1]
def select_active_molecules_from_train(mol, active):
    active_ones = []
    for i in range(len(mol)):
        if active[i] == [1]:
            active_ones.append(mol[i])
    return active_ones


# computing the similarity of active_molecules1 and test_molecules by using topological fingerprints
def topological_fingerprints(active_molecules1, test_molecules):
    similarity = []
    active_molecules_fps = [FingerprintMols.FingerprintMol(p) for p in active_molecules1]
    test_molecules_fps = [FingerprintMols.FingerprintMol(p) for p in test_molecules]
    for i in range(len(test_molecules_fps)):
        num_sim = 0
        for j in range(len(active_molecules_fps)):
            sim = DataStructs.FingerprintSimilarity(test_molecules_fps[i], active_molecules_fps[j])
            if sim > num_sim:
                num_sim = sim
        similarity.append(num_sim)
    return similarity


# computing the similarity of active_molecules1 and test_molecules by using atom pairs similarity
def atom_pairs_similarity(active_molecules1, test_molecules):
    similarity = []
    active_molecules_pairfps = [Pairs.GetAtomPairFingerprint(p) for p in active_molecules1]
    test_molecules_pairsfps = [Pairs.GetAtomPairFingerprint(p) for p in test_molecules]
    for i in range(len(test_molecules_pairsfps)):
        num_sim = 0
        for j in range(len(active_molecules_pairfps)):
            sim = DataStructs.DiceSimilarity(test_molecules_pairsfps[i], active_molecules_pairfps[j])
            if sim > num_sim:
                num_sim = sim
        similarity.append(num_sim)
    return similarity


# computing the similarity of active_molecules1 and test_molecules by using extended connectivity fingerprints
def ecfp_similarity(active_molecules1, test_molecules):
    similarity = []
    active_molecules_ecfpfps = [AllChem.GetMorganFingerprint(p, 3) for p in active_molecules1]
    test_molecules_ecfpfps = [AllChem.GetMorganFingerprint(p, 3) for p in test_molecules]
    for i in range(len(test_molecules_ecfpfps)):
        num_sim = 0
        for j in range(len(active_molecules_ecfpfps)):
            sim = DataStructs.DiceSimilarity(test_molecules_ecfpfps[i], active_molecules_ecfpfps[j])
            if sim > num_sim:
                num_sim = sim
        similarity.append(num_sim)
    return similarity


# it sorts the molecules similarity from highest to lowest and returns their activity in this ordering
def similarity_sorting(molecules_sim, activity):
    for i in range(len(molecules_sim)-1):
        for j in range(len(molecules_sim)-i-1):
            if molecules_sim[j+1] > molecules_sim[j]:
                tmp = molecules_sim[j+1]
                molecules_sim[j+1] = molecules_sim[j]
                molecules_sim[j] = tmp
                tmp = activity[j+1]
                activity[j+1] = activity[j]
                activity[j] = tmp
    return activity


# evaluation of the model
def evaluate(activity_arr):
    auc = Scoring.CalcAUC(activity_arr, 0)
    print("AUC: ", auc)
    ef = Scoring.CalcEnrichment(activity_arr, 0, [0.01])
    print("EF for 1%: ", ef[0])
    ef = Scoring.CalcEnrichment(activity_arr, 0, [0.05])
    print("EF for 5%: ", ef[0])
    rie = Scoring.CalcRIE(activity_arr, 0, 100)
    print("RIE for 100: ", rie)
    bedroc = Scoring.CalcBEDROC(activity_arr, 0, 100)
    print("BEDROC for 100: ", bedroc)


# importing molecules from file
inf = gzip.open('data/actives_5ht3.sdf.gz')
gzsuppl = rdkit.Chem.ForwardSDMolSupplier(inf)
ms = [x for x in gzsuppl if x is not None]       # ms - array where are molecules

# lets say that first 20 molecules are active
number_of_molecules = len(ms)
number_of_active_molecules = 20
number_of_inactive_molecules = number_of_molecules - number_of_active_molecules
act = [[1] for x in range(0, number_of_active_molecules)]
inact = [[0] for x in range(0, number_of_inactive_molecules)]
score = act + inact                            # first 20 elements are [1], the rest are [0]

# random splitting into train and test set 50:50
num_train = number_of_molecules // 2
num_test = number_of_molecules - num_train
random.seed(123)              # for replicating our results
s = random.sample(range(number_of_molecules), k=number_of_molecules)
indices_train = s[0:num_train]
indices_test = s[num_train:number_of_molecules]
train = set_choosing(ms, indices_train)
test = set_choosing(ms, indices_test)
train_score = set_choosing(score, indices_train)
test_score = set_choosing(score, indices_test)

# selecting active molecules from training set
active_molecules = select_active_molecules_from_train(train, train_score)
test_score1 = test_score.copy()

# topological fingerprints
topolog_fgp = topological_fingerprints(active_molecules, test)
sorted_act = similarity_sorting(topolog_fgp, test_score)
print("Topological fingerprints")
evaluate(sorted_act)
print()

# atom pairs
test_score = test_score1.copy()
atom_pairs = atom_pairs_similarity(active_molecules, test)
sorted_act = similarity_sorting(atom_pairs, test_score)
print("Atom pairs: ")
evaluate(sorted_act)
print()

# extended connectivity fingerprints
test_score = test_score1.copy()
ecfp = ecfp_similarity(active_molecules, test)
sorted_act = similarity_sorting(ecfp, test_score)
print("Extended connectivity fingerprints: ")
evaluate(sorted_act)
print()







