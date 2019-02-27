""""
Evaluation of a given model.

Usage:
    python evaluation.py
        -j (input json file or directory with json file that is an output of extract_fragments.py)
        -c (input csv file or directory with csv file that is an output of descriptors.py)
        -e (path to output file)

"""

import rdkit
import csv
import os
import argparse
import json
import random
from rdkit.ML.Scoring import Scoring


# Get and return application settings
def _read_configuration() -> dict:
    parser = argparse.ArgumentParser(
        description="Compute RDKit descriptors for given"
                    "molecules/fragments.")
    parser.add_argument("-j", type=str, dest="input_fragments",
                        help="input json file with fragments",
                        required=True)
    parser.add_argument("-c", type=str, dest="input_descriptors",
                        help="input csv with descriptors", required=True)
    parser.add_argument("-m", type=str, dest="model",
                        help="input model specification", required=False)
    parser.add_argument("-e", type=str, dest="evaluation",
                        help="file with evaluation values", required=True)
    parser.add_argument("-r", type=str, dest="results",
                        help="file with results", required=False)

    return vars(parser.parse_args())


# Create directory if it does not exists.
def create_parent_directory(path: str):
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)


# read inputted csv file
def read_csv_file(file: str) -> list:
    with open(file, "r") as csv_reader:
        reader = csv.reader(csv_reader, delimiter=",")
        descriptors = list(reader)
    return descriptors


# read inputted json file and return array with molecules their activities and the positions of their fragments
def read_json_file(file: str) -> list:
    start = 1
    molecule_smiles = []
    with open(file, "r") as json_reader:
        line_number = 0
        for line in json_reader:
            molecule = json.loads(line)
            molecule_smiles.append([molecule["smiles"]])
            molecule_smiles[line_number].append(int(molecule["active"]))
            molecule_smiles[line_number].append(start)
            count = len(molecule["fragments"])
            start = start + count
            molecule_smiles[line_number].append(start-1)
            line_number += 1
    return molecule_smiles


# return rows from molecules_arr with row number in numbers_arr
def choose_sets(molecules_arr, numbers_arr) -> list:
    arr = []
    for i in range(len(numbers_arr)):
        arr.append(molecules_arr[numbers_arr[i]])
    return arr


# select active molecules from the train set
def select_active_molecules_from_train(molecules: list) -> list:
    active_ones = []
    for i in range(len(molecules)):
        if molecules[i][1] == 1:
            active_ones.append(molecules[i])
    return active_ones


# return fragments from active molecules from training set
def select_active_fragments(active_molecules: list, fragments: list) -> list:
    active_fragments = []
    for molecule in active_molecules:
        for i in range(molecule[2], molecule[3]+1):
            active_fragments.append(fragments[i])
    return active_fragments


# for every molecule computes the ratio of active fragment to total fragments
def compute_similarity(active_molecules: list, test_molecules: list, descriptors: list) -> list:
    similarity = []
    active_fragments = select_active_fragments(active_molecules, descriptors)
    for molecule in test_molecules:
        total = 0
        sim = 0
        for fragment in range(molecule[2], molecule[3]+1):
            for active_fragment in active_fragments:
                if active_fragment == descriptors[fragment]:
                    sim += 1
                    break
            total += 1
        similarity.append(sim / total)
    return similarity


# return array of activity(0,1) sorted by the similarity
def sort_by_similarity(similarity: list, molecules: list) -> list:
    sim_active = []
    for i in range(len(similarity)):
        sim_active.append([similarity[i]])
        sim_active[i].append(molecules[i][1])
    sim_active = sorted(sim_active, reverse=True)
    active = []
    for i in range(len(sim_active)):
        active.append([sim_active[i][1]])
    return active


# evaluation of model
def evaluation(activity_arr: list, output_file: str):
    with open(output_file, "w") as stream:
        auc = Scoring.CalcAUC(activity_arr, 0)
        stream.write("AUC: ")
        stream.write(str(auc))
        ef = Scoring.CalcEnrichment(activity_arr, 0, [0.01])
        stream.write("\nEF for 1%: ")
        stream.write(str(ef[0]))
        ef = Scoring.CalcEnrichment(activity_arr, 0, [0.05])
        stream.write("\nEF for 5%: ")
        stream.write(str(ef[0]))
        rie = Scoring.CalcRIE(activity_arr, 0, 100)
        stream.write("\nRIE for 100: ")
        stream.write(str(rie))
        bedroc = Scoring.CalcBEDROC(activity_arr, 0, 100)
        stream.write("\nBEDROC for 100: ")
        stream.write(str(bedroc))


# divide data into train and test set and sorts the test ones by similarity
def test_model_and_compute_similarity(molecules_smiles: list, descriptors: list):
    num_of_molecules = len(molecules_smiles)
    num_train = num_of_molecules // 2
    random.seed(12)
    s = random.sample(range(num_of_molecules), k=num_of_molecules)
    indices_train = s[0:num_train]
    indices_test = s[num_train:num_of_molecules]
    train = choose_sets(molecules_smiles, indices_train)
    test = choose_sets(molecules_smiles, indices_test)
    active_train = select_active_molecules_from_train(train)
    similarity = compute_similarity(active_train, test, descriptors)
    sorted_sim = sort_by_similarity(similarity, test)
    return sorted_sim


def _main():
    configuration = _read_configuration()
    descriptors = read_csv_file(configuration["input_descriptors"])
    molecules_smiles = read_json_file(configuration["input_fragments"])
    sorted_sim = test_model_and_compute_similarity(molecules_smiles, descriptors)
    create_parent_directory(configuration["evaluation"])
    evaluation(sorted_sim, configuration["evaluation"])


if __name__ == "__main__":
    _main()
