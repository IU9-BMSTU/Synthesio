#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '../')

from LinearRegressionModel import CreateModel, ReadData, descs
from GenAlg import GenAlg

import numpy as np

from prettytable import PrettyTable

# рассчет индекса GLB
glb = {"OS(=O)(=O)[O-][Na+]": 38.7,
       "OS(=O)(=O)[O-][K+]": 38.7,
       "C(=O)[O-][Na+]": 21.1,
       "C(=O)[O-][K+]": 19.1,
       "C(=O)[O-]":  2.1,
       "[N+]": 9.4,
       "O":1.9,
       "C":0.5
}

def count_glb(mol_gr):
        return sum([glb[el] for el in mol_gr])
# функция оценки качества сгенерированного мыла
def create_quality_func(expected_v, predict_model, rules):
    def q(smiles):
        for rule in rules:
            if not rule(smiles):
                return None
        predicted = np.array(predict_model(smiles))
        return np.linalg.norm(predicted- np.array(expected_v))
    return q
# функция проверки безопасности полученного мыла
def check_danger(smiles):
    if 'ON' in smiles or 'NO' in smiles or 'OO' in smiles:
        return False
    return True
# запуск генетического алгоритма
fragments_head = ["C"]
fragments_end = ["OS(=O)(=O)[O-][Na+]", "OS(=O)(=O)[O-][K+]", "C(=O)[O-][Na+]", "C(=O)[O-][K+]", "C(=O)[O-]"]
fragments = ["O", "C", "[N+]", "C"]
results = GenAlg(create_quality_func([13], CreateModel(ReadData('soap-learn.csv'), descs), [check_danger]),
                [fragments_head, fragments, fragments_end],
                ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'OS(=O)(=O)[O-][Na+]'],
                max_iters=900)

# top-5
t = PrettyTable(['smiles', 'predicted quality', 'real GLB index'])
[t.add_row([''.join(mol), predict, count_glb(mol)]) for mol, predict in results[:5]]
print t
