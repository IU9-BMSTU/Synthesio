#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '../')

from LinearRegressionModel import CreateModel, ReadData
from GenAlg import GenAlg

import pybel
import numpy as np

from prettytable import PrettyTable

def create_quality_func(expected_v, predict_model, rules):
    def q(smiles):
        for rule in rules:
            if not rule(smiles):
                return None
        predicted = np.array(predict_model(smiles))
        return np.linalg.norm(predicted- np.array(expected_v))
    return q

descs = [lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MW'])['MW'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['bonds'])['bonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['atoms'])['atoms'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['sbonds'])['sbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['dbonds'])['dbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MR'])['MR']]


fragments_head = ["C(=O)", "C(=O)O", "O", "C", "Cl", "Br", "I", "F", "N", "C#N", "N=O"]
fragments_end = fragments_head
fragments = ["C", "=C", "C=", "(O=)S(=O)", "S", "O(O=)S(=O)", "N=N", "N=C", "C(=O)O", "C(=O)", "O", "C(=O)N", "C#C"]
results = GenAlg(create_quality_func([5], CreateModel(ReadData('logP-learn.csv'), descs), []),
                [fragments_head, fragments, fragments_end],
                ['F', 'N=N', 'C', 'C', 'N'],
                max_iters=500)

# top-5
t = PrettyTable(['smiles', 'predicted quality', 'real GLB index'])
[t.add_row([''.join(mol), predict, pybel.readstring('smi', ''.join(mol)).calcdesc(['logP'])['logP']]) for mol, predict in results[:5]]
print t
