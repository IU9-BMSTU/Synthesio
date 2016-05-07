#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '../')

from LinearRegressionModel import CreateModel, ReadData
from GenAlg import GenAlg

import pybel
import numpy as np

from prettytable import PrettyTable

import time

class Profiler(object):
    def __enter__(self):
        self._startTime = time.time()
    def __exit__(self, type, value, traceback):
        print "Время выполнения: {:.3f} с".format(time.time() - self._startTime)

def create_quality_func(expected_v, predict_model, rules):
    def q(smiles):
        for rule in rules:
            if not rule(smiles):
                return None
        predicted = np.array(predict_model(smiles))
        return np.linalg.norm(predicted- np.array(expected_v))
    return q

def expected_mass(low, high):
    return (lambda smiles: True if low <= pybel.readstring("smi", smiles).molwt <= high else False)

def bad_seq():
    return (lambda smiles: False if "OO" in smiles else True)

descs = [lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MW'])['MW'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['bonds'])['bonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['atoms'])['atoms'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['sbonds'])['sbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['dbonds'])['dbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MR'])['MR']]

def print_results(results):
    t = PrettyTable(['smiles', 'predicted quality', 'logP', 'molwt'])
    [t.add_row([''.join(mol), predict,
                pybel.readstring('smi', ''.join(mol)).calcdesc(['logP'])['logP'],
                pybel.readstring("smi", ''.join(mol)).molwt])
            for mol, predict in results[:5]]
    print t

#===============================================================================
fragments_head = ["N", "C", "Cl", "Br", "I"]
fragments_end = ["CO", "C(C=O)", "C", "Cl", "Br", "I", "O", "C"]
fragments = ["C", "S", "N", "(O=)S(=O)", "C(=O)O", "C(C=O)", "CCOCC", "CC"]
#===============================================================================
print 'TEST #1-A'
print 'Целевое значение logP: -1,3'
print 'Ограничения: Мол. масса в пределах 50 - 70 г/моль'
print 'Количество итераций генетического алгоритма: 500'
with Profiler() as stopwatch:
    results = GenAlg(create_quality_func([-1.3],
                                         CreateModel(ReadData('logP-learn.csv'), descs),
                                         [expected_mass(50, 70)]),
                     [fragments_head, fragments, fragments_end],
                     ['O', 'C', 'CO'],
                     iters=500)
print_results(results)
#===============================================================================
fragments_head = ["C(C=O)", "C", "Cl", "Br", "I"]
fragments_end = ["CO", "C", "Cl", "Br", "I", "O"]
fragments = ["C", "S", "(O=)S(=O)", "C(C=O)", "CCOCC", "CC"]
#===============================================================================
print 'TEST #1-B'
print 'Целевое значение logP: -0,24'
print 'Ограничения: Мол. масса в пределах 30 - 80 г/моль'
print 'Количество итераций генетического алгоритма: 800'
with Profiler() as stopwatch:
    results = GenAlg(create_quality_func([-0.24],
                                         CreateModel(ReadData('logP-learn.csv'), descs),
                                         [expected_mass(30, 80), bad_seq()]),
                    [fragments_head, fragments, fragments_end],
                    ['C', 'C', 'O', 'CO'],
                    iters=800)
print_results(results)
#===============================================================================
fragments_end = fragments_head = ["C", "Cl"]
fragments = ["C"]
#===============================================================================
print 'TEST #1-C'
print 'Целевое значение logP: +0,48'
print 'Ограничения: Мол. масса в пределах 80 - 110 г/моль'
print 'Количество итераций генетического алгоритма: 1000'
with Profiler() as stopwatch:
    results = GenAlg(create_quality_func([0.48],
                                         CreateModel(ReadData('logP-learn.csv'), descs),
                                         [expected_mass(80, 110)]),
                    [fragments_head, fragments, fragments_end],
                    ['Cl', 'C', 'C', 'CO'],
                    iters=1000)
print_results(results)
'''
#===============================================================================

print 'Test #2-Лекарственное вещество с заданным logP'
print 'Целевое значение logP: +0,75'
print 'Ограничения: Мол. масса в пределах 220 - 290 г/моль'

Обучаться можно на этом:
NC1=CC=C(C=C1)S(N)(=O)=O
CC1=NOC(NS(=O)(=O)C2=CC=C(N)C=C2)=C1C
CC1=NN=C(NS(=O)(=O)C2=CC=C(N)C=C2)S1
CC(=O)NS(=O)(=O)C1=CC=C(N)C=C1
CCN1C=CC(NS(=O)(=O)C2=CC=C(N)C=C2)=NC1=O
CC1=NC(NS(=O)(=O)C2=CC=C(N)C=C2)=NC=C1
CC1=CC(C)=NC(NS(=O)(=O)C2=CC=C(N)C=C2)=N1
NC1=CC=C(C=C1)S(=O)(=O)NC1=CC=NN1C1=CC=CC=C1
В выдаче должно присутствовать:
CC1=CC(NS(=O)(=O)C2=CC=C(N)C=C2)=NO1
'''
