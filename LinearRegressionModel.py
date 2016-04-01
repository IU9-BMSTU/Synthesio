#!/usr/bin/python
# -*- coding: utf-8 -*-

from sklearn import linear_model
import pybel

def ReadData(filename):
    # формат входных данных: csv таблица с разделителем ';'
    # первый столбец содержит формулу в формате smiles
    # остальные столбцы содержат числовые значения признаков
    return [line.rstrip().split(';') for line in open(filename)]

def CalcPredicates(smiles, predicates):
    return [f(smiles) for f in predicates]

def CreateModel(data, predicates):
    model= [] # вектор линейных моделей
    for i in range(len(data[0])-1):
        clf = linear_model.LinearRegression()
        clf.fit([CalcPredicates(row[0], predicates) for row in data], [float(row[i+1]) for row in data])
        model.append(clf)
    return (lambda smiles: [clf.predict([CalcPredicates(smiles, predicates)])[0] for clf in model])

def Test(learn, test):
    func = CreateModel(ReadData(learn), descs)
    for row in ReadData(test):
        print 'Smiles:', row[0]
        print 'Expected value:', apply(float, row[1:])
        print 'Predicted value:', func(row[0])
        print

descs = [lambda smiles: pybel.readstring('smi', smiles).calcdesc(['logP'])['logP'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MW'])['MW'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['bonds'])['bonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['atoms'])['atoms'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['sbonds'])['sbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['dbonds'])['dbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MR'])['MR']]
