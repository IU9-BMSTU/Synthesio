#!/usr/bin/python
# -*- coding: utf-8 -*-

import sknn
import pickle

import numpy as np

from sknn.mlp import Regressor, Layer
from sknn.platform import cpu64, threading
#from sknn.platform import gpu32

import pybel

import logging
logging.basicConfig(level = logging.INFO)

import datetime
import sys

def ReadData(filename):
    # формат входных данных: csv таблица с разделителем ';'
    # первый столбец содержит формулу в формате smiles
    # остальные столбцы содержат числовые значения признаков
    return [line.rstrip().split(';') for line in open(filename)]

def CalcPredicates(smiles, predicates):
    return [f(smiles) for f in predicates]

def CreateNetwork(data, predicates):
    # входная размерность
    dim_in = len(predicates)
    # выходная размерность
    dim_out = len(data[0]) - 1
    # конфигурация сети
    neural_network = Regressor(
        layers=[
            Layer("Rectifier", units=50),
            Layer("Linear")],
        learning_rate=0.001,
        n_iter=5000)
    # формирование обучающей выборки
    x_train = np.array([CalcPredicates(row[0], predicates) for row in data])
    y_train = np.array([apply(float, row[1:]) for row in data])
    # обучение
    logging.info('Start training')
    logging.info('\n'+str(x_train))
    logging.info('\n'+str(y_train))
    try:
        neural_network.fit(x_train, y_train)
    except KeyboardInterrupt:
        logging.info('User break')
        pass
    logging.info('Network created successfully')
    logging.info('score = '+str(neural_network.score(x_train, y_train)))
    # сохранение обученной сети
    pickle.dump(neural_network, open(datetime.datetime.now().isoformat()+'.pkl', 'wb'))
    return neural_network

def CreateRegressor(descs, nn=None, learn=None, saved=None):
    if nn == None:
        if saved == None:
            logging.info('Creating neural network')
            nn = CreateNetwork(ReadData(learn), descs)
        else:
            logging.info('Loading network from file '+saved)
            nn = pickle.load(open(saved, 'rb'))
    return (lambda smiles: nn.predict(np.array([CalcPredicates(smiles, descs)]))[0])

descs = [lambda smiles: pybel.readstring('smi', smiles).calcdesc(['logP'])['logP'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MW'])['MW'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['bonds'])['bonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['atoms'])['atoms'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['sbonds'])['sbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['dbonds'])['dbonds'],
         lambda smiles: pybel.readstring('smi', smiles).calcdesc(['MR'])['MR']]

predict = CreateRegressor(learn='learn.csv', descs=descs)
