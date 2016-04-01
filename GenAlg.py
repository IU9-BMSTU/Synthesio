#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from LinearRegressionModel import CreateModel, ReadData, descs
import Queue
import operator

glb = { 'C(=O)O': 2.1, 
        'O'     : 1.9, 
        'C'     : -0.5, 
        'N'     : 9.4, 
        'OS(=O)(=O)[O-][Na+]'   : 38.7, 
        'C(=O)[O-][Na+]'        : 19.1, 
        'C(=O)[O-][K+]'         : 21.1 }
def count_dglb(mol_gr):
        return sum([glb[el] for el in mol_gr])
    

def GenAlg(quality_func, fragment_list, adam, max_iters = 1000, pop_length = 300):
    smi_adam = ''.join(adam)
    all_mols = { smi_adam: quality_func(smi_adam) }
    mol_molgr = { smi_adam: adam }
    
    
    mol_q = Queue.Queue()
    mol_q.put(adam)
    
    def generate(mol_groups):
        newmol0 = mol_groups[:]
        
        indm = np.random.randint(len(mol_groups))
        newmol0[indm] = fragment_list[np.random.randint(len(fragment_list))]
        
        start_ind1 = np.random.randint(len(mol_groups))
        start_ind2 = np.random.randint(len(mol_groups))
        finish_ind1 = start_ind1 + np.random.randint(len(mol_groups)-start_ind1)
        finish_ind2 = start_ind2 + np.random.randint(len(mol_groups)-start_ind2)
        
        newmol1 = mol_groups[:][:start_ind1] + newmol0[:][start_ind2:finish_ind2] + mol_groups[:][finish_ind1:]
        newmol2 = newmol0[:][:start_ind2] + mol_groups[:][start_ind1:finish_ind1] + mol_groups[:][finish_ind2:]
        
        return [newmol0, newmol1, newmol2]
    
    for generation in range(max_iters):
        next_mols = generate(mol_q.get())
        #print next_mols
        for mol in next_mols:
            smi_mol = ''.join(mol)
            res = quality_func(smi_mol)
            if res:
                all_mols[smi_mol] = res
                mol_molgr[smi_mol] = mol
                mol_q.put(mol)
        
        if len(all_mols) > pop_length:
            sorted_mols = sorted(all_mols.items(), key =operator.itemgetter(1))#, reverse = True)
            all_mols = { sorted_mols[i][0] : sorted_mols[i][1] for i in range(len(all_mols)/4)}
            mol_q = Queue.Queue()
            for mol in all_mols:
                mol_q.put(mol_molgr[mol])
                
    #return sorted(all_mols.items(), key =operator.itemgetter(1), reverse = True)
    return [(mol, pr, count_dglb(mol_molgr[mol])) for (mol, pr) in sorted(all_mols.items(), key =operator.itemgetter(1))]
  
def create_quality_func(expected_v, predict_model, rules):
    def q(smiles):
        for rule in rules:
            if not rule(smiles):
                return None
        predicted = np.array(predict_model(smiles))
        return np.linalg.norm(predicted- np.array(expected_v))
    return q

def check_danger(smiles):
    if 'ON' in smiles or 'NO' in smiles or 'OO' in smiles:
        return False
    return True

print GenAlg(create_quality_func([13], CreateModel(ReadData('learn.csv'), descs), [check_danger]),
             ['C(=O)O', 'O', 'OS(=O)(=O)[O-][Na+]', 'C(=O)[O-][Na+]', 'C(=O)[O-][K+]', 'C', 'O', 'N'],
             ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'C', 'C', 'OS(=O)(=O)[O-][Na+]'],
             max_iters=1000)[:10]
        
                