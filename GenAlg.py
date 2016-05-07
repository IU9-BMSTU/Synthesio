#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import operator
import Queue
import ast

def GenAlg(quality_func, fragment_lists, adam, iters = 1000, pop_length = 300):
    fr_head, fr_tail, fr_end = fragment_lists[0], fragment_lists[1], fragment_lists[2]
    # словарь полученных веществ
    substances = { str(adam): quality_func(''.join(adam))}
    # очередь веществ на мутацию
    queue = Queue.Queue();
    queue.put(adam)
    # функция мутации
    def generate(mol_groups):
        # копирование исходного вещества
        newmol0 = mol_groups[:]
        # рассчет индексов и выполнение мутаций
        indm = np.random.randint(len(mol_groups))
        if indm == len(mol_groups) - 1:
            newmol0[indm] = fr_end[np.random.randint(len(fr_end))]
        elif indm == 0:
            newmol0[indm] = fr_head[np.random.randint(len(fr_head))]
        else:
            newmol0[indm] = fr_tail[np.random.randint(len(fr_tail))]
        start_ind1 = np.random.randint(len(mol_groups)-2) + 1
        start_ind2 = np.random.randint(len(mol_groups)-2) + 1
        finish_ind1 = start_ind1 + np.random.randint(len(mol_groups)-start_ind1-1)
        finish_ind2 = start_ind2 + np.random.randint(len(mol_groups)-start_ind2-1)
        newmol1 = mol_groups[:][:start_ind1] + newmol0[:][start_ind2:finish_ind2] + mol_groups[:][finish_ind1:]
        newmol2 = newmol0[:][:start_ind2] + mol_groups[:][start_ind1:finish_ind1] + mol_groups[:][finish_ind2:]
        # список новых веществ
        return [newmol0, newmol1, newmol2]
    # рассчет популяций
    for generation in range(iters):
        if queue.empty():
            print 'bad start parameters'
            print 'population was extinct'
            break
        next_mols = generate(queue.get())
        for mol in next_mols:
            res = quality_func(''.join(mol))
            if res:
                substances[str(mol)] = res
                queue.put(mol)
        # проверка размера популяций
        if len(substances) > pop_length:
            # отбор лучших
            sorted_mols = sorted(substances.items(), key = operator.itemgetter(1))
            substances = { sorted_mols[i][0] : sorted_mols[i][1] for i in range(len(sorted_mols)/4) }
            queue = Queue.Queue()
            for mol in substances.keys():
                queue.put(ast.literal_eval(mol))
    return [(ast.literal_eval(mol), predict) for (mol, predict) in sorted(substances.items(), key =operator.itemgetter(1))]
