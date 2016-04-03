#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import operator
import Queue
import ast

# исходное вещество для первой итерации генетического алгоритма (adam)
# задается списком молекулярных групп
#
# quality_func по формуле smiles вычисляет на сколько полученное вещество близко
# к заданным пользователем параметрам (чем ближе к 0, тем лучше)
# если полученное вещество вообще недопустимо, функция должна вернуть None
#
# возвращаемое значение -- список сгенерированных веществ, отсортированный по
# убыванию результатов функции оценки
def GenAlg(quality_func, fragment_list, adam, max_iters = 1000, pop_length = 300):
    # словарь полученных веществ
    substances = { str(adam): quality_func(''.join(adam))}
    # очередь веществ на мутацию
    queue = Queue.Queue();
    queue.put(adam)
    # функци мутации
    def generate(mol_groups):
        # копирование исходного вещества
        newmol0 = mol_groups[:]
        # рассчет индексов и выполнение мутаций
        indm = np.random.randint(len(mol_groups))
        newmol0[indm] = fragment_list[np.random.randint(len(fragment_list))]
        start_ind1 = np.random.randint(len(mol_groups))
        start_ind2 = np.random.randint(len(mol_groups))
        finish_ind1 = start_ind1 + np.random.randint(len(mol_groups)-start_ind1)
        finish_ind2 = start_ind2 + np.random.randint(len(mol_groups)-start_ind2)
        newmol1 = mol_groups[:][:start_ind1] + newmol0[:][start_ind2:finish_ind2] + mol_groups[:][finish_ind1:]
        newmol2 = newmol0[:][:start_ind2] + mol_groups[:][start_ind1:finish_ind1] + mol_groups[:][finish_ind2:]
        # список новых веществ
        return [newmol0, newmol1, newmol2]
    # рассчет популяций
    for generation in range(max_iters):
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
