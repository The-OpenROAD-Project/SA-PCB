import random, math
import numpy as np
import numba
from numba import jit

import copy, sys, time
import collections
from operator import itemgetter

from tqdm import tqdm
import itertools

import utils

import multiprocessing as mp
import subprocess

try:
    num_cpus = mp.cpu_count()
except NotImplementedError:
    num_cpus = 2

"""
for now, make directories associated with each idx manually.
"""

def f(idx):
    res = subprocess.check_output(["./sa", "-p","designs/bm1","-i", str(2500),"-j",str(25),"-t",str(0.025), "-x", str(idx)])
    return res.strip()[-2:]

def worker(args, output):
	output.put(f(*args))

def multistart():
    K = max(4,mp.cpu_count())
    idx = -1
    for k in tqdm(range(K),desc='multistart'):
        processes = []
        manager = mp.Manager()
        output = manager.Queue()
        for i in range(num_cpus):
            p = mp.Process(target=worker, args=((i,),output))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()
            results = [output.get() for p in processes]
            print(results)
            best_result = max(results,key=itemgetter(1)) # max result by cost
            best_cost = best_result[1]
            best_idx = best_result[0]
    print('best result: ')
    print("idx: " + str(best_idx))
    print("cost: " + str(best_cost))
    return cost

multistart()

def nmultistart():
    K = max(4,mp.cpu_count())
    idx = -1
    for k in tqdm(range(K),desc='multistart'):
        processes = []
        manager = mp.Manager()
        output = manager.Queue()
        for i in range(num_cpus):
            p = mp.Process(target=worker, args=((i,),output))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()
            results = [output.get() for p in processes]
            print(results)
            #best_result = max(results,key=itemgetter(1)) # max result by cost
            #if best_result[1] < cost:
            #	cost = best_result[1]
            #	idx = best_result[0]
            #else:
            #	cost_history.extend([best_cost]*1000)
    print('done')
    return cost
