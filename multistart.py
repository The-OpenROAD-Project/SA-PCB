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

def f(idx):
    res = subprocess.check_output(["./sa", "-x", str(idx)])
    return res.strip()

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
            #best_result = max(results,key=itemgetter(1)) # max result by cost
            #if best_result[1] < cost:
            #	cost = best_result[1]
            #	idx = best_result[0]
            #else:
            #	cost_history.extend([best_cost]*1000)
    print('done')
    return cost

multistart()

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
            #best_result = max(results,key=itemgetter(1)) # max result by cost
            #if best_result[1] < cost:
            #	cost = best_result[1]
            #	idx = best_result[0]
            #else:
            #	cost_history.extend([best_cost]*1000)
    print('done')
    return cost
