"""
///////////////////////////////////////////////////////////////////////////////
// Authors: Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng, Steven Swanson).
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////
"""

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
