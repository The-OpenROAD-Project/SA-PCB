"""
///////////////////////////////////////////////////////////////////////////////
// Authors: Ilgweon Kang and Lutong Wang
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng),
//          based on Dr. Jingwei Lu with ePlace and ePlace-MS
//
//          Many subsequent improvements were made by Mingyu Woo
//          leading up to the initial release.
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

import matplotlib
from matplotlib import animation

import matplotlib.pyplot as plt
import numpy as np
"""
wl = open('wl.txt').read().splitlines()
oa = open('oa.txt').read().splitlines()
cost = open('cost.txt').read().splitlines()
ar = open('accept_ratio.txt').read().splitlines()

wl = [float(w) for w in wl]
oa = [float(o) for o in oa]
cost = [float(c) for c in oa]
ar = [float(a) for a in ar]
"""
#plt_dict = {'wirelength':wl, 'overlap':oa, 'cost':cost, 'acceptance ratio':ar}

#plt_dict = {'acceptance_ratio':ar}

#for subject, sig in plt_dict.items():
#    plt.plot(sig, label=subject)
#plt.legend()
#plt.xlim((-1,1000000))


#plt.show()

import os
import pandas as pd
directory = os.fsencode('./cache_rudy/')
from tqdm import tqdm
routabilities = []
fnames = []

for (dirpath, dirnames, filenames) in os.walk('./cache_rudy'):
    fnames.extend(filenames)
fnames = sorted(fnames, key=lambda x: float(x.split('.')[0]))


for i,ff in enumerate(tqdm(fnames)):
     #mat = np.loadtxt('./cache_rudy/'+str(filename),delimiter='\t')
     mat = pd.read_table('./cache_rudy/'+str(ff),sep='\t',header=None).values[:,:-1]
     routabilities.append(mat)
import seaborn as sns

nx, ny = routabilities[0].shape

fig = plt.figure()
data = np.random.rand(nx, ny)
sns.heatmap(data, vmax=.8, square=True)

def init():
      sns.heatmap(np.zeros((nx, ny)), vmax=.8, square=True,cbar=False)

def animate(i):
    data = routabilities[i]
    ax = sns.heatmap(data, vmax=.8, square=True,cbar=False)
    ax.invert_yaxis()


plt.yticks = []
plt.xticks = []

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=None, repeat = True)
plt.title('routability')
#plt.gca().invert_yaxis()
plt.show()
