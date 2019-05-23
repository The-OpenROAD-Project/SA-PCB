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
