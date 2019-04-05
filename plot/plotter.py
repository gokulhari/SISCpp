#!/usr/bin/python
# %% Example 1
import numpy as np
import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['plot/presentation.mplstyle', 'grayscale'])

import matplotlib.pyplot as plt
plt.style.use(['plot/presentation.mplstyle', 'grayscale'])

mat0 = genfromtxt("data/Ex_01.txt");
plt.plot(mat0[:,0], mat0[:,1],"-k",fillstyle="none")
plt.xlabel("$y$", fontsize=20)
plt.ylabel("$u(y)$", fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig('pics/Ex_01.svg', bbox_inches='tight')

# %% Example 2
import numpy as np
import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['plot/presentation.mplstyle', 'grayscale'])

import matplotlib.pyplot as plt
plt.style.use(['plot/presentation.mplstyle', 'grayscale'])

mat0 = genfromtxt("data/Ex_02.txt");
plt.plot(mat0[:,0], mat0[:,1],"-k",fillstyle="none")
plt.xlabel("$y$", fontsize=20)
plt.ylabel("$u(y)$", fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig('pics/Ex_02.svg', bbox_inches='tight')
