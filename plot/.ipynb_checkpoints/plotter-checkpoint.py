#!/usr/bin/python
# %% Example 1
import numpy as np
import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['presentation.mplstyle', 'grayscale'])

import matplotlib.pyplot as plt
plt.style.use(['presentation.mplstyle', 'grayscale'])

mat0 = genfromtxt("../data/Ex_01.txt");
plt.plot(mat0[:,0], mat0[:,1],"-k",fillstyle="none")
plt.xlabel("$y$", fontsize=20)
plt.ylabel("$u(y)$", fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig('../pics/Ex_01.svg', bbox_inches='tight')

# %% Example 2
import numpy as np
import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['presentation.mplstyle', 'grayscale'])

import matplotlib.pyplot as plt
plt.style.use(['presentation.mplstyle', 'grayscale'])

mat0 = genfromtxt("../data/Ex_02.txt");
plt.plot(mat0[:,0], mat0[:,1],"-k",fillstyle="none")
plt.xlabel("$y$", fontsize=20)
plt.ylabel("$u(y)$", fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.savefig('../pics/Ex_02.svg', bbox_inches='tight')


# %% Example 3
#import numpy as np

import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['presentation.mplstyle', 'grayscale'])
mat0 = genfromtxt("../data/Ex_03.txt")
lineA = plt.plot(mat0[0:2:-1,0], mat0[0:2:-1,1],"-sk",fillstyle="none")
lineB = plt.plot(mat0[0:2:-1,0], mat0[0:2:-1,2],"-xk")
plt.legend(['Analytical','SISC++'],fontsize = 20)
plt.xlabel("$y$", fontsize=20)
plt.ylabel("$u(y)$", fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)

#plt.show()
plt.savefig('../pics/Ex_03.svg', bbox_inches='tight')

# %% Ex 4

import matplotlib.pyplot as plt
line_up, = plt.plot([1,2,3], label='Line 2')
line_down, = plt.plot([3,2,1], label='Line 1')
plt.legend([line_up, line_down])


# %%
import matplotlib
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.serif'] = "Computer Modern Roman"
# %%
