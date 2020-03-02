#!/usr/bin/python

# %%
import matplotlib
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "CMU Serif"
#matplotlib.rcParams['font.serif'] = "Computer Modern Roman"

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
plt.clf()

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
plt.clf()


# %% Example 3
#import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['presentation.mplstyle', 'grayscale'])
mat0 = genfromtxt("../data/Ex_03.txt")
#print(mat0)
lineA = plt.plot(mat0[0:-1:3,0], mat0[0:-1:3,1],"-k",fillstyle="none",markeredgewidth=2)
lineB = plt.plot(mat0[0:-1:3,0], mat0[0:-1:3,2],"ok",markeredgewidth = 2)
plt.legend(['Analytical','SISC++'],fontsize = 20)
plt.xlabel("$y$", fontsize=20)
plt.ylabel("$u(y)$", fontsize=20)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)

#plt.show()
plt.savefig('../pics/Ex_03.svg', bbox_inches='tight')
plt.clf()


# %% Example 12:
import numpy as np
import matplotlib.pyplot as plt
from pylab import genfromtxt;
plt.style.use(['presentation.mplstyle', 'grayscale'])
mat1 = genfromtxt("../data/Ex12.txt")
print(mat1)
lineA = plt.semilogy(mat1[:,0], mat1[:,1],"-ok",fillstyle="none",markeredgewidth=2,markersize=9)
lineB = plt.semilogy(mat1[:,0], mat1[:,2],"-^k",fillstyle="none",markeredgewidth = 2,markersize=9)
plt.legend(['$\sigma_0$','$\sigma_1$'],fontsize = 20)
plt.xlabel("$\omega$", fontsize=20)
#plt.ylabel("$u(y)$", fontsize=20)
plt.title('Singular values')
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
#plt.show()
plt.savefig('../pics/Ex_12.svg', bbox_inches='tight')
plt.clf()
# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

font = FontProperties()
#font.set_family('CMU Serif')
#font.set_name('Roman')
font.set_file('/Users/gokul/Library/Fonts/cmunbx.ttf')
#cmunbx
#cmunrm
#font.set_weight('Bold')
font.set_size(24)
#del matplotlib.font_manager.weight_dict['roman']
#matplotlib.font_manager._rebuild()
from pylab import genfromtxt;
plt.style.use(['presentation.mplstyle', 'grayscale'])
mat1 = genfromtxt("../data/Ex_12_2.txt")
print(mat1)
lineA = plt.semilogy(mat1[:,0], mat1[:,1],"-ok",fillstyle="none",markeredgewidth=2,markersize=9)
#plt.legend(['$\sigma_0$','$\sigma_1$'],fontsize = 20)
plt.xlabel("$\omega$", fontsize=20)
#plt.ylabel("$u(y)$", fontsize=20)
plt.title('Power spectral density',FontProperties=font)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.show()
#plt.savefig('../pics/Ex_12_2.svg', bbox_inches='tight')
#plt.clf()

# %%
