#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import emcee
import sys
from multiprocessing import Pool
from chainconsumer import ChainConsumer


res = []
with open('./ALLres.dat') as f:
    for k, line in enumerate(f):
        if k==0:
            continue
        res.append(line.split())
        #print(res[-1])


res = np.array(res) #, dtype=float)

#print(res)

value = np.array(res[:,1:], dtype=float)

#print(value)



label = res[:, 0]

#print(label)
#z = res[:, 1]
y_1 = value[:, 1]

#print(y_1)
y_1_bot = value[:, 0]
y_1_top = value[:, 2]
y_2 = value[:, 4]
y_2_bot = value[:, 3]
y_2_top = value[:, 5]


yerror_1 = [y_1 - y_1_bot, y_1_top - y_1]
yerror_2 = [y_2 - y_2_bot, y_2_top - y_2]

N = len(label)

fig, ax = plt.subplots(figsize=(14, 10))

ind = [i for i in range(N)]
#print(ind)

x = np.linspace(min(ind) - 0.01, max(ind) + 0.01, 100)

ax.errorbar(ind, y_1, yerr = yerror_1, capsize=3.5, mew=1.5, fmt='^r', alpha=0.5, label='by T0')        
ax.errorbar(ind, y_2, yerr = yerror_2, capsize=3.5, mew=1.5, fmt='xb', alpha=0.5, label='by Tz')

ax.set_xticks(ind)
ax.set_xticklabels(label, fontsize=15, rotation=70) #'vertical'); # use LaTeX formatted labels

yticks = ax.get_yticks()
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.2f$" % y for y in yticks], fontsize=15); # use LaTeX formatted labels


ax.fill_between(ind, y_1_top, y_1_bot, alpha=0.2, color='r', linewidth=1, linestyle='-')
ax.fill_between(ind, y_2_top, y_2_bot, alpha=0.2, color='b', linewidth=1, linestyle='-')
ax.plot(ind, y_1, '--r')
ax.plot(ind, y_2, '--b')

T0_plank = 2.7255
T0_plank_erorr = 0.0006

ax.plot(x, T0_plank + x * 0, '--k', label=r"$T_0$ COBE = $2.7255\pm0.0006$")
ax.fill_between(x, x * 0 + T0_plank + T0_plank_erorr, x * 0 + T0_plank - T0_plank_erorr , alpha=0.2, color='k', linewidth=1, linestyle='-')


#ax.set_title("Reaults", fontsize=18)
fig.suptitle("Reaults", fontsize=20)


ax.legend(frameon=False, fontsize=15)

fig.tight_layout()
fig.savefig('fig_RES' + '.pdf', bbox_inches='tight')

