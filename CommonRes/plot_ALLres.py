#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import emcee
import sys
from multiprocessing import Pool
from chainconsumer import ChainConsumer


res = []
with open('./ALLres2.dat') as f:
    for k, line in enumerate(f):
        if k==0:
            continue
        res.append(line.split())
        #print(res[-1])

print(res)


label = []
value = []
y_1 = []
y_1_bot = []
y_1_top = []
y_2 = []
y_2_bot = []
y_2_top = []
ind_1 = []
ind_2 = []

for i, r in enumerate(res):
    label.append(r[0])
    if r[2] != '-':
        y_1.append(r[2])
        y_1_bot.append(r[1])
        y_1_top.append(r[3])
        ind_1.append(i)
    if r[5] != '-':
        y_2.append(r[5])
        y_2_bot.append(r[4])
        y_2_top.append(r[6])
        ind_2.append(i)    

y_1 = np.array(y_1, dtype=float)
y_1_bot = np.array(y_1_bot, dtype=float)
y_1_top = np.array(y_1_top, dtype=float)
y_2 = np.array(y_2, dtype=float)
y_2_bot = np.array(y_2_bot, dtype=float)
y_2_top = np.array(y_2_top, dtype=float)

'''
res = np.array(res) #, dtype=float)

#print(res)

 = np.array(res[:,1:], dtype=float)

#print(value)



 = res[:, 0]

#print(label)
#z = res[:, 1]
y_1 = value[:, 1]

#print(y_1)
y_1_bot = value[:, 0]
y_1_top = value[:, 2]
y_2 = value[:, 4]
y_2_bot = value[:, 3]
y_2_top = value[:, 5]

'''
yerror_1 = [y_1 - y_1_bot, y_1_top - y_1]
yerror_2 = [y_2 - y_2_bot, y_2_top - y_2]

N = len(label)

fig, ax = plt.subplots(figsize=(14, 5))

ind = [i for i in range(N)]
#print(ind)

x = np.linspace(min(ind) - 0.01, max(ind) + 0.01, 100)

ax.errorbar(ind_1, y_1, yerr = yerror_1, capsize=3.5, mew=1.5, fmt='^r', alpha=0.5, label='by T0')        
ax.errorbar(ind_2, y_2, yerr = yerror_2, capsize=3.5, mew=1.5, fmt='xb', alpha=0.5, label='by Tz')

ax.set_xticks(ind)
ax.set_xticklabels(label, fontsize=15, rotation=0)#70) #'vertical'); # use LaTeX formatted labels

yticks = ax.get_yticks()
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.2f$" % y for y in yticks], fontsize=15); # use LaTeX formatted labels


labele = r'${: .4f}^{{{:+.4f}}}_{{{:+.4f}}}$'

ax.text(ind_1[0] + 0.01, y_1[0] + 0.0005, labele.format(y_1[0], yerror_1[1][0], -yerror_1[0][0]), horizontalalignment='left',
     verticalalignment='bottom', color='red')
     
ax.text(ind_2[0] + 0.01, y_2[0] + 0.0005, labele.format(y_2[0], yerror_2[1][0], -yerror_2[0][0]), horizontalalignment='left',
     verticalalignment='bottom', color='blue')
     
ax.text(ind_1[1] + 0.03, y_1[1] + 0.001, labele.format(y_1[1], yerror_1[1][1], -yerror_1[0][1]), horizontalalignment='left',
     verticalalignment='bottom', color='red')
     
ax.text(ind_1[2] - 0.03, y_1[2] - 0.001, labele.format(y_1[2], yerror_1[1][2], -yerror_1[0][2]), horizontalalignment='right',
     verticalalignment='top', color='red')


ax.fill_between(ind_1, y_1_top, y_1_bot, alpha=0.2, color='r', linewidth=1, linestyle='-')
ax.fill_between(ind_2, y_2_top, y_2_bot, alpha=0.2, color='b', linewidth=1, linestyle='-')
ax.plot(ind_1, y_1, '--r')
ax.plot(ind_2, y_2, '--b')

T0_plank = 2.7255
T0_plank_erorr = 0.0006

ax.plot(x, T0_plank + x * 0, '--k', label=r"$T_0$ COBE = $2.7255$K")
#ax.fill_between(x, x * 0 + T0_plank + T0_plank_erorr, x * 0 + T0_plank - T0_plank_erorr , alpha=0.2, color='k', linewidth=1, linestyle='-')


ax.set_title("Reaults", fontsize=18)
#fig.suptitle("Reaults", fontsize=20)


ax.legend(frameon=False, fontsize=15)

fig.tight_layout()
fig.savefig('fig_RES' + '.pdf', bbox_inches='tight')

