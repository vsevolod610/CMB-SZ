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

print(res)


label = []
value = []
y_1 = []
y_1_bot = []
y_1_top = []
y_2 = []
y_2_bot = []
y_2_top = []
y_3 = []
y_3_bot = []
y_3_top = []
ind_1 = []
ind_2 = []
ind_3 = []


for i, r in enumerate(res):
    label.append('('+r[0]+')')
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

print(y_1)
yerror_1 = [y_1 - y_1_bot, y_1_top - y_1]
yerror_2 = [y_2[:3] - y_2_bot[:3], y_2_top[:3] - y_2[:3]]
yerror_3 = [y_2[3:] - y_2_bot[3:], y_2_top[3:] - y_2[3:]]



N = len(label)

fig, ax = plt.subplots(figsize=(12, 4))

ind = [i for i in range(N)]
#print(ind)

ax.set_xlim(-0.5, N - 0.5)

x = np.linspace(min(ind) - 0.01, max(ind) + 0.01, 100)

ax.errorbar(ind_1, y_1, yerr = yerror_1, capsize=5, mew=1.5, fmt='^r', alpha=0.5, markersize=8, label='по T0')        
ax.errorbar(ind_2[:3], y_2[:3], yerr = yerror_2, capsize=5, mew=1.5, fmt='xb', alpha=0.5, markersize=8, label='по Tz')

ax.errorbar(ind_2[3:], y_2[3:], yerr = yerror_3, capsize=5, mew=1.5, fmt='*g', alpha=0.5, markersize=8, label='Доп. системы.')

ax.set_xticks(ind)
llabele = ["({:})".format(i+1) for i in ind]
ax.set_xticklabels(label, fontsize=15, rotation=+0)#70) #'vertical'); # use LaTeX formatted labels
#ax.set_xticklabels(llabele, fontsize=16, rotation=+0)#70) #'vertical');

yticks = ax.get_yticks()
ax.set_yticks(yticks[1::2])
ax.set_yticklabels(["$%.2f$" % y for y in yticks[1::2]], fontsize=16); # use LaTeX formatted labels
ax.set_ylim([2.70,2.78])

labele = r'${: .4f}^{{{:+.4f}}}_{{{:+.4f}}}$'


fonts = 15
if 0:
    ax.text(ind_1[0] + 0.01, y_1[0] + 0.0005, labele.format(y_1[0], yerror_1[1][0], -yerror_1[0][0]), horizontalalignment='left',
         verticalalignment='bottom', color='red')
         
    ax.text(ind_2[0] + 0.01, y_2[0] + 0.0005, labele.format(y_2[0], yerror_2[1][0], -yerror_2[0][0]), horizontalalignment='left',
         verticalalignment='bottom', color='blue')
         
    ax.text(ind_1[1] + 0.03, y_1[1] + 0.001, labele.format(y_1[1], yerror_1[1][1], -yerror_1[0][1]), horizontalalignment='left',
         verticalalignment='bottom', color='red')
         
    ax.text(ind_1[2] - 0.03, y_1[2] - 0.001, labele.format(y_1[2], yerror_1[1][2], -yerror_1[0][2]), horizontalalignment='right',
         verticalalignment='top', color='red')
elif 1:
    for ind, y, perr, merr in zip(ind_1, y_1, yerror_1[1], -yerror_1[0]): 
        ax.text(ind + 0.01, y + 0.001, labele.format(y, perr, merr), horizontalalignment='center',
         verticalalignment='bottom', color='red', fontsize=fonts)
    for ind, y, perr, merr in zip(ind_2[:3], y_2[:3], yerror_2[1][:3], -yerror_2[0][:3]): 
        ax.text(ind + 0.01, y + 0.002, labele.format(y, perr, merr), horizontalalignment='center',
         verticalalignment='bottom', color='blue', fontsize=fonts)
    for ind, y, perr, merr in zip(ind_2[3:], y_2[3:], yerror_3[1], -yerror_3[0]): 
        ax.text(ind + 0.01, y + 0.002, labele.format(y, perr, merr), horizontalalignment='center',
         verticalalignment='bottom', color='green', fontsize=fonts)


if 1:
    ax.fill_between(ind_1, y_1_top, y_1_bot, alpha=0.2, color='r', linewidth=1, linestyle='-')
    ax.fill_between(ind_2[:3], y_2_top[:3], y_2_bot[:3], alpha=0.2, color='b', linewidth=1, linestyle='-')
    ax.fill_between(ind_2[3:], y_2_top[3:], y_2_bot[3:], alpha=0.2, color='g', linewidth=1, linestyle='-')
    ax.plot(ind_1, y_1, '--r',alpha=0.4)
    ax.plot(ind_2[:3], y_2[:3], '--b',alpha=0.4)
    ax.plot(ind_2[3:], y_2[3:], '--g',alpha=0.4)


T0_plank = 2.7255
T0_plank_erorr = 0.0006

ax.plot(x, T0_plank + x * 0, '--k', label=r"$T_0$ = $2.7255$K")
#ax.fill_between(x, x * 0 + T0_plank + T0_plank_erorr, x * 0 + T0_plank - T0_plank_erorr , alpha=0.2, color='k', linewidth=1, linestyle='-')


ax.set_title(r"Результаты по $T_0$", fontsize=18)
ax.set_ylabel(r"$T_0$, K", fontsize=18)#, rotation='horizontal')
#fig.suptitle("Reaults", fontsize=20)


ax.legend(frameon=False, fontsize=16)

fig.tight_layout()
fig.savefig('fig_RES_RUS' + '.pdf', bbox_inches='tight')

