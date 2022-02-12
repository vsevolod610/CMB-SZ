# -*- coding: utf-8 -*-
""""
set clusters --> set of T0_i temperature

require:    mcmc_analyze.py
            ../OUT/MCchain{1-77}.dat, ../OUT/chainconsum{1-77}.dat
make:       ../Result/MCcain_analyze{1-77}.png, ../Result/chainconsum_analyze{1-77}.png
            ../Result/Result.txt

coments:
"""


from mcmc_analyze import MCcain_analyze, chainconsum_analize


for n in range(1, 77 + 1):
    filename = '../OUT/MCchain{}.dat'.format(n)
    picname = '../Result/MCcain_analyze{}.png'.format(n)
    MCcain_analyze(filename, picname=picname)
    print('MCcain : ', n)

with open("../Result/Result.txt", 'w') as File:
    for n in range(1, 77 + 1):
        filename = '../OUT/chainconsum{}.dat'.format(n)
        picname = '../Result/chainconsum_analyze{}.png'.format(n)
        estim = chainconsum_analize(filename, picname=picname)
        estim = estim.get_summary()[r'$T_0$']
        if None in estim:
            File.write(str(n) + '  None' + '\n')
            continue
        line = str(n) + ' {0: .4} {1: .4} {2: .4}'
        line = line.format(estim[1], estim[2] - estim[1], estim[1] - estim[0])
        line =  line + '\n'
        File.write(line)
        print('chain : ', line)
