import numpy as np
import numpy.random as rand
import ReadANDWrite as RW
from SZfunction import SZfunction
from trans_function import spec_trans, filtrs, spec_data


def getRAND(mean=0.0, sigma=1.0, mod='uniform', shape=1, end=100):
    min = mean - sigma
    max = mean + sigma
    if mod == 'uniform':
        return (max - min) * rand.rand(shape) + min
    elif mod == 'norm':
        k = 0
        while 1:
            k += 1
            RandList = sigma * rand.randn(shape) + mean
            if np.all(min <= RandList <= max):
                return RandList
            if k == end:
                print('The program cannot find value(s) in the range [{}, {}]\nreturned 1/0'.format(min, max))
                return 1/0



if __name__ == "__main__":
    flag = True
    while flag:
        charbool = input('Do you want to create 10 SZ_data files? (yes/no)\nif \'no\' appear, you will choose your own number\n')
        if charbool == 'yes':
            N = 10
            flag = False
        elif (charbool == 'no'):
            charN = input('How many SZ_data files do you want to create?\n')
            N = int(charN)
            flag = False
        else:
            print('Please, enter \'yes\' or \'no\' in lowercase without blanks')

    nuList = [70.0, 100.0, 143.0, 217.0, 353.0]
    relErrNuList = [0.22, 0.08, 0.1, 0.98, 0.30]
    nuStr = '70, 100, 143, 217, 353'

    for i in range(N):

        Name = 'SZdatas\\SZ_data{}.txt'.format(i + 1)
        print(Name)

        meanT0 = 2.7260
        sigmaT = 0.0013
        T0 = getRAND(mean=meanT0, sigma=sigmaT, mod='norm')[0]
        print(T0)

        meanTau = 1.025
        sigmaTau = 0.975
        Tau = getRAND(mean=meanTau, sigma=sigmaTau)[0]
        print(Tau)

        meanBeta = 0.0
        sigmaBeta = 0.003
        Beta = getRAND(mean=meanBeta, sigma=sigmaBeta)[0]
        #perhaps we have to use uniform destribution
        print(Beta)

        meankTe = 5.5
        sigmakTe = 4.5
        kTe = getRAND(mean=meankTe, sigma=sigmakTe)[0]
        #perhaps we have to use uniform destribution
        print(kTe)

        print()
        Str = ['SZ data for freq. {} GHz | T0 = {}, Tau = {:.2}, beta = {:.1e}, kTe = {:.2}'.format(nuStr, meanT0, Tau, Beta, kTe),
               str(len(nuList))]
        for nu0, wave, relErr in zip(nuList, filtrs, relErrNuList):
            nu = spec_trans[wave]['nu']
            tr = spec_trans[wave]['tr']
            SZf = SZfunction(T0=meanT0, kTe=kTe, beta=Beta, Tau=Tau, nu=nu)
            deltaTsz = np.trapz(SZf*tr, nu)
            relErrP = getRAND(mean=relErr, sigma=relErr/5.0, mod='norm')[0] #???????????????????
            relErrN = getRAND(mean=relErr, sigma=relErr/5.0, mod='norm')[0] #???????????????????
            Str.append('{: .2e} {: .2e} {: .2e}'.format(deltaTsz, abs(deltaTsz) * relErrP, abs(deltaTsz) * relErrN))
            print('{: .2e} {: .2e} {: .2e}'.format(deltaTsz, abs(deltaTsz)*relErrP, abs(deltaTsz)*relErrN))


        print()


        RW.WRITING(Name=Name, Strin=Str)