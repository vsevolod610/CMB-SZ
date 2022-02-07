import os

def READING(Name, ListOfCol=[-1]):

    if os.path.exists(Name):
        with open(Name) as inf:
            Strin = []
            i=0
            for line in inf:
                sttr = ''
                if ListOfCol!=[-1]:
                    for j in ListOfCol:
                        sttr += line.split()[j] + ' '
                else:
                    NCol=len(line.split())
                    for j in range(NCol):
                        sttr += line.split()[j] + ' '

                Strin += [sttr]
    else:
        Strin = None

    return (Strin)


'''
LeastOfCel = (0, 1, 3)
Name = '/home/ilya/Desktop/SSW_NIR/Main_Program/Rawlins_ks18_tmp2.ovr_last'
Strin = READING(Name, LeastOfCel)
NumOfStr = len(Strin)
print(NumOfStr)
j = 0
Name = '/home/ilya/Desktop/SSW_NIR/Main_Program/FIRST.txt'
WRITING(Name, Strin)

'''

def WRITING(Name, Strin=[]):
    NumOfStr = len(Strin)
    with open(Name, 'w') as ouf:
        for j in range(NumOfStr):
            ouf.write(Strin[j] + '\n')
