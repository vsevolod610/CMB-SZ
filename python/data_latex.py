import csv
import numpy as np

with open('data.csv', newline='') as File:
    reader = csv.reader(File, delimiter=',')
    rows = [row for row in reader]
    rows = np.array(rows[1:])
    methods = list(rows[:, 0])
    papers = list(rows[:, 1])

with open('data_latex.txt', 'w') as File:
    File.write(r'\[ \begin{array}{cllll}')
    File.write(r'\text{Metod} & \text{References}& \text{Cluster} & z & Tz, K \\')
    pmethod, ppaper = ['none', 'none']
    for method, paper, cluster, *data in rows:
        if pmethod != method:
            count = methods.count(method)
            method_str = r'\multirow{' + str(count) + r'}{*}{\text{' + str(method) + '}}'
            hhline = r'\hline\hline'
            t = True
        else:
            method_str = ' '
            hhline = ''
            t = False
        if ppaper != paper:
            count = papers.count(paper)
            paper_str = r'\multirow{' + str(count) + r'}{*}{\text{' + str(paper) + '}}'
            if t :
                line = ' '
            else:
                line = r'\cline{2 - 5}'
        else:
            paper_str = ' '
            line = ''
        s = hhline + line + method_str + ' & ' + paper_str + ' & ' + cluster + '&' 
        s += data[0] + '&' 
        s += data[1] + '^{+' + data[2] + '}_{-' + data[3] + '}'+ r'\\ ' + '\n'
        pmethod, ppaper = [method, paper]
        File.write(s)
    File.write(r'\hline\hline \end{array} \]')

