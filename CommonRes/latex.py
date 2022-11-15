import numpy as np

def read(path):
    data = []
    with open(path) as file:
        text = file.read().splitlines()
        print('Read: \n',*text, sep='\n')
        for line in text:
            data.append(line.split())
    return data

def write_latex(path, data):
    with open(path, 'w') as file:
        text = ''
        num_row = len(data)
        num_col = 2
        table_align =  'c' + num_col * 'l' 
        title = ['label', '$T_0 by 1$', '$T_0 by 2$']

        # text filling
        text += r'\begin{array}{%s}'
        text = text % table_align + '\n'
        text += r'\mbox{%s}' + num_col * r' & \mbox{%s}'
        text = text % tuple(title) + r'\\' + '\n'
        text += r'\hline' + '\n'
        for line in data[1:]:
            text += r'\mbox{%s}' + num_col * r' & {%s}^{+%s}_{-%s}'
            text = text % tuple(line) + r' \\' + '\n'
        text += r'\hline' + '\n'
        text += '\end{array}'

        # write in file
        print('Write: \n\n', text)
        file.write(text)

path_read = 'ALLres.dat'
path_write = 'ALLres.tex'
data = read(path_read)
write_latex(path_write, data)

