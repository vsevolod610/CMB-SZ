""""
Write file SZ_data in a right way.

require:    data.py
make:       SZ_data.txt

coments:    Pay attention to spaces in rows
"""


from data import nu_data, SZ_data, SZ_data_errors

with open("SZ_data.txt", 'w') as File:
    File.write('SZ data for freq. 70, 100, 143, 217, 353 GHz\n{}\n'.format(len(nu_data)))
    for j in range(len(nu_data)):
        s = '{: .2e}'.format(SZ_data[j]) + '  ' + '{:.2e}'.format(SZ_data_errors[0, j]) + '  ' + \
            '{:.2e}'.format(SZ_data_errors[1, j])
        s = s + '\n'
        File.write(s)