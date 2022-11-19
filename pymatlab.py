import pyperclip as pc

text = ['120', '10', "0", '0.5', '0.5', '0.25', '1.5', '25',
     '4', "0", '0.5', '0.5', '0.25', '1.5', '0.08', '0.000062','0.003', '0.0017', '1.6', '0.05', '0.5', '0.52', '1', '1', '1', '1',
     '0.075', '0.18', '0.002', '0.000052', '3.8', '0.01', '0.0017', '0.05', '0.25', '20.5', '0.006']

variables = ['v1l', 'v2l', 'tdl', 'tfl', 'trl', 'pwl', 'perl', 'v1r', 'v2r', 'tdr',
                            'tfr', 'trr', 'pwr', 'perr', 'C_sas', 'L_sas', 'R_sas', 'L_sat',
                            'C_sat', 'R_sat', 'R_sar', 'R_scp', 'R_brain', 'R_liver',
                            'R_spleen', 'R_kidney', 'R_svn', 'C_pas',
                            'R_pas', 'L_pas', 'C_pat', 'R_pat', 'L_pat', 'R_par', 'R_pcp',
                            'C_pvn', 'R_pvn']

texts = ''

for i in range(len(variables)):
    texts += f"{variables[i]} = '{text[i]}';\n"

pc.copy(texts)
print(texts)