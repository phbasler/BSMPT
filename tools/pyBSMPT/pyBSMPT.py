import subprocess
import os
import pandas as pd

class BSMPT:

    def __init__(self, modelIN, input_fileIN, bsmpt_pathIN):
        self.bsmpt_path  = str(bsmpt_pathIN)
        self.model       = str(modelIN)
        self.input_file  = str(input_fileIN)

    '''
        Call BSMPT executable
    '''
    def calc_BSMPT(self, outfile, firstline, secondline, logfile=None):
        os.makedirs(os.path.abspath(os.path.dirname(outfile)), exist_ok=True)
        if logfile is None:
            logfile = outfile + '.log'
        self.info_print(outfile, firstline, secondline)
        with open(logfile, 'w') as f:
            subprocess.run([os.path.join(self.bsmpt_path, './bin/BSMPT'), self.model, self.input_file, outfile, str(firstline), str(secondline)], stdout=f, stderr=subprocess.STDOUT)

    def info_print(self, outfile, firstline, secondline):
        print(f'Called BSMPT at {self.bsmpt_path}.')
        print(f'Model: {self.model}')
        print(f'Inputfile: {self.input_file}')
        print(f'Output saved to {outfile}.')
        print(f'Firstline: {firstline}')
        print(f'Secondline: {secondline}')

    '''
        Analysis section
    '''
    def get_value(self, file, column):
        df = pd.read_csv(file, sep='\t')
        return df[column].values