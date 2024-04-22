#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import numpy as np
import sys

def convert(InputFile, OutputFile):
    print(f'Reading {InputFile}.')
    print(f'Output is saved to {OutputFile}.')

    HasIndexCol=False
    Separator='\t'
    Type='yuktype'
    Lambda1='L1'
    Lambda2='L2'
    Lambda3='L3'
    Lambda4='L4'
    re_Lambda5='re_L5'
    im_Lambda5='im_L5'
    tanbeta='tbeta'
    re_m12squared='re_m12sq'

    with open(InputFile, 'r') as file:
        df = pd.read_csv(file,index_col=HasIndexCol,sep=Separator)

    frontcol=[Type,Lambda1,Lambda2,Lambda3,Lambda4,re_Lambda5,im_Lambda5,re_m12squared,tanbeta]

    Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]

    df = df[Col]

    with open(OutputFile, 'w') as file:
        df.to_csv(file, index=False, sep='\t')

if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
